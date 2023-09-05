#ifndef ASYNC_DIST_KMER_COUNTER_HPP
#define ASYNC_DIST_KMER_COUNTER_HPP

#include <stdexcept>
#include <pthread.h>

#include <mpi.h>
#include <stop_clock.hpp>



namespace intel {
namespace esc {

  template<typename T>
  struct kmer_pair_cmp_default_t {
    inline bool operator() (const T& a, const T& b) const {
      return (a.seq < b.seq);
    }
  }; // struct kmer_pair_cmp_default_t //

  template<typename kmer_t, typename kmer_pair_t>
  class Async_Distributed_Kmer_Counter {

    public:
      //////////////////////////////////////////////////////////////////////////
      typedef Async_Distributed_Kmer_Counter dist_kmer_counter_t;
      //////////////////////////////////////////////////////////////////////////

      Async_Distributed_Kmer_Counter(size_t processor_count=1, size_t rank=0,
          size_t kmer_len=32UL, size_t max_kmers=100000000UL,
          char delimiter='\n', char skip_symbol='>') // 100mil //
        : processor_count_(processor_count), rank_(rank), kmer_len_(kmer_len),
          max_kmers_(max_kmers),
          delimiter_(delimiter), skip_symbol_(skip_symbol) , reads_beg_(NULL),
          reads_end_(NULL), comm_thread_(), comp_thread_(),
          is_send_buffer_full_(), all_kmers_processed_(), send_counts_kmer_(), 
          kmer_send_buffer_(),
          kmer_freq_map_(), kmer_mask_() 
      {

        kmer_mask_ = ((~0UL)) >> ((sizeof(kmer_t)*8) - (2*kmer_len));
      }

      template<typename kmer_pair_container_t>
      int run(const char *reads_beg, const char *reads_end, 
          kmer_pair_container_t& final_kmer_container) {
        reads_beg_ = reads_beg;
        reads_end_ = reads_end;

        all_kmers_processed_.store(false);
        is_send_buffer_full_[0].store(false);
        is_send_buffer_full_[1].store(false);

#if 0
        if (pthread_create(&comm_thread_, NULL,
                communication_thread_entry, this)) {
          throw std::runtime_error("[pthread_create failed]: comm_thread_");
        }
#endif


        if (pthread_create(&comp_thread_, NULL, compute_thread_entry, this)) {
          throw std::runtime_error("[pthread_create failed]: compute_thread_");
        }
        this->communication_thread_run();
        pthread_join(comp_thread_, NULL);
        
        //pthread_join(comm_thread_, NULL);


        final_kmer_container.resize(kmer_freq_map_.size());
        size_t i=0; 
        // populate the final_kmer_container //
        for (auto itr=kmer_freq_map_.begin(); itr!=kmer_freq_map_.end(); 
              ++itr) {
          final_kmer_container[i++] = kmer_pair_t{itr->first, (int)itr->second};
        }
        return num_batch_transfers_;
      }

      
    private:

      //NOTE: this routine is taken from PakMan code //
      template <typename T> 
      inline long uhash31( uint64_t a, uint64_t b, T x)
      {
#define MOD 2147483647
#define HL 31

        T result;
        long lresult;  

        result=(a * x) + b;
        result = ((result >> HL) + result) & MOD;
        lresult=(long) result; 
        
        return(lresult);
      }

      // NOTE: this routine is taken from PaKMan code //
      int retrieve_proc_id(kmer_t kmer) {
        static const uint64_t hasha = 68111;
        static const uint64_t hashb = 105929;
        int proc_num = 
          (int)((int)(uhash31(hasha, hashb, kmer)) % (int)processor_count_);
        return proc_num;

      }
      //////////////////////////////////////////////////////////////////////////
      //         CONCURRENT COMMUNICATION AND COMPUTE THREADS                 // 

      static void *communication_thread_entry(void *ptr) {
        Async_Distributed_Kmer_Counter *kmer_counter_ptr = 
            (Async_Distributed_Kmer_Counter *)(ptr);
        if (!kmer_counter_ptr) {
          throw std::logic_error("Invalid cast of the raw pointer ");
        }

        kmer_counter_ptr->communication_thread_run();
        return NULL;
      }

      void* communication_thread_run() {
        // STEP-0: scan and find out the total number of kmers and divide that 
        // by the max_kmer_count_ to find the rounds //
        size_t my_rounds = 1UL;
        const char *reads_beg = reads_beg_;
        const char *reads_end = reads_end_;
        const char *read_start = NULL;
        const char *read_end = NULL;

        size_t local_kmers = 0UL, read_len;
        size_t read_count = 0;
        {
          //intel::esc::stop_clock_t clock("[Scan K-mer Time]:");
          while (reads_beg < reads_end) {
            reads_beg = skip_till_symbol(reads_beg, reads_end, skip_symbol_);
            ++reads_beg;
            if (reads_beg >= reads_end) { break; }

            reads_beg = skip_till_symbol(reads_beg, reads_end, delimiter_);
            ++reads_beg;
            if ((reads_beg+kmer_len_) >= reads_end) { break; }

            ++read_count;

            read_end = skip_till_symbol(reads_beg, reads_end, delimiter_);
            read_len = (size_t)(read_end - reads_beg);
            if (read_len >= kmer_len_) {
              local_kmers += ((read_len - kmer_len_) + 1UL);
            }
            ++read_end;
            reads_beg = read_end;
          }
        }

        //printf("[processor-%lu local_kmers=%lu read_count=%lu]\n", rank_, 
         //   local_kmers, read_count);

        size_t global_max_kmers; // max k-mers across all processors //
        {
          //intel::esc::stop_clock_t clock("[MPI_Allreduce time:]:");
          int ret = MPI_Allreduce(&local_kmers, &global_max_kmers, 1,
              MPI_UNSIGNED_LONG, MPI_MAX, MPI_COMM_WORLD);
          if (ret != MPI_SUCCESS) {
            throw std::runtime_error(
                "[MPI_Allreduce failed with code="+std::to_string(ret));
          }
        }

        size_t rem = global_max_kmers%max_kmers_;
        size_t rounds = (global_max_kmers/max_kmers_) + (rem ? 1UL : 0UL);
        rounds = std::max(1UL, rounds);


        num_batch_transfers_ = (int) rounds;
        //printf("[processor-%lu rounds=%lu]\n", rank_, rounds);

        std::vector<int> recv_counts_kmer(processor_count_, 0);
        std::vector<int> recv_displacements_kmer(processor_count_, 0);
        std::vector<int> send_displacements_kmer(processor_count_, 0);

        // create MPI_Type (credit: taken from PakMan transfer_kmers) //
        MPI_Datatype mpi_kmer_pair_type;
        MPI_Type_contiguous(sizeof(kmer_pair_t), MPI_BYTE, &mpi_kmer_pair_type);
        MPI_Type_commit(&mpi_kmer_pair_type);

        std::vector<kmer_pair_t> kmer_recv_buffer;
        int ret;
        size_t total_recv_count; 

        for (size_t i=0; i<rounds; i++) {
          std::atomic<bool> & is_send_buffer_full = is_send_buffer_full_[i%2UL];

          if (!all_kmers_processed_) {
            // Invariant: to start communication of i^th round the corresponding
            // buffer must be full //
            do { } while(!is_send_buffer_full);
          }
          //printf("[Communication]-%lu starting send buffer round=%lu\n",
           //   rank_, i);
          //fflush(stdout);

          // STEP-1.0: compute the send count displacements //
          std::vector<int>& send_counts_kmer = send_counts_kmer_[i%2UL];
          send_displacements_kmer[0] = 0;
          for (size_t i=1; i<processor_count_; i++) {
            send_displacements_kmer[i] = 
                send_counts_kmer[i-1] + send_displacements_kmer[i-1];
          }

          MPI_Barrier(MPI_COMM_WORLD);

          // STEP-1.1: send and receive partition counts //
          ret = MPI_Alltoall(send_counts_kmer.data(), 1, MPI_INT, 
              recv_counts_kmer.data(), 1, MPI_INT, MPI_COMM_WORLD);
          if (ret != MPI_SUCCESS) {
            throw std::runtime_error("[MPI_Alltoall failed with code="+
                  std::to_string(ret)+"]");
          }

          // STEP-1.2: allocate array for kmers this processor needs to 
          // receive from all kmers //
          recv_displacements_kmer[0] = 0;
          for (size_t i=1; i<processor_count_; i++) {
            recv_displacements_kmer[i] = 
                recv_counts_kmer[i-1] + recv_displacements_kmer[i-1];
          }
          total_recv_count = recv_displacements_kmer[processor_count_-1UL] +
                recv_counts_kmer[processor_count_-1UL];
          //printf("[RECV-COUNT]-%lu = %lu\n", rank_, total_recv_count);
          //fflush(stdout);
          kmer_recv_buffer.resize(total_recv_count);

          std::vector<kmer_pair_t>& kmer_send_buffer = kmer_send_buffer_[i%2UL];

          ret = MPI_Alltoallv(
              // send arguments //
              kmer_send_buffer.data(), send_counts_kmer.data(), 
                send_displacements_kmer.data(), mpi_kmer_pair_type,
              // recv arguments //
              kmer_recv_buffer.data(), recv_counts_kmer.data(), 
                recv_displacements_kmer.data(), mpi_kmer_pair_type,
              MPI_COMM_WORLD
              );

          if (ret != MPI_SUCCESS) {
            throw std::runtime_error("[MPI_Alltoallv failed with code="+
                  std::to_string(ret)+"]");
          }

          send_counts_kmer.resize(processor_count_, 0);
          is_send_buffer_full = false;
          //printf("[Comm]-%lu setting send buffer full = false\n", rank_);
          //fflush(stdout);

          for (auto itr=kmer_recv_buffer.begin(); itr!=kmer_recv_buffer.end();
                ++itr) {
            kmer_freq_map_[itr->seq] += itr->k_count;
          }
          //printf("[Accumulate]-%lu kmer_count=%lu\n", 
           //     rank_, kmer_freq_map_.size());
          //fflush(stdout);
        }

        MPI_Type_free(&mpi_kmer_pair_type);
        return NULL;
      }


      static void *compute_thread_entry(void *ptr) {
        Async_Distributed_Kmer_Counter *kmer_counter_ptr = 
            (Async_Distributed_Kmer_Counter *)(ptr);
        if (!kmer_counter_ptr) {
          throw std::logic_error("Invalid cast of the raw pointer ");
        }

        kmer_counter_ptr->compute_thread_run();
        return NULL;
      }

      struct thread_state_t {
        const char *beg_;
        const char *end_;
        bool inside_read_;
        kmer_t kmer_;
        thread_state_t()
          : beg_(NULL), end_(NULL), inside_read_(false), kmer_() {}

        void print(size_t rank, size_t tid) const {
          printf("[p-%lu tid=%lu beg=%lx end=%lx ]\n", rank, tid, beg_, end_); 
        }

      }; // thread_state_t //

      struct KmerPairCompareHashed {
        inline bool operator() (const kmer_pair_t& a ,
              const kmer_pair_t& b) const {

          return (a.k_count == b.k_count) ? 
              (a.seq < b.seq) : (a.k_count < b.k_count);
        }
      }; // struct KmerPairCompareHashed //

      inline uint8_t char_to_el(char ch) {
        return (uint8_t)((((uint8_t)ch)>>1) & 0x7);
      }

      inline kmer_t kmer_shift(kmer_t kmer_in, uint8_t el) {
        return (kmer_t)((kmer_in<<2) | (kmer_t)el) & (kmer_t)kmer_mask_;
      }

      void* compute_thread_run() {

        size_t tcount = 1UL;
        size_t data_size = ((size_t) reads_end_) - ((size_t) reads_beg_);

#pragma omp parallel
        {
          if (!omp_get_thread_num())
            tcount = omp_get_num_threads();
        }

        std::vector<thread_state_t> thread_state(tcount);

        // STEP-0: initialize the thread states // 
        {
#pragma omp parallel
          {
            size_t my_id = omp_get_thread_num(); 
            thread_state_t& my_state = thread_state[my_id];
            size_t chunk_size = data_size/tcount;

            const char *symbol_begin = reads_beg_ + (my_id*chunk_size);
            const char *symbol_end = symbol_begin + chunk_size;

            if (symbol_begin >= reads_end_) { 
              symbol_begin = reads_end_; 
              symbol_end = reads_end_;
            }
            if (symbol_end >= reads_end_) {
              symbol_end = reads_end_;
            }

            // adjust symbol_end //
            symbol_end = skip_till_symbol(symbol_end, reads_end_, skip_symbol_);
            ++symbol_end;

            if (symbol_end >= reads_end_) {
              symbol_end = reads_end_;
            }

            my_state.beg_ = symbol_begin;
            my_state.end_ = symbol_end;
            my_state.inside_read_ = false;
            my_state.kmer_ = kmer_t(0);
          }
        }

        std::vector<kmer_pair_t> kmer_pair_buffer; 
        std::atomic<size_t> next_avail_idx={0UL};

        kmer_pair_buffer.resize(max_kmers_);

        size_t round_id = 0UL;
        do {
          // STEP-0: reset the buffer index //
          next_avail_idx.store(0UL);

          // STEP-1: start filling kmers into the kmer_pair_buffer in parallel//
#pragma omp parallel
          {
            size_t my_id = omp_get_thread_num(); 
            thread_state_t& my_state = thread_state[my_id];
            size_t next_idx;

            const char *beg = my_state.beg_;
            const char *end = my_state.end_;
            kmer_t kmer = my_state.kmer_; // default is 0 // 

            if (my_state.inside_read_) {
              my_state.inside_read_ = false;
              goto RESTORE_THREAD_STATE; 
            }

            while (beg < end) {
              for (; (beg<end) && (*beg != '>'); ++beg) { }
              for(; (beg<end) && (*beg !='\n'); beg++) {/*noop*/}
              if ((beg+kmer_len_) >= end) { break; }

              ++beg;

              kmer=0;
              for (int i=0; (*beg != '\n') && (i<kmer_len_-1); ++i, ++beg) {
                kmer = kmer_shift(kmer, char_to_el(*beg));
              }

              while ((beg < end) && (*beg != '\n')) {
                kmer = kmer_shift(kmer, char_to_el(*beg));

RESTORE_THREAD_STATE:
                next_idx = 
                    next_avail_idx.fetch_add(1UL, std::memory_order_relaxed);
                if (next_idx >= max_kmers_) {
                  // save state cannot fill the buffer try in next round. //
                  my_state.inside_read_ = true;
                  goto SAVE_THREAD_STATE;
                }
                // put the kmer into kmer pair buffer //
                kmer_pair_buffer[next_idx] = 
                    kmer_pair_t{kmer, retrieve_proc_id(kmer)};
                ++beg;
              }
              ++beg;
            }

SAVE_THREAD_STATE:
            my_state.beg_ = beg;
            my_state.end_ = end;
            my_state.kmer_ = kmer;
          }


          // STEP-2: now sort and remove duplicates in parallel from the
          // kmer_pair_buffer //
          {

            size_t size = processor_count_;
            size_t num_kmers = std::min((size_t) next_avail_idx, max_kmers_);
            KmerPairCompareHashed kmer_pair_compare;
            oneapi::tbb::parallel_sort(kmer_pair_buffer.begin(), 
                  kmer_pair_buffer.begin()+num_kmers, kmer_pair_compare);
            
            std::atomic<int> partition_counts[size] = { {0} };
            std::vector<size_t> uniq_kmer_counts(tcount, 0UL);
            std::vector<size_t> kmer_starts(tcount, 0UL);
            

            // STEP-2.1: parallel duplicate removal, frequency count 
            // and partition count 
#pragma omp parallel
            {
              size_t my_id = omp_get_thread_num(); 
              size_t data_size = num_kmers; 
              size_t chunk_size = data_size/tcount;

              size_t start_idx = my_id*chunk_size;
              size_t end_idx = (my_id == tcount-1UL) ? data_size :
                  start_idx + chunk_size;

              if (start_idx >= data_size) {
                start_idx = data_size; end_idx = data_size;
              }

              if (end_idx >= data_size) {
                end_idx = data_size;
              }

              if (start_idx > 0){
                kmer_t prev_kmer = kmer_pair_buffer[start_idx-1UL].seq;
                while ( (start_idx < end_idx) && 
                        (prev_kmer == kmer_pair_buffer[start_idx].seq )) {
                  ++start_idx;
                }
              }
              kmer_starts[my_id] = start_idx;
              int curr_freq = 1, kmer_part_id;
              size_t uniq_idx = start_idx;

              for (size_t i=start_idx+1; i<end_idx; i++) {
                if (kmer_pair_buffer[uniq_idx].seq == kmer_pair_buffer[i].seq) {
                  curr_freq++;
                } else {
                  // update the partition counts //
                  kmer_part_id = kmer_pair_buffer[uniq_idx].k_count;
                  if ((kmer_part_id < 0) || (kmer_part_id > processor_count_)) {
                    throw std::runtime_error("kmer_part_id="+
                          std::to_string(kmer_part_id)+"  processor_count_="+
                          std::to_string(processor_count_));
                  }
                  partition_counts[kmer_part_id].fetch_add(1,
                        std::memory_order_relaxed);
                  // update the frequency //
                  kmer_pair_buffer[uniq_idx].k_count = curr_freq;

                  // copy the new unique kmer into next slot //
                  kmer_pair_buffer[++uniq_idx] = kmer_pair_buffer[i];
                  curr_freq = 1;
                }
              }

              // the k-mer at the boundary may have partial frequency //
              for (size_t i=end_idx; (i<data_size) && 
                    (kmer_pair_buffer[uniq_idx].seq == kmer_pair_buffer[i].seq); 
                      i++) 
              {
                curr_freq++; 
              }

              // trailing case //
              {
                kmer_part_id = kmer_pair_buffer[uniq_idx].k_count;
                partition_counts[kmer_part_id].fetch_add(1,
                    std::memory_order_relaxed);
                // update the frequency //
                kmer_pair_buffer[uniq_idx].k_count = curr_freq;
              }
              uniq_kmer_counts[ my_id ] = (uniq_idx - start_idx) + 1UL;
            }



            // TODO(vamsikku): try to see if this can be done in parallel
            std::vector<size_t> prefix_sum(tcount, 0UL);
            size_t cum = uniq_kmer_counts[0UL]; 
            prefix_sum[0UL] = 0UL;
            for (size_t i=1; i<tcount; i++) {
              prefix_sum[i] = cum;
              cum += uniq_kmer_counts[i];
            }


            // STEP-2.3: WAIT until the corresponding send buffer is full //
            // the communication thread will set is_send_buffer_full=false after
            // its done doing an all2allv //
            // TODO(vamsikku): may be we need to replace this with mutex to 
            // avoid starvation.
            std::atomic<bool> &is_send_buffer_full =
                is_send_buffer_full_[round_id%2UL]; 
            do { } while(is_send_buffer_full);


            //STEP-2.2: fill the partition counts //
            std::vector<int> &scounts_kmer = send_counts_kmer_[round_id%2UL];
            scounts_kmer.resize(size, 0);
#pragma omp parallel for
            for (size_t i=0; i<size; i++) {
              scounts_kmer[i] = int(partition_counts[i]);
            }

            // STEP-2.2: we now have a kmer_send_buffer 
            std::vector<kmer_pair_t> &kmer_send_buffer = 
                kmer_send_buffer_[round_id%2UL];

            kmer_send_buffer.resize(cum);
            // do a parallel write into kmer_send_buffer //
#pragma omp parallel 
            {
              size_t my_id = omp_get_thread_num(); 
              size_t start_idx = kmer_starts[my_id];
              size_t end_idx = start_idx + uniq_kmer_counts[my_id];
              size_t kmer_send_buffer_idx = prefix_sum[my_id];

              for (size_t i=start_idx; i<end_idx; i++) {
                if (kmer_send_buffer_idx >= kmer_send_buffer.size()) {
                  throw std::logic_error("kmer_send_buffer invariant failed " 
                      " tcount="+std::to_string(tcount) +
                      "cum="+std::to_string(cum)
                      +" kmer_send_buffer_idx="+
                        std::to_string(kmer_send_buffer_idx)
                      +" kmer_send_buffer.size()="+
                        std::to_string(kmer_send_buffer.size())
                     );
                }
                kmer_send_buffer[kmer_send_buffer_idx++] = kmer_pair_buffer[i];
              }
            }
            ////////////////////////////////////////////////////////////////////
            is_send_buffer_full = true;
          }

          ++round_id;
        } while (next_avail_idx >= max_kmers_);

        all_kmers_processed_.store(true);
        return NULL;
      }


      //////////////////////////////////////////////////////////////////////////
      inline const char* skip_till_symbol(const char *beg, const char *end, 
          char symbol) {
        for (; (beg != end) && (*beg != symbol); ++beg) { }
        return beg;
      }


      size_t processor_count_;
      size_t rank_;
      size_t kmer_len_;
      size_t max_kmers_;
      char delimiter_;
      char skip_symbol_;
      const char * reads_beg_;
      const char * reads_end_;
      pthread_t comm_thread_;
      pthread_t comp_thread_;
      std::atomic<bool> is_send_buffer_full_[2UL];
      std::atomic<bool> all_kmers_processed_; 
      std::vector<int> send_counts_kmer_[2UL];
      std::vector<kmer_pair_t> kmer_send_buffer_[2UL];
      std::unordered_map<kmer_t, int> kmer_freq_map_;
      kmer_t kmer_mask_;
      int num_batch_transfers_;
  }; // class Async_Distributed_Kmer_Counter //

} // namespac esc //
} //namespace intel //
#endif
