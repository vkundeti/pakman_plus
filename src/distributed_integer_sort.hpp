#ifndef DISTRIBUTED_INTEGER_SORT_HPP
#define DISTRIBUTED_INTEGER_SORT_HPP

#include <algorithm>
#include <cstring>
#include <cassert>
#include <list>
#include <map>
#include <mpi.h>
#include <record_traits.hpp>
#include <most_significant_record_bit_and_byte.hpp>

namespace intel {
namespace esc {

// Distributed Integer Sort: does radix-sort on fixed-width strings (integers)
// from some alphabet. Note we process the alphabets from left->right since this
// has a non-recursive structure and can be repeated for the next alphabet/byte 
// without worrying about sorting a sub-group.
template<typename RecordTraits>
class Distributed_Integer_Sort {

  public:
    ////////////////////////////////////////////////////////////////////////////
    typedef RecordTraits record_traits;
    typedef typename record_traits::record_t record_t;
    typedef Most_Significant_Record_Bit_Byte<record_traits> msbb_finder_t;
    typedef typename msbb_finder_t::msbb_t msbb_t;

    // TODO(vamsikku): PersistentBucket: hold fixed number of elements incore
    // and rest out-of-core, this guards against skewed distribution in
    // the buckets //
    struct PersistentBucket {
      //////////////////////////////////////////////////////////////////////////
      typedef std::vector<record_t> incore_elements_t;
      typedef typename incore_elements_t::iterator incore_element_iterator_t;
      typedef typename incore_elements_t::const_iterator
          const_incore_element_iterator_t;

      //////////////////////////////////////////////////////////////////////////
      void add_element(const record_t& r) { elements_.push_back(r); }
      size_t size() const { return elements_.size(); }
      void clear() { elements_.clear(); }
      record_t const * data() const { return elements_.data(); }

      const incore_elements_t& incore_elements() const { return elements_; }

      const_incore_element_iterator_t begin() const {
        return elements_.begin();
      }
      bool empty() const { return elements_.empty(); }

      const_incore_element_iterator_t end() const {
        return elements_.end();
      }

      void remove_duplicates_in_bucket() {
        // TODO(vamsikku): need to expand traits to support this //
      }

      incore_elements_t elements_; 
      size_t element_count_;
      size_t max_incore_elements_;
    }; // class PersistentBucket //
    typedef PersistentBucket bucket_t;
    typedef std::map<size_t, bucket_t> dynamic_buckets_t;
    ////////////////////////////////////////////////////////////////////////////

    //NOTE: currently the max_bucket_bytes_=4 this means the number of
    //processors should be <= 2^32 (4 billion processors) //
    Distributed_Integer_Sort(size_t N=1, size_t my_rank=0,
        size_t bucket_bytes=1,
        size_t msb_byte_index=std::numeric_limits<size_t>::max(), 
        size_t incore_bucket_size_in_mb=1000,
        const char *persistence_root="./")
      : N_(N), my_rank_(my_rank), bucket_bytes_(bucket_bytes), 
        incore_bucket_size_in_mb_(incore_bucket_size_in_mb),
        persistence_root_(persistence_root), max_bucket_bytes_(7UL),
        buckets_(), received_keys_(), max_buckets_per_processor_(),
        current_round_(), max_rounds_(), shuffle_mode_(false), msbb_(),
        effective_record_len_(record_traits::record_length_in_bytes()),
        N_bits_(0)
        {
          init(msb_byte_index);
          //compute_processor_bits();
        }

    void init(size_t msb_byte_index) {
      size_t record_len_in_bytes = record_traits::record_length_in_bytes();

      effective_record_len_ = record_len_in_bytes;
      //////////////////////////////////////////////////////////////////////////
      // setup MIN and MAX bucket_bytes_ //
      bucket_bytes_ = std::max(bucket_bytes_, 1UL);
      max_bucket_bytes_ = std::min(record_len_in_bytes, max_bucket_bytes_); 
      bucket_bytes_ = std::min(bucket_bytes_, max_bucket_bytes_);
      //////////////////////////////////////////////////////////////////////////

      size_t bits = 8UL*bucket_bytes_; 
      bits = std::max(bits, 8UL);

      // 1. Number of buckets //
      number_of_buckets_ = (1UL << bits);

      // increase bucket_bytes_ so that N_ > number_of_buckets_ //
      while( (N_ > number_of_buckets_) && 
             (bucket_bytes_ < max_bucket_bytes_) ) {
        ++bucket_bytes_;
        bits = 8UL*bucket_bytes_;
        number_of_buckets_ = (1ULL << bits);
      }

      
      if (N_ > number_of_buckets_) {
        throw std::logic_error("INVARIANT VIOLATION: bucket processor "
              "invariant failed processors="+std::to_string(N_)+
              " buckets="+std::to_string(number_of_buckets_));
      }

      if (!my_rank_) {
        printf("[Distributed_Integer_Sort]: number_of_buckets_=%lu "
              "bucket_bytes_=%lu\n", number_of_buckets_, bucket_bytes_);
      }

      // 2. Max buckets per processor (distribute buckets uniformly)//
      max_buckets_per_processor_ = number_of_buckets_/N_;

      //3. Max rounds (atleast 1) //
      max_rounds_ = record_len_in_bytes/bucket_bytes_;
      if (record_len_in_bytes%bucket_bytes_) { ++max_rounds_; }

      max_rounds_ = std::min( msb_byte_index, max_rounds_);
      N_bits_ = 8*bucket_bytes_;
      N_bits_ = std::max(N_bits_, 8UL);

      if (!my_rank_) {
        printf("[INFO]: total-buckets=%lu max-buckets-per-processor=%lu "
            " max_rounds=%lu N_bits_=%lu\n",
            number_of_buckets_, max_buckets_per_processor_, max_rounds_,
            N_bits_);
      }
    }

  private:

    size_t current_byte_start() const {
      size_t v = current_round_*bucket_bytes_;
      if (!shuffle_mode_ && (v >= effective_record_len_)) {
        throw std::logic_error("[current_byte_start]: invariant violation "
            " effective_rlen="+std::to_string(effective_record_len_)+
            " max_rounds="+std::to_string(max_rounds_)+
            " curr_round="+std::to_string(current_round_));
      }
      return v;
    }

    size_t current_byte_end() const {
      // for shuffle and sort we should not go beyond the effective_record_len_
      size_t v = ((current_round_*bucket_bytes_) + bucket_bytes_);
      return shuffle_mode_ ? v : std::min( v, effective_record_len_);
    }

    template<typename PartIterator>
    void do_local_distribution(PartIterator& part_records_begin, 
        PartIterator& part_records_end) {
      size_t bucket_idx;
      for ( ;part_records_begin!=part_records_end; ++part_records_begin) {
        const record_t& record = *part_records_begin;
        if (!shuffle_mode_) {
          bucket_idx = record_traits::bytes_to_index(record,
              current_byte_start(), current_byte_end());
        } else {
          if (!msbb_.is_valid()) {
            throw std::logic_error("Invalid msbb for shuffle");
          }
          bucket_idx = msbb_finder_t::extract_bits_into_size_t(record, msbb_,
              N_bits_);
        }
        (buckets_[current_round_%2UL])[ bucket_idx ].add_element(record);
      }
    }


  public:

    template<typename PartitionRecordIterator, typename PartitionOutputIterator>
    void sort(
        PartitionRecordIterator part_records_begin,
        PartitionRecordIterator part_records_end, 
        PartitionOutputIterator output, bool remove_duplicates=false,
        size_t start_round_idx=0UL, bool no_output=false) {


      // TODO(vamsikku): assert the underlying type of the iterator is same
      // as the record type //
      clear_buckets();
      current_round_ = start_round_idx;
      // initial distribution //
      do_local_distribution(part_records_begin, part_records_end);

      do {
#if 0
        //TODO(vamsikku): we can get rid of duplicates 
        if (remove_duplicates) {
          remove_duplicate_keys();
        }
#endif
        this->send_and_receive_keys_across_processors();
        this->distribute_received_keys_locally();
        // send // 
      } while (++current_round_ < max_rounds_);

      if (!no_output)
        report_output(output, remove_duplicates);
    }


    template<typename PartitionRecordIterator, typename PartitionOutputIterator>
    void sort_and_shuffle_using_msbb(
        const msbb_t& msbb,
        PartitionRecordIterator part_records_begin,
        PartitionRecordIterator part_records_end, 
        PartitionOutputIterator output, bool remove_duplicates=false,
        size_t start_round_idx=0UL) {
      size_t record_len_in_bytes = record_traits::record_length_in_bytes();
      msbb_ = msbb;
      if ((size_t) msbb.byte_index_ >= record_len_in_bytes) {
        throw std::logic_error("[MSBB Invariant violation]: record_len="+
            std::to_string(record_len_in_bytes)+" msbb.byte_idx="+
            std::to_string(msbb.byte_index_));
      }


      // We need to sort between (msbb.byte_index_ .... 0] which is a
      // record of length msbb.byte_index_ //
      {
        effective_record_len_ = msbb.byte_index_;
        //Invariant: effective_record_len_ >= bucket_bytes_;
        bucket_bytes_ = std::min(effective_record_len_, bucket_bytes_);
        if (bucket_bytes_) {
          max_rounds_ = effective_record_len_/bucket_bytes_;
          if (effective_record_len_%bucket_bytes_){ ++max_rounds_; }
        } else {
          max_rounds_ = 0;
        }
        N_bits_ = 8*bucket_bytes_;
        N_bits_ = std::max(N_bits_, 8UL);
      }

      if (!my_rank_) {
        printf("[shuffle_sort]: max_rounds=%lu\n", max_rounds_);
      }

      if (max_rounds_) { // atleast one rounds of sorting//
        shuffle_mode_ = false;
        sort(part_records_begin, part_records_end, output, remove_duplicates,
          start_round_idx, true);
      } else { // no rounds just shuffle on the byte //
        // just distribute keys //
        shuffle_mode_ = true;
        do_local_distribution(part_records_begin, part_records_end);
      }

      shuffle_mode_ = true;
      this->send_and_receive_keys_across_processors();
      this->distribute_received_keys_locally();
      ++current_round_;
      report_output(output, remove_duplicates);
    }


    template<typename OutputIterator>
    void report_output(OutputIterator output, bool remove_dups=false) {
      dynamic_buckets_t& current_buckets = buckets_[(current_round_)%2];


      for (auto itr=current_buckets.begin(); itr!=current_buckets.end();
            ++itr) {
        bucket_t &bucket = itr->second;
        auto kitr=bucket.begin(), prev_kitr=kitr;
        for (++kitr; kitr!=bucket.end(); ++kitr, ++prev_kitr) {
          if (!remove_dups || 
                (!record_traits::equivalent_record(*prev_kitr, *kitr)) ) {
            *output = *prev_kitr;
            ++output;
          }
        }
        if (prev_kitr != bucket.end()) {
          *output = *prev_kitr;
          ++output;
        }
      }
    }

    void send_and_receive_keys_across_processors() {
#ifdef DEBUG_MPI
      printf("processor-%lu starting communication round=%lu\n",
          my_rank_, current_round_);
#endif

      received_keys_.clear();
      // take the keys from buckets_[current_round_%2] and put them
      // in buckets_[(current_round_+1)%2] //
      typedef typename std::list<MPI_Request>::iterator send_request_iterator_t;
      std::list<MPI_Request> send_requests;
      dynamic_buckets_t send_data; 

      dynamic_buckets_t& current_buckets = buckets_[(current_round_)%2];

      // STEP-1: distribute the keys in current bucket //
      for (auto itr=current_buckets.begin(); itr!=current_buckets.end();
            ++itr) {
        bucket_t &bucket = itr->second;
        // scan and distribute the keys in the current bucket//
        for (auto kitr=bucket.begin(); kitr!=bucket.end(); ++kitr) {
          auto map_result = map_record_to_processor(*kitr);
          size_t processor_id = map_result.first;
          size_t bucket_idx = map_result.second;

          if (processor_id >= N_) {
            throw std::logic_error("Invalid processord_id=" +
                  std::to_string(processor_id) + " bucket_idx="+
                  std::to_string(bucket_idx) );
          }

          if (processor_id == my_rank_) {
            received_keys_[processor_id].add_element(*kitr);
          } else {
            send_data[processor_id].add_element(*kitr);
          }
        }
        bucket.clear();
      }
      current_buckets.clear();


      //STEP-2: create isend requests for each of the processors //
      for (auto itr=send_data.begin(); itr!=send_data.end(); ++itr) {
        if ((itr->second).empty()) {
          throw std::logic_error("[INVARIANT VIOLATION]: send data "
                "cannot be empty");
        }
        size_t processor_id = itr->first; 
        send_records_using_fixed_buffer_width( (itr->second).incore_elements(),
            (int)processor_id, (int) current_round_, send_requests);
      }

      //STEP-3: use sparse data exchange to receive data and distribute keys//
      MPI_Request ibarrier_request;
      bool entered_barrier_phase = false;
      do { // communication protocol: sparse data exchange  //
        //////////////////////////////////////////////////////////////////////
        // 1. Probe for a message and receive it if there is one            //
        //////////////////////////////////////////////////////////////////////
        MPI_Status status;
        int flag = 0;
        int ret = MPI_Iprobe(MPI_ANY_SOURCE, (int) current_round_,
              MPI_COMM_WORLD, &flag, &status);
        if ((ret == MPI_SUCCESS) && flag) {
          //  1.1 If there is some message then receive it //
          int recv_byte_count=0;

          //  1.2 Get the count of bytes //
          ret = MPI_Get_count(&status, MPI_BYTE, &recv_byte_count);
          if ((ret != MPI_SUCCESS) || (recv_byte_count <0) ) {
            throw std::logic_error("MPI_Get_count: failed " 
                "or invalid recv_byte_count " +
                  std::to_string(recv_byte_count));
          }

          // 1.3 Receive //
          std::vector<unsigned char> receive_buffer(recv_byte_count, '\0');
          ret = MPI_Recv(receive_buffer.data(), recv_byte_count, MPI_BYTE,
                status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD,
                  MPI_STATUS_IGNORE);
          
          if (ret != MPI_SUCCESS) {
            throw std::logic_error("MPI_Recv: failed ");
          }


          add_to_received_keys(status.MPI_SOURCE,
              receive_buffer.data(),
              receive_buffer.data()+receive_buffer.size());

#if 0
          // 1.4 Distribute the records of this processor to local buckets // 
          read_and_distribute_records_in_buffer_locally(
              receive_buffer.data(),
              receive_buffer.data()+receive_buffer.size());
#endif

        } else if (ret != MPI_SUCCESS) {
          throw std::logic_error("MPI_Iprobe: failed");
        }
        //////////////////////////////////////////////////////////////////////

        if (!entered_barrier_phase) {
          //////////////////////////////////////////////////////////////////////
          // 2. Check if all my sent messages are received by some one.       //
          //////////////////////////////////////////////////////////////////////
          if (!send_requests.empty()) {
            std::vector<send_request_iterator_t> received_requests;
            for (send_request_iterator_t sitr=send_requests.begin();
                  sitr!=send_requests.end(); ++sitr) {
              int ret = MPI_Test( &(*sitr), &flag, MPI_STATUS_IGNORE);

              if ((ret == MPI_SUCCESS) && flag) {
                received_requests.push_back(sitr);
              } else if (ret != MPI_SUCCESS) {
                throw std::logic_error("MPI_Test: failed ");
              }
            }
            for (auto itr=received_requests.begin();
                  itr!=received_requests.end(); ++itr) {
              send_requests.erase(*itr);
            }
          }

          ////////////////////////////////////////////////////////////////////
          // 3. Go into barrier phase                                       //
          ////////////////////////////////////////////////////////////////////
          if (send_requests.empty()) {
            int ret = MPI_Ibarrier(MPI_COMM_WORLD, &ibarrier_request);
            if (ret != MPI_SUCCESS) {
              throw std::logic_error("MPI_Ibarrier: failed");
            }
            entered_barrier_phase = true;
#ifdef DEBUG_MPI
            printf("processor-%lu touched barrier\n", my_rank_);
#endif

          } 
        }


        if (entered_barrier_phase) {
          // 3. Test for completion of the barrier call //
          int flag = 0;
          int ret = MPI_Test( &ibarrier_request, &flag, MPI_STATUS_IGNORE);
          if ((ret == MPI_SUCCESS) && (flag==1)) {
            break; // terminate communication //
          } else if (ret != MPI_SUCCESS) {
            throw std::logic_error("MPI_Test: failed");
          }
        }

      } while(1); // sparse data exchange protocol //

#ifdef DEBUG_MPI
      printf("processor-%lu done communication round=%lu\n", my_rank_,
          current_round_);
#endif


    }

    size_t map_bucket_to_processor(size_t i) const {
      if (i > number_of_buckets_) { 
        throw std::logic_error("[ERROR]: invalid bucket index\n");
      }
      return std::min( (i/max_buckets_per_processor_), N_-1);
    }

  private:

    void remove_duplicate_keys() {
      dynamic_buckets_t& current_buckets = buckets_[(current_round_)%2];
      for (auto itr=current_buckets.begin(); itr!=current_buckets.end();
            ++itr) {
        (itr->second).remove_duplicates_in_bucket();
      }
    }

    void distribute_received_keys_locally() {
      dynamic_buckets_t& next_buckets = buckets_[(current_round_+1UL)%2];
      next_buckets.clear();
      size_t bucket_idx;

      for (auto itr=received_keys_.begin(); itr!=received_keys_.end(); ++itr) {
        bucket_t & bucket = itr->second;
        for (auto kitr=bucket.begin(); kitr!=bucket.end(); ++kitr) {
          if (!shuffle_mode_) {
            bucket_idx = record_traits::bytes_to_index(*kitr,
                current_byte_start(), current_byte_end());
          } else {
            if (!msbb_.is_valid()) {
              throw std::logic_error("Invalid msbb for shuffle");
            }
            bucket_idx = msbb_finder_t::extract_bits_into_size_t(*kitr, msbb_,
                N_bits_);
          }
          next_buckets[bucket_idx].add_element(*kitr);
        }
      }
    }
  public:

    void print_byte_string(const record_t& v, FILE *fptr=stdout) {
      size_t rlen = record_traits::record_length_in_bytes(); 
      const unsigned char * rbytes = (const unsigned char *)(&v);
      fprintf(fptr, " ");
      for (size_t i=0; i<rlen; i++) {
        size_t val =  (size_t) rbytes[i];
        char bit=(1 << 7);
        fprintf(fptr, " %lu[", v);
        for (size_t k=0; k<8; k++, bit >>= 1) {
          fprintf(fptr, "%d", bit&val ? 1 : 0);
        }
        fprintf(fptr, "] ");
      }
      fprintf(fptr, " ");
    }

    void print_dynamic_buckets(FILE *fptr=stdout) {
      fprintf(fptr, "curr_round=%lu\n", current_round_);
      dynamic_buckets_t& next_buckets = buckets_[current_round_%2];

      for (auto itr=next_buckets.begin(); itr!=next_buckets.end(); ++itr) {
        fprintf(fptr, "bucket[%lu]: {", itr->first);
        for (auto r : itr->second) {
          //this->print_byte_string(r, fptr);
          fprintf(fptr, "%lu ", r);
        }
        fprintf(fptr, "}\n");
      }

    }

  private:


    template<typename records_t, typename send_request_container_t>
    void send_records_using_fixed_buffer_width(const records_t& records_in,
        int processor_id, int tag,
        send_request_container_t& send_request_container,
        size_t max_buffer_size_mb=1024) {
      record_t const *records = records_in.data();
      //TODO(vamsikku): check for overflows //
      size_t mb = (1UL<<20);
      size_t byte_buffer_size =  max_buffer_size_mb*mb;
      byte_buffer_size = std::min((size_t)(std::numeric_limits<int>::max()),
            byte_buffer_size);

      // max number of records in the buffer //
      size_t rec_len = record_traits::record_length_in_bytes();
      if (rec_len > byte_buffer_size) {
        throw std::logic_error("[INVARIANT VIOLATION]: byte_buffer must be"
              " atleast the record size");
      }

      size_t max_records_per_message = byte_buffer_size/rec_len;
      max_records_per_message =
          std::min(max_records_per_message, records_in.size());

      int mesg = 0;
      for (size_t i=0; i<records_in.size(); i+=max_records_per_message) {
        size_t records_in_this_mesg = std::min(max_records_per_message,
            (records_in.size()-i));
        int send_count = (int)(records_in_this_mesg)*(int)(rec_len);
        ////////////////////////////////////////////////////////////////////////
        send_request_container.push_back(MPI_REQUEST_NULL);
        int ret = MPI_Issend((void *)(records+i), send_count, MPI_BYTE,
            (int)processor_id, (int) tag, MPI_COMM_WORLD,
            &send_request_container.back());
        if (ret != MPI_SUCCESS) {
          send_request_container.pop_back();
          throw std::logic_error("[MPI_Issend] failed");
        }
        ////////////////////////////////////////////////////////////////////////
        ++mesg;
      }

#ifdef DEBUG_MPI
      printf("p-%lu->p-%d max_records=%lu bucket=%lu mesgs=%d \n",
          my_rank_, processor_id, max_records_per_message, records_in.size(),
            mesg);
#endif

    }

    template<typename RecordType>
    std::pair<size_t, size_t> map_record_to_processor(const RecordType& record){
      size_t processor_id, bucket_idx;
      if (!shuffle_mode_) {
        bucket_idx = record_traits::bytes_to_index(record,
            current_byte_start(), current_byte_end());
      } else {
        if (!msbb_.is_valid()) {
          throw std::logic_error("Invalid msbb for shuffle");
        }
        bucket_idx = msbb_finder_t::extract_bits_into_size_t(record, msbb_,
            N_bits_);
      }

      processor_id = map_bucket_to_processor(bucket_idx);
      return std::make_pair(processor_id, bucket_idx);
    }

    template<typename ByteIterator>
    void add_to_received_keys(size_t processor_id, 
          ByteIterator bbegin, ByteIterator bend) {
      if (bbegin > bend) {
        throw std::logic_error("[INVARIANT VIOLATION]: bbegin > bend");
      }

      auto itr = received_keys_.find(processor_id);
      if (itr == received_keys_.end()) {
        itr = received_keys_.insert(
                std::make_pair(processor_id, bucket_t())
              ).first;
      }
      bucket_t &bucket = itr->second;
      size_t rlen = record_traits::record_length_in_bytes();
      record_t record;
      size_t bucket_idx;
      for (; bbegin!=bend; bbegin+=rlen ) {
        record_traits::create_record(record, bbegin, bbegin+rlen);
        if (!shuffle_mode_) {
          bucket_idx = record_traits::bytes_to_index(record,
              current_byte_start(), current_byte_end());
        } else {
          if (!msbb_.is_valid()) {
            throw std::logic_error("Invalid msbb for shuffle");
          }
          bucket_idx = msbb_finder_t::extract_bits_into_size_t(record, msbb_,
              N_bits_);
        }
        size_t pid = map_bucket_to_processor(bucket_idx);


        if ((pid != my_rank_) && !shuffle_mode_) {
          std::string mesg;
          mesg += "[INVARIANT VIOLATION]: processor_id : " +
              std::to_string(my_rank_) + " got key for : " +
                std::to_string(pid) + " round="+std::to_string(current_round_)
                + " record="+std::to_string(record) + " from processor:"+
                  std::to_string(processor_id); 
          fflush(stdout);
          throw std::logic_error(mesg);
        }


        bucket.add_element(record);
      }
    }

    void compute_processor_bits() {
      size_t pow_two = 1;
      N_bits_ = 1;
      while (pow_two < N_) {
        pow_two <<= 1;
        N_bits_++;
      }
      N_bits_ = std::max(8UL, N_bits_);
      printf("N_bits=%lu\n", N_bits_);
    }

    void clear_buckets() { buckets_[0].clear(); buckets_[1].clear(); }

    size_t N_; // N is the total number of processors
    size_t my_rank_; // current processor //
    size_t bucket_bytes_; 
    size_t incore_bucket_size_in_mb_;
    const char *persistence_root_;
    size_t max_bucket_bytes_;
    dynamic_buckets_t buckets_[2UL];
    dynamic_buckets_t received_keys_;
    size_t number_of_buckets_;
    size_t max_buckets_per_processor_;
    size_t current_round_;
    size_t max_rounds_;
    bool shuffle_mode_;
    msbb_t msbb_;
    size_t effective_record_len_;
    size_t N_bits_;
}; // class Distributed_Integer_Sort //


typedef Distributed_Integer_Sort<size_t_record_traits>
  distributed_sorter_size_t;


typedef Distributed_Integer_Sort<int_record_traits> distributed_sorter_int_t;

} // namespace esc   //
} // namespace intel //



#endif
