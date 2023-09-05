#ifndef CONCURRENT_KMER_COUNTER
#define CONCURRENT_KMER_COUNTER

#include <memory>
#include <unordered_map>
#include <algorithm>
#include <map>

#include "stop_clock.hpp"

#include <oneapi/tbb/concurrent_hash_map.h>
#include <oneapi/tbb/concurrent_unordered_map.h>
#include <oneapi/tbb/concurrent_vector.h>
#include <oneapi/tbb/parallel_sort.h>

#include <omp.h>

namespace intel {
namespace esc {

template<typename kmer_t,
         template<typename T> class allocator_t = std::allocator>
class Concurrent_Kmer_Partial_Counter {

  private:
    ////////////////////////////////////////////////////////////////////////////
    template<typename T1>
    using alloc_t = allocator_t<T1>;

    struct kmer_pair_t {
      kmer_t kmer_;
      int freq_; 
    };  // struct kmer_pair_t //

    struct default_processor_t {
      template<typename kmer_pair_iterator_t>
      inline void operator() (kmer_pair_iterator_t begin,
            kmer_pair_iterator_t end) const { }
    }; // struct default_processor_t //
    typedef std::vector< kmer_pair_t, alloc_t<kmer_pair_t> >
        partial_count_buffer_t;
    typedef oneapi::tbb::concurrent_hash_map<kmer_t, int> concurrent_hash_map_t;
    typedef std::unordered_map<kmer_t, int> hash_map_t;
    typedef typename concurrent_hash_map_t::accessor accessor_t;
    ////////////////////////////////////////////////////////////////////////////

  public:

    Concurrent_Kmer_Partial_Counter(size_t klen=32,
        const char *symbol_begin=NULL, const char *symbol_end=NULL,
        size_t read_len=100, size_t max_kmer_count=1000000,
        char delimiter='\n', char skip_symbol='>')
      : klen_(klen), symbol_begin_(symbol_begin), symbol_end_(symbol_end),
        read_len_(read_len), max_kmer_count_(max_kmer_count),
        delimiter_(delimiter), skip_symbol_(skip_symbol),
        partial_count_buffer_(), concurrent_hash_map_(NULL) {}

    ~Concurrent_Kmer_Partial_Counter() {
      if (concurrent_hash_map_) { 
        delete concurrent_hash_map_; 
        concurrent_hash_map_ = NULL;
      }
    }

    void reset(char const *begin, char const *end) {
      symbol_begin_ = begin;
      symbol_end_ = end;
    }


    template<typename processor_t=default_processor_t>
    size_t concurrent_count_and_process_omp_v1(processor_t& processor) {
      if (!concurrent_hash_map_) {
        // TODO(vamsikku): change the starting concurrent hash size //
        concurrent_hash_map_ = new concurrent_hash_map_t(40000000UL);
        if (!concurrent_hash_map_) {
          throw std::logic_error(
              "ERROR: unable to create a concurrent hash map");
        }
      }

      //TODO(vamsikku): min-kmers per thread should be a parameter //
      size_t min_len=3UL;
      size_t rid = 0;
#pragma omp parallel
      {
        size_t my_id = omp_get_thread_num(); 
        size_t tcount = omp_get_num_threads();
        if (!my_id) {
          printf("tcount=%lu\n", tcount);
        }
        size_t chunk_size = std::max(read_len_/tcount, min_len);
        const char *symbol_begin = symbol_begin_;
        const char *symbol_end = symbol_end_;
        while (symbol_begin != symbol_end) {
          // STEP-1: set the start of the 
          symbol_begin =
              skip_till_symbol(symbol_begin, symbol_end, skip_symbol_);
          if (symbol_begin == symbol_end) { break; }
          ++symbol_begin;

          // STEP-2: skip till delimiter //  
          symbol_begin =
              skip_till_symbol(symbol_begin, symbol_end, delimiter_);
          if (symbol_begin == symbol_end) { break;}
          ++symbol_begin;

          size_t read_len = read_len_;
          const char *read_beg = symbol_begin;
          const char *read_end = symbol_begin+read_len;

          if (read_end > symbol_end) {
            read_end = symbol_end; 
            read_len = symbol_end-read_beg;
          }
          if (read_len < klen_) { continue; } // read small //

          // adjust sliding window indicies //
          const char *window_end = read_beg + (read_len -klen_ +1UL);
          {
            //Invariant: wbeg+k <= read_end //
            do {
              const char *wbeg = read_beg + (my_id * chunk_size);
              const char *wend = wbeg+chunk_size;

              if (wbeg >= window_end) { break; }
              if (wend >= window_end) { wend = window_end; }


              // setup first k-mer //
              kmer_t kmer=kmer_t(0UL);
              const char *wbeg_plus_k = wbeg;
              for (size_t i=0; i<(klen_-1UL); i++, wbeg_plus_k++) {
                kmer <<= 2UL;
                kmer |= symbol_to_kmer(*wbeg_plus_k);
              }

              // slide //
              while (wbeg < wend) {
                kmer <<= 2UL;
                kmer |= symbol_to_kmer(*wbeg_plus_k);
                {
                  accessor_t accessor;
                  concurrent_hash_map_->insert(accessor, kmer);
                  accessor->second++;
                  accessor.release();
                }
                ++wbeg; ++wbeg_plus_k; 
              }
            } while (0);
          }

#if 0
          //TODO(vamsikku): if the size of hashmap exceeds threshold all
          //threads must wait until its flushed //
          if (concurrent_hash_map_->size() > max_kmer_count_) {
            processor(concurrent_hash_map_->begin(),
                  concurrent_hash_map_->end());
            concurrent_hash_map_->clear();
          }
#endif


          if (!my_id) {
            ++rid;
          }
        } // foreach read //

      }


      if (concurrent_hash_map_) {
        // trailing case //
        processor(concurrent_hash_map_->begin(), concurrent_hash_map_->end());
        concurrent_hash_map_->clear();
        delete concurrent_hash_map_;
        concurrent_hash_map_ = NULL;
      }
      return rid;
    }

    //ALGO_v0: partition on the k-mer space //
    template<typename processor_t=default_processor_t>
    size_t concurrent_count_and_process_omp_v0(processor_t& processor) {
      if (!concurrent_hash_map_) {
        // TODO(vamsikku): change the starting concurrent hash size //
        concurrent_hash_map_ = new concurrent_hash_map_t(1000000UL);
        if (!concurrent_hash_map_) {
          throw std::logic_error(
              "ERROR: unable to create a concurrent hash map");
        }
      }

      //TODO(vamsikku): min-kmers per thread should be a parameter //
      size_t tcount, min_len=3UL;
#pragma omp parallel
      {
        if (!omp_get_thread_num() ) {
          tcount = omp_get_num_threads();
        }
      }
      printf("tcount=%lu\n", tcount);
      size_t chunk_size = std::max(read_len_/tcount, min_len);

      size_t rid = 0;
      while (symbol_begin_ != symbol_end_) {
        // STEP-1: set the start of the 
        symbol_begin_ =
            skip_till_symbol(symbol_begin_, symbol_end_, skip_symbol_);
        if (symbol_begin_ == symbol_end_) { break; }
        ++symbol_begin_;

        // STEP-2: skip till delimiter //  
        symbol_begin_ =
            skip_till_symbol(symbol_begin_, symbol_end_, delimiter_);
        if (symbol_begin_ == symbol_end_) { break;}
        ++symbol_begin_;

        size_t read_len = read_len_;
        const char *read_beg = symbol_begin_;
        const char *read_end = symbol_begin_+read_len;

        if (read_end > symbol_end_) {
          read_end = symbol_end_; 
          read_len = symbol_end_-read_beg;
        }
        if (read_len < klen_) { continue; } // read small //

        // adjust sliding window indicies //
        const char *window_end = read_beg + (read_len -klen_ +1UL);

#pragma omp parallel
        {
          //Invariant: wbeg+k <= read_end //
          size_t my_kmers = 0UL;
          do {
            size_t my_id = omp_get_thread_num(); 
            const char *wbeg = read_beg + (my_id * chunk_size);
            const char *wend = wbeg+chunk_size;

            if (wbeg >= window_end) { break; }
            if (wend >= window_end) { wend = window_end; }


            // setup first k-mer //
            kmer_t kmer=kmer_t(0);
            const char *wbeg_plus_k = wbeg;
            for (size_t i=0; i<(klen_-1UL); i++, wbeg_plus_k++) {
              kmer <<= 2UL;
              kmer |= symbol_to_kmer(*wbeg_plus_k);
            }

            // slide //
            while (wbeg < wend) {
              kmer <<= 2UL;
              kmer |= symbol_to_kmer(*wbeg_plus_k);
              {
                accessor_t accessor;
                concurrent_hash_map_->insert(accessor, kmer);
                accessor->second++;
                accessor.release();
              }
              ++wbeg; ++wbeg_plus_k; ++my_kmers;
            }
          } while (0);
        }

#if 0
        if (concurrent_hash_map_->size() > max_kmer_count_) {
          processor(concurrent_hash_map_->begin(), concurrent_hash_map_->end());
          concurrent_hash_map_->clear();
        }
#endif
        ++rid;
      } // foreach read //

      if (concurrent_hash_map_) {
        // trailing case //
        processor(concurrent_hash_map_->begin(), concurrent_hash_map_->end());
        concurrent_hash_map_->clear();
        delete concurrent_hash_map_;
        concurrent_hash_map_ = NULL;
      }
      return rid;
    }

    struct implicit_kmer_compare_t {

      implicit_kmer_compare_t(const char *read_head,
          size_t kcol_size, size_t read_len, size_t klen)
        : read_head_(read_head), kcol_size_(kcol_size), read_len_(read_len),
          klen_(klen) {}

      inline bool operator()(
            const unsigned int& a, const unsigned int& b) const {

        size_t row_a = size_t(a)/kcol_size_;
        size_t col_a = size_t(a)%kcol_size_; 

        size_t row_b = size_t(b)/kcol_size_;
        size_t col_b = size_t(b)%kcol_size_; 

        const char *kmer_a = read_head_ + ((read_len_*row_a) + col_a);
        const char *kmer_b = read_head_ + ((read_len_*row_b) + col_b);

        for (size_t k=0; k<klen_; k++) {
          if (kmer_a[k] == kmer_b[k]) { continue;}
          return (kmer_a[k] < kmer_b[k]);
        }
        return false;
      }

      inline void get_kmer(size_t a, std::string& kmer) const {
        kmer.clear();

        size_t row_a = size_t(a)/kcol_size_;
        size_t col_a = size_t(a)%kcol_size_; 

        const char *kmer_a = read_head_ + ((read_len_*row_a) + col_a);

        for (size_t i=0; i<klen_; i++) {
          kmer.push_back(kmer_a[i]);
        }
      }

      inline bool equiv(unsigned int a, unsigned int b) {
        size_t row_a = size_t(a)/kcol_size_;
        size_t col_a = size_t(a)%kcol_size_; 

        size_t row_b = size_t(b)/kcol_size_;
        size_t col_b = size_t(b)%kcol_size_; 

        const char *kmer_a = read_head_ + ((read_len_*row_a) + col_a);
        const char *kmer_b = read_head_ + ((read_len_*row_b) + col_b);
        return !strncmp(kmer_a, kmer_b, klen_);
      }

      const char * const read_head_;
      size_t kcol_size_;
      size_t read_len_;
      size_t klen_;
    }; // struct implicit_kmer_compare_t //

    struct implicit_kmer_hash_t {
      implicit_kmer_hash_t(size_t klen) : klen_(klen) {}

      inline size_t operator() (const char *kmer_ptr) const {
        size_t hash_value = 0UL;
        for (size_t i=0UL; i<klen_; ++i, ++kmer_ptr) {
          hash_value += (((size_t) *kmer_ptr) + 1UL)%5UL;
          hash_value <<= 2UL;
        }
        return hash_value;
      }
      size_t klen_;
    }; // struct implicit_kmer_hash_t //




    struct kmer_strncmp_t {
      kmer_strncmp_t(size_t klen) : klen_(klen) {}

      inline bool operator() (const char* a, const char*b) const {
        return strncmp(a, b, klen_) < 0; 
      }

      inline bool operator() (const std::pair<const char *, size_t>& a,
          const std::pair<const char *, size_t>& b) const {
        return strncmp(a.first, b.first, klen_) < 0;
      }

      size_t klen_;
    }; // struct kmer_strncmp_t //

    struct implicit_kmer_equal_t {
      implicit_kmer_equal_t(size_t klen) : klen_(klen) {}

      inline bool operator() (const char *a, const char *b) const {
        return !strncmp(a, b, klen_);
      }
      size_t klen_;
    }; // struct kmer_strncmp_equiv_t //


    template<typename BitContainer>
    struct bit_kmer_strncmp_t {
      typedef typename BitContainer::const_iterator const_iterator_t;
      bit_kmer_strncmp_t(const BitContainer& bit_container_even, 
          const BitContainer& bit_container_odd, size_t klen)
          : bit_container_even_(bit_container_even), 
            bit_container_odd_(bit_container_odd), klen_(klen) {}

      inline bool operator() (size_t a, size_t b) const {
        const_iterator_t even_itr_a = bit_container_even_.cbegin() + a;
        const_iterator_t even_itr_b = bit_container_even_.cbegin() + b;

        const_iterator_t odd_itr_a = bit_container_odd_.cbegin() + a;
        const_iterator_t odd_itr_b = bit_container_odd_.cbegin() + b;

        for (size_t i=0; i<klen_; i++) {
          // [2*i, 2*i+1] (odd-bit, even-bit) //
          if ( *even_itr_a != *even_itr_b) {
            return (*even_itr_a) < (*even_itr_b);
          } else if (*odd_itr_a != *odd_itr_b) {
            return (*odd_itr_a) < (*odd_itr_b); 
          }

          ++even_itr_a; ++odd_itr_a;
          ++even_itr_b; ++odd_itr_b;
        }
        return false;
      }

      inline bool equiv(size_t a, size_t b) const {
        const_iterator_t even_itr_a = bit_container_even_.cbegin() + a;
        const_iterator_t even_itr_b = bit_container_even_.cbegin() + b;

        const_iterator_t odd_itr_a = bit_container_odd_.cbegin() + a;
        const_iterator_t odd_itr_b = bit_container_odd_.cbegin() + b;

        for (size_t i=0; i<klen_; i++) {

          // [2*i, 2*i+1] (odd-bit, even-bit) //
          if ( *even_itr_a != *even_itr_b) {
            return false;
          } else if (*odd_itr_a != *odd_itr_b) {
            return false;
          }

          ++even_itr_a; ++odd_itr_a;
          ++even_itr_b; ++odd_itr_b;
        }
        return true;
      }


      const BitContainer &bit_container_even_;
      const BitContainer &bit_container_odd_;
      size_t klen_;
    }; // struct bit_kmer_strncmp_t //

    template<typename processor_t=default_processor_t>
    size_t concurrent_count_and_process_omp_implicit_map(
        processor_t& processor) {
      ////////////////////////////////////////////////////////
      // STEP-0: find the number of bytes in all the reads //

      const char *symbol_begin = symbol_begin_;
      const char *symbol_end = symbol_end_;

      size_t data_size = ((size_t) symbol_end )-((size_t)symbol_begin);
      std::vector<char> read_data;

      read_data.resize( data_size, '\0' );
      char *read_buffer = read_data.data();
      char *curr_read_buffer_head = read_buffer;
      const char *read_buffer_end = read_data.data() + read_data.size();

      size_t read_count =0;
      while (symbol_begin < symbol_end) {
        // STEP-1: set the start of the 
        symbol_begin =
            skip_till_symbol(symbol_begin, symbol_end, skip_symbol_);
        if (symbol_begin >= symbol_end) { break; }
        ++symbol_begin;

        // STEP-2: skip till delimiter //  
        symbol_begin =
            skip_till_symbol(symbol_begin, symbol_end, delimiter_);
        if (symbol_begin >= symbol_end) { break;}
        ++symbol_begin;

        size_t read_len = read_len_;
        const char *read_beg = symbol_begin;
        const char *read_end = symbol_begin+read_len;

        if (read_end > symbol_end_) {
          read_end = symbol_end_; 
          read_len = symbol_end_-read_beg;
        }

        if (read_len != read_len_) {
          throw std::logic_error("Variable size read\n");
        }

        if (curr_read_buffer_head >= read_buffer_end) {
          throw std::logic_error("Invalid buffer size");
        }
        memcpy(curr_read_buffer_head, read_beg, read_len);
        curr_read_buffer_head += read_len;
        read_count++;
      }


      size_t kcol_size = (read_len_ - klen_ + 1UL); 
      implicit_kmer_compare_t kmer_compare(read_buffer, kcol_size, read_len_,
            klen_);
      std::map<unsigned int, size_t, implicit_kmer_compare_t> 
          kmer_freq(kmer_compare);

      size_t kmer_count = read_count*kcol_size;

      for (unsigned int i=0; i<kmer_count; i++) {
        auto itr = kmer_freq.find(i);
        if (itr == kmer_freq.end()) {
          itr = kmer_freq.insert(std::make_pair(i, 0UL)).first;
        }
        itr->second++;
      }

      printf("...done sorting...\n");
      printf("unique kmers=%lu\n", kmer_freq.size());
      return 0;

    }


    template<typename processor_t=default_processor_t>
    size_t concurrent_count_and_process_omp_explicit(processor_t& processor) {
      size_t tcount;
#pragma omp parallel
      {
        if (!omp_get_thread_num() ) {
          tcount = omp_get_num_threads();
        }
      }

      oneapi::tbb::concurrent_vector<kmer_t> kmer_buffer;
      std::atomic<size_t> total_reads = {0UL};
#pragma omp parallel
      {
        size_t my_id = omp_get_thread_num(); 
        size_t data_size = ((size_t) symbol_end_ )-((size_t)symbol_begin_);
        size_t chunk_size = data_size/tcount;

        const char *symbol_begin = symbol_begin_ + (my_id*chunk_size);
        const char *symbol_end = symbol_begin + chunk_size;

        if (symbol_begin >= symbol_end_) { 
          symbol_begin = symbol_end_; 
          symbol_end = symbol_end_;
        }
        if (symbol_end >= symbol_end_) {
          symbol_end = symbol_end_;
        }

        // adjust symbol_end //
        symbol_end = skip_till_symbol(symbol_end, symbol_end_, skip_symbol_);
        ++symbol_end;

        if (symbol_end >= symbol_end_) {
          symbol_end = symbol_end_;
        }


        size_t rid = 0;
        while (symbol_begin < symbol_end) {
          // STEP-1: set the start of the 
          symbol_begin =
              skip_till_symbol(symbol_begin, symbol_end, skip_symbol_);
          if (symbol_begin >= symbol_end) { break; }
          ++symbol_begin;

          // STEP-2: skip till delimiter //  
          symbol_begin =
              skip_till_symbol(symbol_begin, symbol_end, delimiter_);
          if (symbol_begin >= symbol_end) { break;}
          ++symbol_begin;

          size_t read_len = read_len_;
          const char *read_beg = symbol_begin;
          const char *read_end = symbol_begin+read_len;

          if (read_end > symbol_end_) {
            read_end = symbol_end_; 
            read_len = symbol_end_-read_beg;
          }
          if (read_len < klen_) { continue; } // read small //

          // adjust sliding window indicies //
          const char *window_end = read_beg + (read_len -klen_ +1UL);
          {
            //Invariant: wbeg+k <= read_end //
            const char *wbeg = read_beg; 
            const char *wend = window_end;

            // setup first k-mer //
            kmer_t kmer=kmer_t(0UL);
            const char *wbeg_plus_k = wbeg;
            for (size_t i=0; i<(klen_-1UL); i++, wbeg_plus_k++) {
              kmer <<= 2UL;
              kmer |= symbol_to_kmer(*wbeg_plus_k);
            }

            // slide //
            while (wbeg < wend) {
              kmer <<= 2UL;
              kmer |= symbol_to_kmer(*wbeg_plus_k);
              kmer_buffer.push_back(kmer);
              ++wbeg; ++wbeg_plus_k;
            }
          }
          ++rid;
          total_reads++;
        } // foreach read //
      }

      oneapi::tbb::parallel_sort(kmer_buffer.begin(), kmer_buffer.end());

      std::vector< oneapi::tbb::concurrent_vector<size_t> > kmer_freq;
      std::vector< size_t > kmer_starts;

      kmer_starts.resize( tcount );
      kmer_freq.resize( tcount );

#pragma omp parallel
      {
        size_t my_id = omp_get_thread_num(); 
        size_t data_size = kmer_buffer.size(); 
        size_t chunk_size = data_size/tcount;

        size_t start_idx = my_id*chunk_size;
        size_t end_idx = (my_id == tcount-1UL) ? data_size :
            start_idx + chunk_size;
        auto &my_kmer_freq = kmer_freq[my_id];

        my_kmer_freq.clear();
        if (start_idx >= data_size) {
          start_idx = data_size; end_idx = data_size;
        }

        if (end_idx >= data_size) {
          end_idx = data_size;
        }

        if (start_idx > 0){
          size_t prev_idx = start_idx-1UL;
          while ( (start_idx < end_idx) &&
              (kmer_buffer[prev_idx] == kmer_buffer[start_idx])) { 
            ++start_idx;
          }
        }
        kmer_starts[my_id] = start_idx;
        size_t curr_freq = 1;
        size_t uniq_idx = start_idx;

        for (size_t i=start_idx+1; i<end_idx; i++) {
          if (kmer_buffer[uniq_idx] == kmer_buffer[i]) { curr_freq++; } else {
            my_kmer_freq.push_back(curr_freq);
            kmer_buffer[++uniq_idx] = kmer_buffer[i];
            curr_freq = 1;
          }
        }
        my_kmer_freq.push_back(curr_freq);
      }

      size_t total_uniq_kmer_count = 0;
      for (size_t i=0; i<tcount; i++) {
        total_uniq_kmer_count += kmer_freq[i].size();
      }

      printf("unique kmers=%lu\n", total_uniq_kmer_count);
      printf("concurrent_vector size=%lu\n", kmer_buffer.size());
      return (size_t) total_reads;
    }


    template<typename processor_t=default_processor_t>
    size_t concurrent_count_and_process_omp_implicit_bit_vector(
        processor_t& processor) {
      size_t tcount = 0;

#pragma omp parallel
      {
        if (!omp_get_thread_num())
          tcount = omp_get_num_threads();
      }

      const char *symbol_begin = symbol_begin_;
      const char *symbol_end = symbol_end_;

      // 0 - even 1 - odd //
      //oneapi::tbb::concurrent_vector<bool> bit_encoded_reads[2UL];
      //oneapi::tbb::concurrent_vector<size_t> bit_ordered_permutation;

      std::vector<bool> bit_encoded_reads[2UL];
      std::vector<size_t> bit_ordered_permutation;
      bit_kmer_strncmp_t< std::vector<bool> > bit_kmer_strncmp(
          bit_encoded_reads[0], bit_encoded_reads[1], klen_);

      //////////////////////////////////////////////////////////////////////////
      std::vector<std::pair<const char *, const char *> > 
          threads_begin_end(tcount); 
      std::vector<size_t> threads_part_count(tcount, 0UL);
      std::vector<size_t> threads_read_count(tcount, 0UL);

#pragma omp parallel
      {
        size_t my_id = omp_get_thread_num(); 
        size_t &my_reads = threads_read_count[my_id];

        size_t &my_part_count = threads_part_count[my_id];
        std::pair<const char *, const char *> &my_begin_end = 
          threads_begin_end[my_id];

        size_t data_size = ((size_t) symbol_end_ )-((size_t)symbol_begin_);
        size_t chunk_size = data_size/tcount;

        const char *symbol_begin = symbol_begin_ + (my_id*chunk_size);
        const char *symbol_end = symbol_begin + chunk_size;

        if (symbol_begin >= symbol_end_) { 
          symbol_begin = symbol_end_; 
          symbol_end = symbol_end_;
        }
        if (symbol_end >= symbol_end_) {
          symbol_end = symbol_end_;
        }

        // adjust symbol_end //
        symbol_end = skip_till_symbol(symbol_end, symbol_end_, skip_symbol_);
        ++symbol_end;

        if (symbol_end >= symbol_end_) {
          symbol_end = symbol_end_;
        }


        my_begin_end.first = symbol_begin;
        my_begin_end.second = symbol_end;
        size_t rid = 0;
        while (symbol_begin < symbol_end) {
          // STEP-1: set the start of the 
          symbol_begin =
              skip_till_symbol(symbol_begin, symbol_end, skip_symbol_);
          if (symbol_begin >= symbol_end) { break; }
          ++symbol_begin;

          // STEP-2: skip till delimiter //  
          symbol_begin =
              skip_till_symbol(symbol_begin, symbol_end, delimiter_);
          if (symbol_begin >= symbol_end) { break;}
          ++symbol_begin;

          size_t read_len = read_len_;
          const char *read_beg = symbol_begin;
          const char *read_end = symbol_begin+read_len;

          if (read_end > symbol_end_) {
            read_end = symbol_end_; 
            read_len = symbol_end_-read_beg;
          }
          my_part_count += (read_len-klen_+1UL);
          my_reads += read_len;
        }
      }

      size_t cum_total = 0, cum_total_reads = 0;
      size_t value = threads_part_count[0];
      size_t value_reads = threads_read_count[0];

      threads_part_count[0] = cum_total;
      threads_read_count[0] = cum_total_reads;

      cum_total += value;
      cum_total_reads += value_reads;

      for (size_t i=1; i<tcount; i++) {
        value = threads_part_count[i];
        value_reads = threads_read_count[i];

        threads_part_count[i] = cum_total;
        threads_read_count[i] = cum_total_reads;

        cum_total += value;
        cum_total_reads += value_reads;
      }

      // pre-allocate vectors to sort //
      bit_ordered_permutation.resize( cum_total,  size_t(0));
      bit_encoded_reads[0].resize( cum_total_reads, false);
      bit_encoded_reads[1].resize( cum_total_reads, false);

      printf("rcount=%lu ordered_permutation.size() = %lu\n",
          (size_t) cum_total_reads , bit_ordered_permutation.size());
      std::atomic<size_t> write_count = {0UL};

#pragma omp parallel
      {
        size_t my_id = omp_get_thread_num();
        size_t my_start_idx = threads_part_count[my_id];
        size_t my_read_start_idx = threads_read_count[my_id];

        const char *symbol_begin = threads_begin_end[my_id].first;
        const char *symbol_end = threads_begin_end[my_id].second;

        if ( !((symbol_begin >= symbol_begin_) && (symbol_end <=symbol_end_)) )
        {
          throw std::logic_error("[Invalid symbol begin and end pointers]");
        }

        while (symbol_begin < symbol_end) {
          symbol_begin =
              skip_till_symbol(symbol_begin, symbol_end, skip_symbol_);
          if (symbol_begin >= symbol_end) { break; }
          ++symbol_begin;

          // STEP-2: skip till delimiter //  
          symbol_begin =
              skip_till_symbol(symbol_begin, symbol_end, delimiter_);
          if (symbol_begin >= symbol_end) { break;}
          ++symbol_begin;

          size_t read_len = read_len_;
          const char *read_beg = symbol_begin;
          const char *read_end = symbol_begin+read_len;

          if (read_end > symbol_end_) {
            read_end = symbol_end_; 
            read_len = symbol_end_-read_beg;
          }

          size_t i=0;
          for (; i<(read_len-klen_+1); i++) {
            bit_ordered_permutation[my_start_idx++] = my_read_start_idx;

            set_symbol_bits_in_container(read_beg[i], my_read_start_idx,
                  bit_encoded_reads[0UL], 
                  bit_encoded_reads[1UL] );
            ++my_read_start_idx;
          }

          for (; i<read_len; i++) {
            set_symbol_bits_in_container(read_beg[i], my_read_start_idx,
                bit_encoded_reads[0UL], bit_encoded_reads[1UL] );
            ++my_read_start_idx;
          }
          write_count += (read_len-klen_+1UL);
        }
      }
      //////////////////////////////////////////////////////////////////////////
      {
        intel::esc::stop_clock_t clock("[Sorting Time]:");
        oneapi::tbb::parallel_sort(bit_ordered_permutation, bit_kmer_strncmp);
      }

      std::vector< oneapi::tbb::concurrent_vector<size_t> > kmer_freq;
      std::vector< size_t > kmer_starts;

      kmer_starts.resize( tcount );
      kmer_freq.resize( tcount );

#pragma omp parallel
      {
        size_t my_id = omp_get_thread_num(); 
        size_t data_size = bit_ordered_permutation.size(); 
        size_t chunk_size = data_size/tcount;

        size_t start_idx = my_id*chunk_size;
        size_t end_idx = (my_id == tcount-1UL) ? data_size :
            start_idx + chunk_size;
        auto &my_kmer_freq = kmer_freq[my_id];

        my_kmer_freq.clear();
        if (start_idx >= data_size) {
          start_idx = data_size; end_idx = data_size;
        }

        if (end_idx >= data_size) {
          end_idx = data_size;
        }

        if (start_idx > 0){

          size_t prev_idx = start_idx-1UL;
          while ( (start_idx < end_idx) &&
                bit_kmer_strncmp.equiv( bit_ordered_permutation[prev_idx],
                    bit_ordered_permutation[start_idx ] ) ) 
          { ++start_idx; }
        }

        kmer_starts[my_id] = start_idx;
        size_t curr_freq = 1;
        size_t uniq_idx = start_idx;

        for (size_t i=start_idx+1; i<end_idx; i++) {
          if (bit_kmer_strncmp.equiv(bit_ordered_permutation[uniq_idx],
                  bit_ordered_permutation[i])) { 
            curr_freq++; 
          } else {
            my_kmer_freq.push_back(curr_freq);
            bit_ordered_permutation[++uniq_idx] = bit_ordered_permutation[i];
            curr_freq = 1;
          }
        }
        my_kmer_freq.push_back(curr_freq);
      }

      size_t total_uniq_kmer_count = 0;
      for (size_t i=0; i<tcount; i++) {
        total_uniq_kmer_count += kmer_freq[i].size();
      }

      printf("unique kmers=%lu\n", total_uniq_kmer_count);
      printf("concurrent_vector size=%lu\n", bit_ordered_permutation.size());
      return  cum_total_reads;
    }

    template<typename processor_t=default_processor_t>
    size_t concurrent_count_and_process_omp_implicit_v0(
				processor_t& processor) {
      size_t tcount = 0;
#pragma omp parallel
      {
        if (!omp_get_thread_num())
          tcount = omp_get_num_threads();
      }
      ////////////////////////////////////////////////////////
      // STEP-0: find the number of bytes in all the reads //

      const char *symbol_begin = symbol_begin_;
      const char *symbol_end = symbol_end_;

      kmer_strncmp_t kmer_strncmp(klen_);
      oneapi::tbb::concurrent_vector<const char *> ordered_permutation;
      ordered_permutation.reserve(1000000);
      size_t read_count =0;
      while (symbol_begin < symbol_end) {
        // STEP-1: set the start of the 
        symbol_begin =
            skip_till_symbol(symbol_begin, symbol_end, skip_symbol_);
        if (symbol_begin >= symbol_end) { break; }
        ++symbol_begin;

        // STEP-2: skip till delimiter //  
        symbol_begin =
            skip_till_symbol(symbol_begin, symbol_end, delimiter_);
        if (symbol_begin >= symbol_end) { break;}
        ++symbol_begin;

        size_t read_len = read_len_;
        const char *read_beg = symbol_begin;
        const char *read_end = symbol_begin+read_len;

        if (read_end > symbol_end_) {
          read_end = symbol_end_; 
          read_len = symbol_end_-read_beg;
        }

        if (read_len != read_len_) {
          throw std::logic_error("Variable size read\n");
        }

        for (size_t i=0; i<(read_len-klen_+1); i++) {
          ordered_permutation.push_back(read_beg + i);
        }
        ++read_count;
      }

      oneapi::tbb::parallel_sort(ordered_permutation, kmer_strncmp);

      std::vector< oneapi::tbb::concurrent_vector<size_t> > kmer_freq;
      std::vector< size_t > kmer_starts;

      kmer_starts.resize( tcount );
      kmer_freq.resize( tcount );

#pragma omp parallel
      {
        size_t my_id = omp_get_thread_num(); 
        size_t data_size = ordered_permutation.size(); 
        size_t chunk_size = data_size/tcount;

        size_t start_idx = my_id*chunk_size;
        size_t end_idx = (my_id == tcount-1UL) ? data_size :
            start_idx + chunk_size;
        auto &my_kmer_freq = kmer_freq[my_id];

        my_kmer_freq.clear();
        if (start_idx >= data_size) {
          start_idx = data_size; end_idx = data_size;
        }

        if (end_idx >= data_size) {
          end_idx = data_size;
        }

        if (start_idx > 0){
          size_t prev_idx = start_idx-1UL;
          while ( (start_idx < end_idx) &&
                !strncmp(ordered_permutation[prev_idx], 
                  ordered_permutation[start_idx], klen_) ) {
            ++start_idx;
          }
        }
        kmer_starts[my_id] = start_idx;
        size_t curr_freq = 1;
        size_t uniq_idx = start_idx;

        for (size_t i=start_idx+1; i<end_idx; i++) {
          if (strncmp(ordered_permutation[uniq_idx], ordered_permutation[i],
              klen_) == 0) { curr_freq++; } else {
            my_kmer_freq.push_back(curr_freq);
            ordered_permutation[++uniq_idx] = ordered_permutation[i];
            curr_freq = 1;
          }
        }
        my_kmer_freq.push_back(curr_freq);
      //  printf("t-%lu start=%lu end=%lu uniq=%lu total=%lu\n", 
       //    my_id, start_idx, end_idx, my_kmer_freq.size(), end_idx-start_idx);
      }

      size_t total_uniq_kmer_count = 0;
      for (size_t i=0; i<tcount; i++) {
        total_uniq_kmer_count += kmer_freq[i].size();
      }

      printf("unique kmers=%lu\n", total_uniq_kmer_count);
      printf("concurrent_vector size=%lu\n", ordered_permutation.size());
      return read_count;
    }



    template<typename processor_t=default_processor_t>
    size_t concurrent_count_and_process_omp_implicit_v0_bad(
        processor_t& processor) {
      size_t tcount = 0;
#pragma omp parallel
      {
        if (!omp_get_thread_num())
          tcount = omp_get_num_threads();
      }

      kmer_strncmp_t kmer_strncmp(klen_);
      oneapi::tbb::concurrent_vector<const char *> ordered_permutation;
      ordered_permutation.reserve(1000000);

      std::atomic<size_t> read_count = {0UL};
#pragma omp parallel
      {
        size_t my_id = omp_get_thread_num(); 
        size_t data_size = ((size_t) symbol_end_ )-((size_t)symbol_begin_);
        size_t chunk_size = data_size/tcount;

        const char *symbol_begin = symbol_begin_ + (my_id*chunk_size);
        const char *symbol_end = symbol_begin + chunk_size;

        if (symbol_begin >= symbol_end_) { 
          symbol_begin = symbol_end_; 
          symbol_end = symbol_end_;
        }
        if (symbol_end >= symbol_end_) {
          symbol_end = symbol_end_;
        }

        // adjust symbol_end //
        symbol_end = skip_till_symbol(symbol_end, symbol_end_, skip_symbol_);
        ++symbol_end;

        if (symbol_end >= symbol_end_) {
          symbol_end = symbol_end_;
        }


        size_t rid = 0;
        while (symbol_begin < symbol_end) {
          // STEP-1: set the start of the 
          symbol_begin =
              skip_till_symbol(symbol_begin, symbol_end, skip_symbol_);
          if (symbol_begin >= symbol_end) { break; }
          ++symbol_begin;

          // STEP-2: skip till delimiter //  
          symbol_begin =
              skip_till_symbol(symbol_begin, symbol_end, delimiter_);
          if (symbol_begin >= symbol_end) { break;}
          ++symbol_begin;

          size_t read_len = read_len_;
          const char *read_beg = symbol_begin;
          const char *read_end = symbol_begin+read_len;

          if (read_end > symbol_end_) {
            read_end = symbol_end_; 
            read_len = symbol_end_-read_beg;
          }
          for (size_t i=0; i<(read_len-klen_+1); i++) {
            ordered_permutation.push_back(read_beg + i);
          }
        }
        ++read_count;
      }

      oneapi::tbb::parallel_sort(ordered_permutation, kmer_strncmp);

      std::vector< oneapi::tbb::concurrent_vector<size_t> > kmer_freq;
      std::vector< size_t > kmer_starts;

      kmer_starts.resize( tcount );
      kmer_freq.resize( tcount );

#pragma omp parallel
      {
        size_t my_id = omp_get_thread_num(); 
        size_t data_size = ordered_permutation.size(); 
        size_t chunk_size = data_size/tcount;

        size_t start_idx = my_id*chunk_size;
        size_t end_idx = (my_id == tcount-1UL) ? data_size :
            start_idx + chunk_size;
        auto &my_kmer_freq = kmer_freq[my_id];

        my_kmer_freq.clear();
        if (start_idx >= data_size) {
          start_idx = data_size; end_idx = data_size;
        }

        if (end_idx >= data_size) {
          end_idx = data_size;
        }

        if (start_idx > 0){
          size_t prev_idx = start_idx-1UL;
          while ( (start_idx < end_idx) &&
                !strncmp(ordered_permutation[prev_idx], 
                  ordered_permutation[start_idx], klen_) ) {
            ++start_idx;
          }
        }
        kmer_starts[my_id] = start_idx;
        size_t curr_freq = 1;
        size_t uniq_idx = start_idx;

        for (size_t i=start_idx+1; i<end_idx; i++) {
          if (strncmp(ordered_permutation[uniq_idx], ordered_permutation[i],
              klen_) == 0) { curr_freq++; } else {
            my_kmer_freq.push_back(curr_freq);
            ordered_permutation[++uniq_idx] = ordered_permutation[i];
            curr_freq = 1;
          }
        }
        my_kmer_freq.push_back(curr_freq);
      }

      size_t total_uniq_kmer_count = 0;
      for (size_t i=0; i<tcount; i++) {
        total_uniq_kmer_count += kmer_freq[i].size();
      }

      printf("unique kmers=%lu\n", total_uniq_kmer_count);
      printf("concurrent_vector size=%lu\n", ordered_permutation.size());
      return (size_t) read_count;
    }


    template<typename processor_t=default_processor_t>
    size_t concurrent_count_and_process_omp_implicit_v2_bad(
        processor_t& processor) {
      std::atomic<size_t> read_count = {0UL};

      size_t tcount = 0;
#pragma omp parallel
      {
        if (!omp_get_thread_num())
          tcount = omp_get_num_threads();
      }

      kmer_strncmp_t kmer_strncmp(klen_);
      std::vector< std::pair<const char *, int> >
          ordered_permutation;

      std::atomic<size_t> ordered_permutation_length = {0UL};
      std::atomic<size_t> write_count = {0UL};

      std::vector<std::pair<const char *, const char *> > 
          threads_begin_end(tcount); 
      std::vector<size_t> threads_part_count(tcount, 0UL);

      typedef std::unordered_map<const char*, int,
              implicit_kmer_hash_t, implicit_kmer_equal_t> local_hash_t; 

      implicit_kmer_hash_t implicit_kmer_hash(klen_);
      implicit_kmer_equal_t implicit_kmer_equal(klen_);

      local_hash_t local_hash_template(100000UL, implicit_kmer_hash,
            implicit_kmer_equal);
      std::vector<local_hash_t> local_hashes(tcount, local_hash_template);

      {
        intel::esc::stop_clock_t pre_sort_clock("[Scan+PreSort]");
#pragma omp parallel
        {
          size_t my_id = omp_get_thread_num(); 
          local_hash_t& my_local_hash = local_hashes[my_id];

          size_t &my_part_count = threads_part_count[my_id];
          std::pair<const char *, const char *> &my_begin_end = 
            threads_begin_end[my_id];

          size_t data_size = ((size_t) symbol_end_ )-((size_t)symbol_begin_);
          size_t chunk_size = data_size/tcount;

          const char *symbol_begin = symbol_begin_ + (my_id*chunk_size);
          const char *symbol_end = symbol_begin + chunk_size;

          if (symbol_begin >= symbol_end_) { 
            symbol_begin = symbol_end_; 
            symbol_end = symbol_end_;
          }
          if (symbol_end >= symbol_end_) {
            symbol_end = symbol_end_;
          }

          // adjust symbol_end //
          symbol_end = skip_till_symbol(symbol_end, symbol_end_, skip_symbol_);
          ++symbol_end;

          if (symbol_end >= symbol_end_) {
            symbol_end = symbol_end_;
          }


          my_begin_end.first = symbol_begin;
          my_begin_end.second = symbol_end;
          size_t rid = 0;
          while (symbol_begin < symbol_end) {
            // STEP-1: set the start of the 
            symbol_begin =
                skip_till_symbol(symbol_begin, symbol_end, skip_symbol_);
            if (symbol_begin >= symbol_end) { break; }
            ++symbol_begin;

            // STEP-2: skip till delimiter //  
            symbol_begin =
                skip_till_symbol(symbol_begin, symbol_end, delimiter_);
            if (symbol_begin >= symbol_end) { break;}
            ++symbol_begin;

            size_t read_len = read_len_;
            const char *read_beg = symbol_begin;
            const char *read_end = symbol_begin+read_len;

            if (read_end > symbol_end_) {
              read_end = symbol_end_; 
              read_len = symbol_end_-read_beg;
            }

            for (size_t i=0; i<(read_len-klen_+1UL); i++) {
              my_local_hash[ (read_beg + i) ]++;
            }
            ++read_count;
          }
          my_part_count = my_local_hash.size();
        }

        size_t cum_total = 0;
        size_t value = threads_part_count[0];

        threads_part_count[0] = cum_total;
        cum_total += value;
        for (size_t i=1; i<tcount; i++) {
          value = threads_part_count[i];
          threads_part_count[i] = cum_total;
          cum_total += value;
        }

        ordered_permutation.resize( cum_total);
        printf("rcount=%lu ordered_permutation.size() = %lu\n",
            (size_t) read_count, ordered_permutation.size());


#pragma omp parallel
        {
          size_t my_id = omp_get_thread_num();
          size_t my_start_idx = threads_part_count[my_id];
          local_hash_t& my_local_hash = local_hashes[my_id];

          for (auto itr=my_local_hash.begin(); itr!=my_local_hash.end(); ++itr) {
            ordered_permutation[my_start_idx++] = *itr;
          }
          my_local_hash.clear();
        }

      }


      printf("...starting to sort write_count=%lu size=%lu\n",
          (size_t) write_count, ordered_permutation.size());
      {
        intel::esc::stop_clock_t clock("[parallel sort time]:");
        oneapi::tbb::parallel_sort(ordered_permutation, kmer_strncmp);
      }
      printf("..done sorting...\n");

      std::vector< size_t > kmer_starts;
      std::vector< size_t > uniq_kmer_counts;
      kmer_starts.resize( tcount );
      uniq_kmer_counts.resize( tcount );

      {
        intel::esc::stop_clock_t merge_clock("[Merge Clock]:");
#pragma omp parallel
        {
          size_t my_id = omp_get_thread_num(); 
          size_t data_size = ordered_permutation.size(); 
          size_t chunk_size = data_size/tcount;
          size_t &my_uniq_kmers = uniq_kmer_counts[my_id];
          my_uniq_kmers = 0UL;

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
            size_t prev_idx = start_idx-1UL;
            while ( (start_idx < end_idx) &&
                  !strncmp(ordered_permutation[prev_idx].first, 
                    ordered_permutation[start_idx].first, klen_) ) {
              ++start_idx;
            }
          }
          kmer_starts[my_id] = start_idx;
          size_t uniq_idx = start_idx;
          int curr_freq = ordered_permutation[uniq_idx].second;
          ++my_uniq_kmers;

          for (size_t i=start_idx+1; i<end_idx; i++) {
            if (strncmp(ordered_permutation[uniq_idx].first, 
                  ordered_permutation[i].first, klen_) == 0) { 
              curr_freq += ordered_permutation[i].second; 
            } else {
              ordered_permutation[uniq_idx].second = curr_freq;
              ordered_permutation[++uniq_idx] = ordered_permutation[i];
              ++my_uniq_kmers;
              curr_freq = ordered_permutation[uniq_idx].second;
            }
          }
          ordered_permutation[uniq_idx].second = curr_freq;
        } // pragma omp parallel //

        size_t total_uniq_kmer_count = 0;
        for (size_t i=0; i<tcount; i++) {
          total_uniq_kmer_count += uniq_kmer_counts[i]; 
        }

        printf("unique kmers=%lu\n", total_uniq_kmer_count);
        printf("concurrent_vector size=%lu\n", ordered_permutation.size());
      }
      return (size_t) read_count;
    }


    template<typename processor_t=default_processor_t>
    size_t concurrent_count_and_process_omp_implicit_v1(
        processor_t& processor) {

      size_t tcount = 0;
#pragma omp parallel
      {
        if (!omp_get_thread_num())
          tcount = omp_get_num_threads();
      }

      kmer_strncmp_t kmer_strncmp(klen_);
      std::vector< std::pair<const char *, size_t> > ordered_permutation;

      std::atomic<size_t> read_count = {0UL};
      std::atomic<size_t> ordered_permutation_length = {0UL};

      std::vector<std::pair<const char *, const char *> > 
          threads_begin_end(tcount); 
      std::vector<size_t> threads_part_count(tcount, 0UL);

#pragma omp parallel
      {
        size_t my_id = omp_get_thread_num(); 

        size_t &my_part_count = threads_part_count[my_id];
        std::pair<const char *, const char *> &my_begin_end = 
          threads_begin_end[my_id];

        size_t data_size = ((size_t) symbol_end_ )-((size_t)symbol_begin_);
        size_t chunk_size = data_size/tcount;

        const char *symbol_begin = symbol_begin_ + (my_id*chunk_size);
        const char *symbol_end = symbol_begin + chunk_size;

        if (symbol_begin >= symbol_end_) { 
          symbol_begin = symbol_end_; 
          symbol_end = symbol_end_;
        }
        if (symbol_end >= symbol_end_) {
          symbol_end = symbol_end_;
        }

        // adjust symbol_end //
        symbol_end = skip_till_symbol(symbol_end, symbol_end_, skip_symbol_);
        ++symbol_end;

        if (symbol_end >= symbol_end_) {
          symbol_end = symbol_end_;
        }


        my_begin_end.first = symbol_begin;
        my_begin_end.second = symbol_end;
        size_t rid = 0;
        while (symbol_begin < symbol_end) {
          // STEP-1: set the start of the 
          symbol_begin =
              skip_till_symbol(symbol_begin, symbol_end, skip_symbol_);
          if (symbol_begin >= symbol_end) { break; }
          ++symbol_begin;

          // STEP-2: skip till delimiter //  
          symbol_begin =
              skip_till_symbol(symbol_begin, symbol_end, delimiter_);
          if (symbol_begin >= symbol_end) { break;}
          ++symbol_begin;

          size_t read_len = read_len_;
          const char *read_beg = symbol_begin;
          const char *read_end = symbol_begin+read_len;

          if (read_end > symbol_end_) {
            read_end = symbol_end_; 
            read_len = symbol_end_-read_beg;
          }
          my_part_count += (read_len-klen_+1UL);
          ++read_count;
        }
      }

      size_t cum_total = 0;
      size_t value = threads_part_count[0];

      threads_part_count[0] = cum_total;
      cum_total += value;
      for (size_t i=1; i<tcount; i++) {
        value = threads_part_count[i];
        threads_part_count[i] = cum_total;
        cum_total += value;
      }

      ordered_permutation.resize( cum_total);
      printf("rcount=%lu ordered_permutation.size() = %lu\n",
          (size_t) read_count, ordered_permutation.size());

      std::atomic<size_t> write_count = {0UL};

#pragma omp parallel
      {
        size_t my_id = omp_get_thread_num();
        size_t my_start_idx = threads_part_count[my_id];
        const char *symbol_begin = threads_begin_end[my_id].first;
        const char *symbol_end = threads_begin_end[my_id].second;

        if ( !((symbol_begin >= symbol_begin_) && (symbol_end <=symbol_end_)) )
        {
          throw std::logic_error("[Invalid symbol begin and end pointers]");
        }

        while (symbol_begin < symbol_end) {
          symbol_begin =
              skip_till_symbol(symbol_begin, symbol_end, skip_symbol_);
          if (symbol_begin >= symbol_end) { break; }
          ++symbol_begin;

          // STEP-2: skip till delimiter //  
          symbol_begin =
              skip_till_symbol(symbol_begin, symbol_end, delimiter_);
          if (symbol_begin >= symbol_end) { break;}
          ++symbol_begin;

          size_t read_len = read_len_;
          const char *read_beg = symbol_begin;
          const char *read_end = symbol_begin+read_len;

          if (read_end > symbol_end_) {
            read_end = symbol_end_; 
            read_len = symbol_end_-read_beg;
          }

          for (size_t i=0; i<(read_len-klen_+1UL); i++) {
            if (my_start_idx >= ordered_permutation.size()) {
              throw std::logic_error("[Invalid Bounds]: error");
            }
            ordered_permutation[my_start_idx++] = 
              std::make_pair((read_beg + i), 1UL);
          }

          write_count += (read_len-klen_+1UL);
        }
      }


      printf("...starting to sort write_count=%lu size=%lu\n",
          (size_t) write_count, ordered_permutation.size());
      {
        intel::esc::stop_clock_t clock("[parallel sort time]:");
        oneapi::tbb::parallel_sort(ordered_permutation, kmer_strncmp);
      }
      printf("..done sorting...\n");

      std::vector< size_t > kmer_starts;
      std::vector< size_t > uniq_kmer_counts;
      kmer_starts.resize( tcount );
      uniq_kmer_counts.resize( tcount );

      {
        intel::esc::stop_clock_t merge_time("[Merge Time]:");
#pragma omp parallel
        {
          size_t my_id = omp_get_thread_num(); 
          size_t data_size = ordered_permutation.size(); 
          size_t chunk_size = data_size/tcount;
          size_t &my_uniq_kmers = uniq_kmer_counts[my_id];
          my_uniq_kmers = 0UL;

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
            size_t prev_idx = start_idx-1UL;
            while ( (start_idx < end_idx) &&
                  !strncmp(ordered_permutation[prev_idx].first, 
                    ordered_permutation[start_idx].first, klen_) ) {
              ++start_idx;
            }
          }
          kmer_starts[my_id] = start_idx;
          size_t curr_freq = 1;
          size_t uniq_idx = start_idx;
          ++my_uniq_kmers;

          for (size_t i=start_idx+1; i<end_idx; i++) {
            if (strncmp(ordered_permutation[uniq_idx].first, 
                  ordered_permutation[i].first, klen_) == 0) { 
              curr_freq++; 
            } else {
              ordered_permutation[uniq_idx].second = curr_freq;
              ordered_permutation[++uniq_idx] = ordered_permutation[i];
              ++my_uniq_kmers;
              curr_freq = 1;
            }
          }
          ordered_permutation[uniq_idx].second = curr_freq;
        }

        size_t total_uniq_kmer_count = 0;
        for (size_t i=0; i<tcount; i++) {
          total_uniq_kmer_count += uniq_kmer_counts[i]; 
        }

        printf("unique kmers=%lu\n", total_uniq_kmer_count);
        printf("concurrent_vector size=%lu\n", ordered_permutation.size());
      }
      return (size_t) read_count;
    }


    template<typename processor_t=default_processor_t>
    size_t concurrent_count_and_process_omp(processor_t& processor) {
      return concurrent_count_and_process_omp_implicit_bit_vector(processor);
    }

    template<typename processor_t=default_processor_t>
    size_t concurrent_count_and_process_omp_implicit(processor_t& processor) {
      size_t tcount = 0;
#pragma omp parallel
      {
        if (!omp_get_thread_num())
          tcount = omp_get_num_threads();
      }
      ////////////////////////////////////////////////////////
      // STEP-0: find the number of bytes in all the reads //


      kmer_strncmp_t kmer_strncmp(klen_);
      std::vector<const char *> ordered_permutation;

      size_t read_count = 0UL;
      size_t kmer_count = 0UL;

      {
        const char *symbol_begin = symbol_begin_;
        const char *symbol_end = symbol_end_;
        while (symbol_begin < symbol_end) {
          // STEP-1: set the start of the 
          symbol_begin =
              skip_till_symbol(symbol_begin, symbol_end, skip_symbol_);
          if (symbol_begin >= symbol_end) { break; }
          ++symbol_begin;

          // STEP-2: skip till delimiter //  
          symbol_begin =
              skip_till_symbol(symbol_begin, symbol_end, delimiter_);
          if (symbol_begin >= symbol_end) { break;}
          ++symbol_begin;

          size_t read_len = read_len_;
          const char *read_beg = symbol_begin;
          const char *read_end = symbol_begin+read_len;

          if (read_end > symbol_end_) {
            read_end = symbol_end_; 
            read_len = symbol_end_-read_beg;
          }

          if (read_len != read_len_) {
            throw std::logic_error("Variable size read\n");
          }

          kmer_count += (read_len-klen_+1UL);
          ++read_count;
        }
      }

      ordered_permutation.resize(kmer_count, NULL);

      {
        size_t start_idx = 0;
        const char *symbol_begin = symbol_begin_;
        const char *symbol_end = symbol_end_;
        while (symbol_begin < symbol_end) {
          // STEP-1: set the start of the 
          symbol_begin =
              skip_till_symbol(symbol_begin, symbol_end, skip_symbol_);
          if (symbol_begin >= symbol_end) { break; }
          ++symbol_begin;

          // STEP-2: skip till delimiter //  
          symbol_begin =
              skip_till_symbol(symbol_begin, symbol_end, delimiter_);
          if (symbol_begin >= symbol_end) { break;}
          ++symbol_begin;

          size_t read_len = read_len_;
          const char *read_beg = symbol_begin;
          const char *read_end = symbol_begin+read_len;

          if (read_end > symbol_end_) {
            read_end = symbol_end_; 
            read_len = symbol_end_-read_beg;
          }

          if (read_len != read_len_) {
            throw std::logic_error("Variable size read\n");
          }

          for (size_t i=0; i<(read_len-klen_+1UL); i++) {
            ordered_permutation[ start_idx++ ] = read_beg + i;
          }
        }
      }

      {
        intel::esc::stop_clock_t clock("[Parallel Sort Time]:");
        oneapi::tbb::parallel_sort(ordered_permutation, kmer_strncmp);
      }

      std::vector< oneapi::tbb::concurrent_vector<size_t> > kmer_freq;
      std::vector< size_t > kmer_starts;

      kmer_starts.resize( tcount );
      kmer_freq.resize( tcount );

#pragma omp parallel
      {
        size_t my_id = omp_get_thread_num(); 
        size_t data_size = ordered_permutation.size(); 
        size_t chunk_size = data_size/tcount;

        size_t start_idx = my_id*chunk_size;
        size_t end_idx = (my_id == tcount-1UL) ? data_size :
            start_idx + chunk_size;
        auto &my_kmer_freq = kmer_freq[my_id];

        my_kmer_freq.clear();
        if (start_idx >= data_size) {
          start_idx = data_size; end_idx = data_size;
        }

        if (end_idx >= data_size) {
          end_idx = data_size;
        }

        if (start_idx > 0){
          size_t prev_idx = start_idx-1UL;
          while ( (start_idx < end_idx) &&
                !strncmp(ordered_permutation[prev_idx], 
                  ordered_permutation[start_idx], klen_) ) {
            ++start_idx;
          }
        }
        kmer_starts[my_id] = start_idx;
        size_t curr_freq = 1;
        size_t uniq_idx = start_idx;

        for (size_t i=start_idx+1; i<end_idx; i++) {
          if (strncmp(ordered_permutation[uniq_idx], ordered_permutation[i],
              klen_) == 0) { curr_freq++; } else {
            my_kmer_freq.push_back(curr_freq);
            ordered_permutation[++uniq_idx] = ordered_permutation[i];
            curr_freq = 1;
          }
        }
        my_kmer_freq.push_back(curr_freq);
      }

      size_t total_uniq_kmer_count = 0;
      for (size_t i=0; i<tcount; i++) {
        total_uniq_kmer_count += kmer_freq[i].size();
      }

      printf("unique kmers=%lu\n", total_uniq_kmer_count);
      printf("concurrent_vector size=%lu\n", ordered_permutation.size());
      return read_count;
    }

    template<typename processor_t=default_processor_t>
    size_t concurrent_count_and_process_omp_good(processor_t& processor) {
      ////////////////////////////////////////////////////////
      // STEP-0: find the number of bytes in all the reads //

      const char *symbol_begin = symbol_begin_;
      const char *symbol_end = symbol_end_;

      size_t data_size = ((size_t) symbol_end )-((size_t)symbol_begin);
      std::vector<char> read_data;

      read_data.resize( data_size, '\0' );
      char *read_buffer = read_data.data();
      char *curr_read_buffer_head = read_buffer;
      const char *read_buffer_end = read_data.data() + read_data.size();

      size_t read_count =0;
      while (symbol_begin < symbol_end) {
        // STEP-1: set the start of the 
        symbol_begin =
            skip_till_symbol(symbol_begin, symbol_end, skip_symbol_);
        if (symbol_begin >= symbol_end) { break; }
        ++symbol_begin;

        // STEP-2: skip till delimiter //  
        symbol_begin =
            skip_till_symbol(symbol_begin, symbol_end, delimiter_);
        if (symbol_begin >= symbol_end) { break;}
        ++symbol_begin;

        size_t read_len = read_len_;
        const char *read_beg = symbol_begin;
        const char *read_end = symbol_begin+read_len;

        if (read_end > symbol_end_) {
          read_end = symbol_end_; 
          read_len = symbol_end_-read_beg;
        }

        if (read_len != read_len_) {
          throw std::logic_error("Variable size read\n");
        }

        if (curr_read_buffer_head >= read_buffer_end) {
          throw std::logic_error("Invalid buffer size");
        }
        memcpy(curr_read_buffer_head, read_beg, read_len);
        curr_read_buffer_head += read_len;
        read_count++;
      }


      oneapi::tbb::concurrent_vector<unsigned int> ordered_permutation;
      size_t kcol_size = (read_len_ - klen_ + 1UL); 
      size_t kmer_count = read_count*kcol_size;
      ordered_permutation.resize(kmer_count);

      implicit_kmer_compare_t kmer_compare(read_buffer, kcol_size, read_len_,
            klen_);

#pragma omp parallel for
      for (size_t i=0; i<kmer_count; i++) { ordered_permutation[i] = i; }


      oneapi::tbb::parallel_sort(ordered_permutation, kmer_compare);
      printf("...done sorting...\n");

      std::vector<size_t> kmer_freq;
      kmer_freq.reserve( (size_t) (double(0.05)*double(kmer_count)) );
      size_t curr_freq = 1;
      size_t uniq_idx = 0;

      for (size_t i=1; i<kmer_count; i++) {
        if (kmer_compare.equiv( ordered_permutation[uniq_idx],
                ordered_permutation[i])) {
            curr_freq++;
        } else {
          kmer_freq.push_back(curr_freq);
          ordered_permutation[++uniq_idx] = ordered_permutation[i];
          curr_freq = 1;
        }
      }
      kmer_freq.push_back(curr_freq);

      printf("unique kmers=%lu\n", kmer_freq.size());
      printf("concurrent_vector size=%lu\n", ordered_permutation.size());
      return 0;

    }


    template<typename processor_t=default_processor_t>
    size_t concurrent_count_and_process_omp_v3(processor_t& processor) {

      concurrent_hash_map_t concurrent_hash_map(1000000);
      std::atomic<size_t> total_reads = {0UL};
#pragma omp parallel
      {
        size_t tcount = omp_get_num_threads();
        size_t my_id = omp_get_thread_num(); 
        size_t data_size = ((size_t) symbol_end_ )-((size_t)symbol_begin_);
        size_t chunk_size = data_size/tcount;

        const char *symbol_begin = symbol_begin_ + (my_id*chunk_size);
        const char *symbol_end = symbol_begin + chunk_size;

        if (symbol_begin >= symbol_end_) { 
          symbol_begin = symbol_end_; 
          symbol_end = symbol_end_;
        }
        if (symbol_end >= symbol_end_) {
          symbol_end = symbol_end_;
        }

        // adjust symbol_end //
        symbol_end = skip_till_symbol(symbol_end, symbol_end_, skip_symbol_);
        ++symbol_end;

        if (symbol_end >= symbol_end_) {
          symbol_end = symbol_end_;
        }


        size_t rid = 0;
        while (symbol_begin < symbol_end) {
          // STEP-1: set the start of the 
          symbol_begin =
              skip_till_symbol(symbol_begin, symbol_end, skip_symbol_);
          if (symbol_begin >= symbol_end) { break; }
          ++symbol_begin;

          // STEP-2: skip till delimiter //  
          symbol_begin =
              skip_till_symbol(symbol_begin, symbol_end, delimiter_);
          if (symbol_begin >= symbol_end) { break;}
          ++symbol_begin;

          size_t read_len = read_len_;
          const char *read_beg = symbol_begin;
          const char *read_end = symbol_begin+read_len;

          if (read_end > symbol_end_) {
            read_end = symbol_end_; 
            read_len = symbol_end_-read_beg;
          }
          if (read_len < klen_) { continue; } // read small //

          // adjust sliding window indicies //
          const char *window_end = read_beg + (read_len -klen_ +1UL);
          {
            //Invariant: wbeg+k <= read_end //
            const char *wbeg = read_beg; 
            const char *wend = window_end;

            // setup first k-mer //
            kmer_t kmer=kmer_t(0UL);
            const char *wbeg_plus_k = wbeg;
            for (size_t i=0; i<(klen_-1UL); i++, wbeg_plus_k++) {
              kmer <<= 2UL;
              kmer |= symbol_to_kmer(*wbeg_plus_k);
            }

            // slide //
            while (wbeg < wend) {
              kmer <<= 2UL;
              kmer |= symbol_to_kmer(*wbeg_plus_k);
              {
                accessor_t accessor;
                concurrent_hash_map.insert(accessor, kmer);
                accessor->second++;
                accessor.release();
              }
              ++wbeg; ++wbeg_plus_k;
            }
          }
          ++rid;
          total_reads++;
          if (!my_id) {
            if (((size_t)total_reads)%100000) {
              printf("processed %lu reads\r", (size_t) total_reads);
            }
          }
        } // foreach read //
        //printf("tid-%lu reads=%lu\n", my_id, rid);
      }

      // trailing case //
      processor(concurrent_hash_map.begin(), concurrent_hash_map.end());
      concurrent_hash_map.clear();
      return total_reads;
    }

    template<typename processor_t=default_processor_t>
    size_t default_count_and_process_sort(processor_t& processor) {

      std::vector<kmer_t> kmer_buffer;
      size_t rid = 0;

      while (symbol_begin_ != symbol_end_) {
        // STEP-1: set the start of the 
        symbol_begin_ =
            skip_till_symbol(symbol_begin_, symbol_end_, skip_symbol_);
        if (symbol_begin_ == symbol_end_) { break; }
        ++symbol_begin_;

        // STEP-2: skip till delimiter //  
        symbol_begin_ =
            skip_till_symbol(symbol_begin_, symbol_end_, delimiter_);
        if (symbol_begin_ == symbol_end_) { break;}
        ++symbol_begin_;

        size_t read_len = read_len_;
        const char *read_beg = symbol_begin_;
        const char *read_end = symbol_begin_+read_len;

        if (read_end > symbol_end_) {
          read_end = symbol_end_; 
          read_len = symbol_end_-read_beg;
        }
        if (read_len < klen_) { continue; } // read small //

        // adjust sliding window indicies //
        const char *window_end = read_beg + (read_len -klen_ +1UL);
        {
          //Invariant: wbeg+k <= read_end //
          const char *wbeg = read_beg; 
          const char *wend = window_end;

          // setup first k-mer //
          kmer_t kmer=kmer_t(0UL);
          const char *wbeg_plus_k = wbeg;
          for (size_t i=0; i<(klen_-1UL); i++, wbeg_plus_k++) {
            kmer <<= 2UL;
            kmer |= symbol_to_kmer(*wbeg_plus_k);
          }

          // slide //
          while (wbeg < wend) {
            kmer <<= 2UL;
            kmer |= symbol_to_kmer(*wbeg_plus_k);
            kmer_buffer.push_back(kmer);
            ++wbeg; ++wbeg_plus_k;
          }
        }

#if 0 //TODO(vamsikku): //
        if (kmer_hash_map.size() > max_kmer_count_) {
          processor(kmer_hash_map.begin(), kmer_hash_map.end());
          kmer_hash_map.clear();
        }
#endif
        ++rid;
      } // foreach read //

      std::vector< std::pair<kmer_t, size_t> > kmer_freq;
      std::sort(kmer_buffer.begin(), kmer_buffer.end());

      if (kmer_buffer.empty()) { return rid; }

      kmer_t prev = kmer_buffer[0];
      size_t freq = 1UL;

      for (size_t i=1; i<kmer_buffer.size(); i++) {
        if (prev == kmer_buffer[i]) { ++freq; }
        else {
          kmer_freq.push_back(std::make_pair(prev, freq));
          prev = kmer_buffer[i];
          freq = 1UL;
        }
      }
      kmer_buffer.clear();
      kmer_freq.push_back( std::make_pair(prev, freq));
      printf("Reads=%lu , Unique Kmers=%lu\n", rid, kmer_freq.size());
      //processor(kmer_freq.begin(), kmer_freq.end());
      return rid;
    }


    template<typename processor_t=default_processor_t>
    size_t default_count_and_process(processor_t& processor) {
      hash_map_t kmer_hash_map;
      size_t rid = 0;

      while (symbol_begin_ != symbol_end_) {
        // STEP-1: set the start of the 
        symbol_begin_ =
            skip_till_symbol(symbol_begin_, symbol_end_, skip_symbol_);
        if (symbol_begin_ == symbol_end_) { break; }
        ++symbol_begin_;

        // STEP-2: skip till delimiter //  
        symbol_begin_ =
            skip_till_symbol(symbol_begin_, symbol_end_, delimiter_);
        if (symbol_begin_ == symbol_end_) { break;}
        ++symbol_begin_;

        size_t read_len = read_len_;
        const char *read_beg = symbol_begin_;
        const char *read_end = symbol_begin_+read_len;

        if (read_end > symbol_end_) {
          read_end = symbol_end_; 
          read_len = symbol_end_-read_beg;
        }
        if (read_len < klen_) { continue; } // read small //

        // adjust sliding window indicies //
        const char *window_end = read_beg + (read_len -klen_ +1UL);
        {
          //Invariant: wbeg+k <= read_end //
          const char *wbeg = read_beg; 
          const char *wend = window_end;

          // setup first k-mer //
          kmer_t kmer=kmer_t(0UL);
          const char *wbeg_plus_k = wbeg;
          for (size_t i=0; i<(klen_-1UL); i++, wbeg_plus_k++) {
            kmer <<= 2UL;
            kmer |= symbol_to_kmer(*wbeg_plus_k);
          }

          // slide //
          while (wbeg < wend) {
            kmer <<= 2UL;
            kmer |= symbol_to_kmer(*wbeg_plus_k);
            kmer_hash_map[kmer]++;
            ++wbeg; ++wbeg_plus_k;
          }
        }

        if (kmer_hash_map.size() > max_kmer_count_) {
          processor(kmer_hash_map.begin(), kmer_hash_map.end());
          kmer_hash_map.clear();
        }
        ++rid;
      } // foreach read //

      // trailing case //
      processor(kmer_hash_map.begin(), kmer_hash_map.end());
      printf("Reads=%lu , Unique Kmers=%lu\n", rid, kmer_hash_map.size());
      kmer_hash_map.clear();
      return rid;
    }



  private:

    const char* skip_till_symbol(const char *beg, const char *end, 
        char symbol) const {
      while (beg != end) {
        if (*beg == symbol) { return beg; }
        ++beg;
      }
      return beg;
    }

    kmer_t symbol_to_kmer(char c) const {
      switch(c) {
        case 'A':
          return kmer_t(0UL);
        case 'C':
          return kmer_t(1UL);
        case 'T':
          return kmer_t(2UL);
        case 'G':
          return kmer_t(3UL);
        default:
          return kmer_t(0UL);
      }
    }

    template<typename PushBackContainer>
    void push_symbol_bits_to_container(char c,
          PushBackContainer& container_even,
          PushBackContainer& container_odd) const {

      bool zero_bit = false;
      bool one_bit = true;

      switch (c) {
        case 'A':
          // 00 //
          container_odd.push_back(zero_bit); container_even.push_back(zero_bit);
          break;

        case 'C':
          // 01 //
          container_odd.push_back(zero_bit); container_even.push_back(one_bit);
          break;

        case 'T':
          // 10 //
          container_odd.push_back(one_bit); container_even.push_back(zero_bit);
          break;

        case 'G':
          // 11 //
          container_odd.push_back(one_bit); container_even.push_back(one_bit);
          break;

        default:
          container_odd.push_back(zero_bit); container_even.push_back(zero_bit);
      }
    }

    template<typename PushBackContainer>
    inline void set_symbol_bits_in_container(char c, size_t i,
          PushBackContainer& container_even,
          PushBackContainer& container_odd) const {
      bool zero_bit = false;
      bool one_bit = true;

      switch (c) {
        case 'A':
          // 00 //
          container_odd[i]= zero_bit;
          container_even[i] = zero_bit;
          break;

        case 'C':
          // 01 //
          container_odd[i]= zero_bit;
          container_even[i] = one_bit;
          break;

        case 'T':
          // 10 //
          container_odd[i]= one_bit;
          container_even[i] = zero_bit;
          break;

        case 'G':
          // 11 //
          container_odd[i]= one_bit;
          container_even[i] = one_bit;
          break;

        default:
          container_odd[i]= zero_bit;
          container_even[i] = zero_bit;
      }
    }

    template<typename BitContainer>
    inline char get_symbol_bits_in_container(size_t i, 
        const BitContainer& container_even, 
        const BitContainer& container_odd) {
      static char table[4] = { 'A', 'C', 'T', 'G'};
      size_t v = ((size_t)( container_odd[i]))*2UL + 
          ((size_t) (container_even[i]));
      return table[v];
    }

  private:

    size_t klen_;
    char const* symbol_begin_;
    char const* symbol_end_;
    size_t read_len_;
    size_t max_kmer_count_; // max k-mers in partial_count_buffer_ //
    char delimiter_;
    char skip_symbol_;
    partial_count_buffer_t partial_count_buffer_;
    concurrent_hash_map_t *concurrent_hash_map_;
}; // class Concurrent_Kmer_Partial_Counter //


} // namespace esc //
} // namespace intel //
#endif
