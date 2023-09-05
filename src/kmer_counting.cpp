// ***********************************************************************
//
// PaKman: Algorithm for generating genomic contigs on distributed-memory machines
// 
// Priyanka Ghosh (Pacific Northwest National Laboratory)
// Sriram Krishnamoorthy (Pacific Northwest National Laboratory)
// Ananth Kalyanaraman (Washington State University)
//               
//
// ***********************************************************************
//
//       Copyright (2020) Battelle Memorial Institute
//                      All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the copyright holder nor the names of its
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
// FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
// COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
// INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
// BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
// ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// ************************************************************************

#include <atomic>
#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <inttypes.h>
#include <vector>
#include <set>
#include <algorithm>
#include <unordered_map>
#include <parallel/algorithm>
#include <numeric>
#include <omp.h>
#include "distribute_kmers.h"
#include "timers.h"

#include "concurrent_kmer_counter.hpp"
#include "asynchronous_distributed_kmer_counting.hpp"

extern bool USE_NEW_CONCURRENT_KMER_COUNTER;
extern long int MAX_KMER_COUNT;
extern int rank, size;
extern int coverage;
extern int num_threads;
extern int num_batch_transfers;
uint64_t num_recalculate_lmer, global_num_recalculate_lmer;

extern std::vector<lmer_t> lmer_frequency;
extern std::vector<lmer_t> global_lmer_frequency;

extern std::vector<KmerPairs> kmer_proc_buf;


void sort_recv_buffer(std::vector<KmerPairs>& kmer_recv_buf, 
                      std::vector<int>& rcounts_kmer, 
                      std::vector<int>& rdisp_kmer)
{
     size_t rsize=0;
     std::vector<int> offsets(rdisp_kmer);
     offsets.push_back(rcounts_kmer[size-1]+rdisp_kmer[size-1]);
     for (int t=0; t<rcounts_kmer.size(); t++) rsize += rcounts_kmer[t];
     assert(kmer_recv_buf.size() == rsize);
     assert(offsets.size() == (size+1));

     while(offsets.size()>2) {
            //assert(offsets.back() == kmer_recv_buf.size());
            assert(offsets.front() == 0);
            std::vector<int> new_offsets;
            int x = 0;
            while(x+2 < offsets.size()) {
                    // mergesort (offsets[x],offsets[x+1]) and (offsets[x+1],offsets[x+2])
                    std::inplace_merge(kmer_recv_buf.begin()+offsets[x]
                                 ,kmer_recv_buf.begin()+offsets[x+1]
                                 ,kmer_recv_buf.begin()+offsets[x+2] // this *might* be at the end
                                 ,[](const auto& i, const auto& j) {return i.seq < j.seq;} 
                                 );
                    // now they are sorted, we just put offsets[x] and offsets[x+2] into the new offsets.
                    // offsets[x+1] is not relevant any more
                    new_offsets.push_back(offsets[x]);
                    new_offsets.push_back(offsets[x+2]);
                    x += 2;
            }
            // if the number of offsets was odd, there might be a dangling offset
            // which we must remember to include in the new_offsets
            if(x+2==offsets.size()) {
                    new_offsets.push_back(offsets[x+1]);
            }
            // assert(new_offsets.front() == 0);
            //assert(new_offsets.back() == kmer_recv_buf.size());
            offsets.swap(new_offsets);

    }
    offsets.clear();
    offsets.shrink_to_fit();
}  


/*void SortAndAggregate(std::vector<kmer_t>& arr, std::vector<int>& count)
{

     std::vector<int> add;
     std::vector<size_t> indices(count.size());
     std::iota(indices.begin(), indices.end(), 0);

     std::sort(indices.begin(), indices.end(), Comp(arr));
     sort(arr.begin(), arr.end());

     assert(indices.size() == arr.size());
     std::vector<int> temp(count.size());
     for (int it = 0; it < (int)indices.size(); it++){
             temp[it] = count[ indices[it] ];
     }

     count = temp;
     temp.clear();
     indices.clear();

     int aggr=count[0];
     kmer_t prev=arr[0];
     for(int i = 1; i < (int)(arr.size()); i++)
     {
          if (arr[i] == prev)
              aggr += count[i];
           else {
                add.push_back(aggr);
                aggr = count[i];
                prev = arr[i];
           }
     }
     add.push_back(aggr);

     arr.erase( unique( arr.begin(), arr.end() ), arr.end() );
     assert(arr.size() == add.size());
     count = add;

}
*/

extern bool USE_NEW_SORT_AGGREGATE;

// Inplace sort and aggreate
void SortAndAggregate_InPlace(std::vector<KmerPairs>& arr)
{
  if (arr.empty()) { return; }

  std::vector<KmerPairs>::iterator uniq_key_itr=arr.begin(), itr=arr.begin();
  int cum_sum = uniq_key_itr->k_count;

  for(++itr; itr!=arr.end(); ++itr)
  {
    if (itr->seq != uniq_key_itr->seq) {
      // update previous count //
      uniq_key_itr->k_count = cum_sum;
      // reset //
      ++uniq_key_itr;
      uniq_key_itr->seq = itr->seq;
      cum_sum = 0;
    }
    cum_sum += itr->k_count;
  }

  // trailing case //
  uniq_key_itr->k_count = cum_sum;
  arr.resize((uniq_key_itr - arr.begin())+1);
}

void SortAndAggregate(std::vector<KmerPairs>& arr)
{
    if (USE_NEW_SORT_AGGREGATE) {
      SortAndAggregate_InPlace(arr);
      return;
    }

    std::vector<KmerPairs>::iterator low,up, it;
    //auto cmp = [](const KmerPairs& i, const KmerPairs& j) { return i.seq < j.seq; };
    //size_t new_size = std::set<KmerPairs, decltype(cmp)>(arr.begin(), arr.end(), cmp).size();
    std::vector<KmerPairs> new_arr;
    //new_arr.reserve(new_size);

/*    
#ifdef MEMOPT
    auto cmp = [](const KmerPairs& i, const KmerPairs& j) { return i.seq < j.seq; };
    size_t new_size = std::set<KmerPairs, decltype(cmp)>(arr.begin(), arr.end(), cmp).size();
    std::vector<KmerPairs> new_arr (new_size);
    size_t pos=0;
#else
    std::vector<KmerPairs> new_arr;
#endif
*/ 

    for( it = arr.begin(); it != arr.end(); )
          {
              kmer_t key = (*it).seq;
              low=std::lower_bound (arr.begin(), arr.end(), key,
                                   [] (const KmerPairs& lhs, kmer_t rhs) {
                             return (lhs.seq < rhs);
                             });

              up= std::upper_bound (arr.begin(), arr.end(), key,
                                   [] (kmer_t rhs, const KmerPairs& lhs) {
                             return (rhs < lhs.seq);
                             });

              int sum=0;
              for (auto itr=low; itr!= up; itr++) {
                        sum += (*itr).k_count;
              }
/*
#ifdef MEMOPT
              new_arr[pos].seq= key;
              new_arr[pos].k_count=sum;
              pos++;
*/
              new_arr.push_back(KmerPairs{key, sum});

              it = up;
          }

     arr=new_arr;
     new_arr.clear();
     new_arr.shrink_to_fit();

}



void transfer_kmers (std::vector<int>& scounts_kmer, 
                     std::vector<KmerPairs> &kmer_send_buf) 
{

    int ssize=0, rsize=0;// disp=0;
    std::vector<int> rcounts_kmer (size,0);
    std::vector<int> rdisp_kmer (size,0);
    std::vector<int> sdisp_kmer (size,0);

    for (int t=0; t<size; t++) ssize += scounts_kmer[t];

    sdisp_kmer[0] = 0;
    for (int i=1; i<size; i++) sdisp_kmer[i] = scounts_kmer[i-1] + sdisp_kmer[i-1];

     //create contiguous derived data type
     MPI_Datatype rowtype;
     MPI_Type_contiguous(sizeof(KmerPairs), MPI_BYTE, &rowtype);
     MPI_Type_commit(&rowtype);

    MPI_Barrier(MPI_COMM_WORLD);

    double t7 = MPI_Wtime ();
    MPI_Alltoall (scounts_kmer.data(), 1, MPI_INT, rcounts_kmer.data(), 1, MPI_INT, MPI_COMM_WORLD);
    //MPI_Barrier(MPI_COMM_WORLD);
    double t8 = MPI_Wtime ();
    alltoall_time += (t8 - t7);

    for (int t=0; t<size; t++) rsize += rcounts_kmer[t];
    rdisp_kmer[0] = 0;
    for (int i=1; i<size; i++) rdisp_kmer[i] = rcounts_kmer[i-1] + rdisp_kmer[i-1];

    std::vector<KmerPairs> kmer_recv_buf (rsize);

    MPI_Barrier(MPI_COMM_WORLD);

    double t9 = MPI_Wtime ();
     MPI_Alltoallv(kmer_send_buf.data(), scounts_kmer.data(), sdisp_kmer.data(), rowtype,
                   kmer_recv_buf.data(), rcounts_kmer.data(), rdisp_kmer.data(), rowtype, 
                   MPI_COMM_WORLD);
    //MPI_Barrier(MPI_COMM_WORLD);
    double t10 = MPI_Wtime ();
    alltoallv_time += (t10 - t9);

    kmer_send_buf.clear();
    kmer_send_buf.shrink_to_fit();

    double t21 = MPI_Wtime ();

    // sort the recv buffer
    double tm1 = MPI_Wtime ();
    //sort_recv_buffer(kmer_recv_buf, rcounts_kmer, rdisp_kmer);
    std::sort(kmer_recv_buf.begin(), kmer_recv_buf.end(),
              [](const auto& i, const auto& j) {return i.seq < j.seq;}
             );

    double tm2 = MPI_Wtime ();
    unpack_rbuf_sort += (tm2 - tm1);

    // Insert and sort
    double tm3 = MPI_Wtime ();
    if (kmer_proc_buf.size())
    {
       size_t offset = kmer_proc_buf.size();
       kmer_proc_buf.insert(kmer_proc_buf.end(), kmer_recv_buf.begin(), kmer_recv_buf.end());

       std::inplace_merge(kmer_proc_buf.begin(),
                  kmer_proc_buf.begin()+offset,
                  kmer_proc_buf.end(), // this *might* be at the end
                  [](const auto& i, const auto& j) { return i.seq < j.seq; }
                 );
    }
    else
      kmer_proc_buf.insert(kmer_proc_buf.end(), kmer_recv_buf.begin(), kmer_recv_buf.end());
       
    double tm4 = MPI_Wtime ();
    unpack_rbuf_insert += (tm4 - tm3);

    kmer_recv_buf.clear();
    kmer_recv_buf.shrink_to_fit();

    double tm5 = MPI_Wtime ();
    SortAndAggregate (kmer_proc_buf);
    double tm6 = MPI_Wtime ();
    unpack_rbuf_acc += (tm6 - tm5);


    double t22 = MPI_Wtime ();
    unpack_rbuf_time += (t22-t21);
    num_batch_transfers++;
  
    // free datatype
    MPI_Type_free(&rowtype);
    
    rcounts_kmer.clear();
    rdisp_kmer.clear();
    sdisp_kmer.clear();
     
 
}


void recalculate_min_lmer (kmer_t kmer_in, lmer_t *m_lmer, lmer_t *m_lmer_freq, int *m_pos)
{
    lmer_t min_lmer=0, tmp_lmer=0;
    lmer_t min_lmer_freq=0, tmp_lmer_freq=0;
    int min_pos=0, k=0;

    for (k=0; ((KMER_LENGTH-1) - k) >= (LMER_LENGTH-1); k++) {
        lmer_t lmer_out=0;
        for(int j=k; j<LMER_LENGTH+k; j++) {
            lmer_out = kmer_to_lmer (kmer_in, j, lmer_out);
        }

        tmp_lmer = lmer_out;
        tmp_lmer_freq = global_lmer_frequency[tmp_lmer];

        if (k == 0) {
            min_lmer = tmp_lmer;
            min_lmer_freq = tmp_lmer_freq;
            min_pos = 0;
        }
        else {
           if (tmp_lmer_freq < min_lmer_freq) {
               min_lmer = tmp_lmer;
               min_lmer_freq = tmp_lmer_freq;
               min_pos = k;
           }
        }
    }
    assert (k == (KMER_LENGTH-LMER_LENGTH+1));

    *m_lmer = min_lmer;
    *m_lmer_freq = min_lmer_freq;
    *m_pos = min_pos;
}


void Sliding_window_l (const char *ptr, size_t length) {

#ifdef LMER_DEBUG2
      ElType this_alpha;
      char kmer_out[lmer_len+1];
#endif 

  size_t p=0;
  //lmer_frequency.reserve(pow(4, LMER_LENGTH));

  /*find start of a read*/
  for(; ptr[p]!='>' && p<length; p++) {/*noop*/ }

  lmer_t kmer = 0;

  while(p<length) {
    assert(ptr[p]=='>'); /*this will be true*/

    /*skip till newline*/
    for(; p<length && ptr[p]!='\n'; p++) {/*noop*/ }
    p++; /*skip the newline*/

    if(p+LMER_LENGTH > length) break; /*too short a read*/
    kmer = 0;
    int i;
    for(i=0; ptr[p]!='\n' && i<LMER_LENGTH-1; i++) {
      kmer = lmer_shift(kmer, char_to_el(ptr[p++]));
      //kmer = kmer_cons(kmer, i, char_to_el(ptr[p++]));
    }

    while(p<length && ptr[p]!='\n') {
      kmer = lmer_shift(kmer, char_to_el(ptr[p++]));

#ifdef LMER_DEBUG2
      printf("lmer: %lu,");
      for(int j=0; j<LMER_LENGTH; j++) {
          this_alpha = lmerel (kmer, j);
          lmer_out[j] = el_to_char(this_alpha);
      }
      lmer_out[lmer_len] = '\0';
      printf(" lmer_out: %s \n", lmer_out);
#endif

      lmer_frequency[kmer]++;
    }
    p++; /*skip the newline*/
  }
}

////////////////////////////////////////////////////////////////////////////////
static inline const char* skip_till_symbol(const char *beg, const char *end, 
    char symbol) {
  for (; (beg != end) && (*beg != symbol); ++beg) { }
  return beg;
}

////////////////////////////////////////////////////////////////////////////////

struct kmer_strncmp_t {
  kmer_strncmp_t(size_t klen=32UL) : klen_(klen) {}

  inline bool operator() (const char* a, const char*b) const {
    return strncmp(a, b, klen_) < 0; 
  }

  inline bool operator() (const std::pair<const char *, size_t>& a,
      const std::pair<const char *, size_t>& b) const {
    return strncmp(a.first, b.first, klen_) < 0;
  }
  size_t klen_;
}; // struct kmer_strncmp_t //


struct KmerPairsImplicitCompare {

  KmerPairsImplicitCompare(size_t klen=32UL) : klen_(klen) {}

  inline bool operator() (const KmerPairsImplicit& a,
        const KmerPairsImplicit& b) const {
    return strncmp(a.kmer_ptr_, b.kmer_ptr_, klen_) < 0; 
  }

  size_t klen_;
}; // struct kmer_strncmp_t //

struct KmerPairsImplicitComparePartID {

  KmerPairsImplicitComparePartID(size_t klen=32UL) : klen_(klen) {}

  inline bool operator() (const KmerPairsImplicit& a,
        const KmerPairsImplicit& b) const {
    return (a.k_count == b.k_count) ? 
      (strncmp(a.kmer_ptr_, b.kmer_ptr_, klen_) < 0) : (a.k_count < b.k_count); 
  }

  size_t klen_;
}; // struct kmer_strncmp_t //



struct KmerPairsImplicitUtils {
  inline static int get_part_id(const char *kmer_ptr, int klen) {
    kmer_t kmer = 0;
    for (int i=0; i<KMER_LENGTH; i++) {
      kmer = kmer_shift(kmer, char_to_el(kmer_ptr[i])); 
    }
    return retrieve_proc_id(kmer);
  }
};



struct KmerPairCompare {
  inline bool operator() (const KmerPairs& a ,const KmerPairs& b) const {
    return a.seq < b.seq;
  }
}; // struct KmerPairCompare //

struct KmerPairCompareHashed {
  inline bool operator() (const KmerPairs& a ,const KmerPairs& b) const {
    return (a.k_count == b.k_count) ? (a.seq < b.seq) : (a.k_count < b.k_count);
  }
}; // struct KmerPairCompareHashed //
////////////////////////////////////////////////////////////////////////////////



//TODO(vamsikku): algorithm
//
// 0. use a atomic STATE_VARIABLE which
// 1. use a kmer_pair_buffer[ ][odd_even] of MAX_KMER_COUNT
// 2. use atomic variable to fill this array in parallel


////////////////////////////////////////////////////////////////////////////////

void Sliding_window_concurrent_v1_partition_and_sort_extra_space(const char *ptr,
      size_t length, size_t *n_kmers, 
      std::vector<std::vector<kmer_t>> &partial_kmer_counts) {
  size_t p=0;
  std::vector<int> scounts_kmer (size,0);
  std::vector<KmerPairs> kmer_pair_buffer; //k_count as hash place holder//
  std::vector<KmerPairs> kmer_pair_buffer_partition;
  size_t kpos=0;
  size_t num_kmers;
  {
      const char *symbol_begin_ = ptr;
      const char *symbol_end_ = ptr + length;
      const char skip_symbol_ = '>';
      const char delimiter_ = '\n';
      size_t read_len_ = 100UL;
      size_t klen_ = (size_t) KMER_LENGTH;
      // Find the k-mers a thread would process //
      size_t tcount = 0;
#pragma omp parallel
      {
        if (!omp_get_thread_num())
          tcount = omp_get_num_threads();
      }

      std::vector<size_t> threads_part_count(tcount, 0UL);
      std::vector<std::pair<const char *, const char *> > 
          threads_begin_end(tcount); 
      std::atomic<int> partition_counts[size] = { {0} };
      int kmer_part_id;

      // STEP-1 find part counts of each thread in parallel //
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
        }
      }

      {
        size_t cum = 0, aux;
        for (size_t i=0; i<tcount; i++) {
          aux = threads_part_count[i];
          threads_part_count[i] = cum;
          cum += aux;
        }
        kmer_pair_buffer.resize(cum);
        kmer_pair_buffer_partition.resize(cum);
        num_kmers = cum;
      }
      

      // STEP-2: extract k-mers and compute partition counts //
#pragma omp parallel
      {
        size_t my_id = omp_get_thread_num();
        const auto &my_begin_end = threads_begin_end[my_id];

        const char *beg = my_begin_end.first;
        const char *end = my_begin_end.second;
        size_t start_idx = threads_part_count[my_id];   
        size_t next_idx;
        int kmer_part_id;

        while (beg < end) {
          for (; (beg<end) && (*beg != '>'); ++beg) { }
          for(; (beg<end) && (*beg !='\n'); beg++) {/*noop*/}
          if ((beg+KMER_LENGTH) >= end) { break; }

          ++beg;

          kmer_t kmer=0;
          for (int i=0; (*beg != '\n') && (i<KMER_LENGTH-1); ++i, ++beg) {
            kmer = kmer_shift(kmer, char_to_el(*beg));
          }

          while ((beg < end) && (*beg != '\n')) {
            if (start_idx >= kmer_pair_buffer.size()) {
              throw std::logic_error("[start_idx]: invariant failed start_idx="
                  +std::to_string(start_idx)+" kmer_pair_buffer.size()="
                  +std::to_string(kmer_pair_buffer.size()));
            }
            kmer = kmer_shift(kmer, char_to_el(*beg));
            kmer_part_id = retrieve_proc_id(kmer);
            kmer_pair_buffer[start_idx++] = KmerPairs{kmer, kmer_part_id};
            partition_counts[ kmer_part_id ].fetch_add(1,
                  std::memory_order_relaxed);
            ++beg;
          }
          ++beg;
        }
      }

      // STEP-3: partition the k-mers into kmer_pair_buffer_partition //

      std::vector<size_t> partition_begin(size, 0);
      {
        std::atomic<size_t> next_avail_part_idx[size] = { {0UL} };
        {
          size_t part_cum = 0UL;
          for (size_t i=0; i<size; i++) {
            next_avail_part_idx[i].store(part_cum);
            partition_begin[i] = part_cum;
            part_cum += partition_counts[i];
          }
        }

#pragma omp parallel for 
        for (size_t i=0; i<kmer_pair_buffer.size(); i++) {
          int kmer_part_id = kmer_pair_buffer[i].k_count;
          size_t next_idx = next_avail_part_idx[ kmer_part_id ].fetch_add(1,
              std::memory_order_relaxed);
          kmer_pair_buffer_partition[ next_idx ] = kmer_pair_buffer[i];
        }
      }

      // STEP-4: sort each of the partitions in parallel //
      {

        // parallel sort each partitions rather than entire k-mer array //
#pragma omp parallel for
        for (size_t i=0; i<size; i++) {
          KmerPairCompare kmer_pair_compare;
          auto itr_beg = kmer_pair_buffer_partition.begin() + partition_begin[i];
          auto itr_end = itr_beg + partition_counts[i];
          oneapi::tbb::parallel_sort(itr_beg, itr_end, kmer_pair_compare);
          partition_counts[i].store(0);
        }
      }
      kmer_pair_buffer.clear();

      // STEP-5: parallel duplicate removal and compute skmer_counts //
      {

        std::vector<size_t> uniq_kmer_counts(tcount, 0UL);
        std::vector<size_t> kmer_starts(tcount, 0UL);


        // parallel duplicate removal, frequency count and partition count 
#pragma omp parallel
        {
          size_t my_id = omp_get_thread_num(); 
          size_t data_size = kmer_pair_buffer_partition.size(); 
          size_t chunk_size = data_size/tcount;
          int curr_freq = 1, kmer_part_id;

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
            kmer_t prev_kmer = kmer_pair_buffer_partition[start_idx-1UL].seq;
            while ( (start_idx < end_idx) && 
                    (prev_kmer == kmer_pair_buffer_partition[start_idx].seq )) {
              ++start_idx;
            }
          }
          kmer_starts[my_id] = start_idx;
          size_t uniq_idx = start_idx;

          for (size_t i=start_idx+1; i<end_idx; i++) {
            if (kmer_pair_buffer_partition[uniq_idx].seq == 
                  kmer_pair_buffer_partition[i].seq) {
              curr_freq++;
            } else {
              kmer_part_id = 
                retrieve_proc_id(kmer_pair_buffer_partition[uniq_idx].seq);
              partition_counts[kmer_part_id].fetch_add(1UL,
                      std::memory_order_relaxed);
              // update the frequency //
              kmer_pair_buffer_partition[uniq_idx].k_count = curr_freq;
              // copy the new unique kmer into next slot //
              kmer_pair_buffer_partition[++uniq_idx] = 
                  kmer_pair_buffer_partition[i];
              curr_freq = 1;
            }
          }

          // the k-mer at the boundary may have partial frequency //
          for (size_t i=end_idx; (i<data_size) && 
                (kmer_pair_buffer_partition[uniq_idx].seq == 
                  kmer_pair_buffer_partition[i].seq); i++) 
          {
            curr_freq++; 
          }

          // trailing case //
          {
            kmer_part_id = 
                retrieve_proc_id(kmer_pair_buffer_partition[uniq_idx].seq);
            partition_counts[kmer_part_id].fetch_add(1UL,
                  std::memory_order_relaxed);
            kmer_pair_buffer_partition[uniq_idx].k_count = curr_freq;
          }
          uniq_kmer_counts[ my_id ] = (uniq_idx - start_idx) + 1UL;
        }

#pragma omp parallel for
        for (size_t i=0; i<size; i++) {
          scounts_kmer[i] = int(partition_counts[i]);
        }

        // prefix sum //
        // TODO(vamsikku): 
        std::vector<size_t> prefix_sum(tcount, 0UL);
        size_t uniq_cum = uniq_kmer_counts[0UL]; 
        prefix_sum[0UL] = 0UL;
        for (size_t i=1; i<tcount; i++) {
          prefix_sum[i] = uniq_cum;
          uniq_cum += uniq_kmer_counts[i];
        }
        kmer_pair_buffer.resize(uniq_cum);


        // do a parallel write into kmer_send_buffer //
#pragma omp parallel 
        {
          size_t my_id = omp_get_thread_num(); 
          size_t start_idx = kmer_starts[my_id];
          size_t end_idx = start_idx + uniq_kmer_counts[my_id];
          size_t kmer_send_buffer_idx = prefix_sum[my_id];

          for (size_t i=start_idx; i<end_idx; i++) {
            if (kmer_send_buffer_idx >= kmer_pair_buffer_partition.size()) {
              throw std::logic_error("kmer_send_buffer invariant failed " 
                  " tcount="+std::to_string(tcount) +
                  +" kmer_send_buffer_idx="+
                    std::to_string(kmer_send_buffer_idx)
                  +" kmer_send_buffer.size()="+
                    std::to_string(kmer_pair_buffer_partition.size())
                 );
            }
            kmer_pair_buffer[kmer_send_buffer_idx++] = 
                kmer_pair_buffer_partition[i];
          }
        }
        ////////////////////////////////////////////////////////////////////

        size_t total_scounts = 0;

        for (size_t i=0; i<size; i++) {
          total_scounts += scounts_kmer[i];
        }

        if (total_scounts != kmer_pair_buffer.size()) {
          printf("total_soucnts = %lu kmer_send_buffer.size()=%lu\n",
              total_scounts, kmer_pair_buffer.size());
          throw std::logic_error("[total_scounts != kmer_send_buffer.size()]\n");
        }

        {
          // TRANSFER and reset //
          transfer_kmers (scounts_kmer, kmer_pair_buffer);
          num_kmers = 0;
          kmer_pair_buffer.clear();
          kmer_pair_buffer_partition.clear();
        }
      }
  }
  *n_kmers = 0;
}


void Sliding_window_concurrent_v1_partition_and_sort(const char *ptr,
      size_t length, size_t *n_kmers, 
      std::vector<std::vector<kmer_t>> &partial_kmer_counts) {
  size_t p=0;
  std::vector<int> scounts_kmer (size,0);
  std::vector<KmerPairs> kmer_pair_buffer; //k_count as hash place holder//
  std::vector<KmerPairs> kmer_send_buffer;
  size_t kpos=0;
  size_t num_kmers;
  {
      const char *symbol_begin_ = ptr;
      const char *symbol_end_ = ptr + length;
      const char skip_symbol_ = '>';
      const char delimiter_ = '\n';
      size_t read_len_ = 100UL;
      size_t klen_ = (size_t) KMER_LENGTH;
      // Find the k-mers a thread would process //
      size_t tcount = 0;
#pragma omp parallel
      {
        if (!omp_get_thread_num())
          tcount = omp_get_num_threads();
      }

      std::vector<size_t> threads_part_count(tcount, 0UL);
      std::vector<std::pair<const char *, const char *> > 
          threads_begin_end(tcount); 
      std::atomic<int> partition_counts[size] = { {0} };
      int kmer_part_id;

      // STEP-1 find part counts of each thread in parallel //
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

          { // determine the partition counts //
            const char *beg = read_beg;
            const char *end = read_end;

            do {
              kmer_t kmer=0;
              for (int i=0; (*beg != '\n') && (i<KMER_LENGTH-1); ++i, ++beg) {
                kmer = kmer_shift(kmer, char_to_el(*beg));
              }

              while ((beg < end) && (*beg != '\n')) {
                kmer = kmer_shift(kmer, char_to_el(*beg));
                kmer_part_id = retrieve_proc_id(kmer);
                partition_counts[kmer_part_id].fetch_add(1UL,
                      std::memory_order_relaxed);
                ++beg;
              }
              ++beg;
              symbol_begin = beg;
            } while (0);
          }

          my_part_count += (read_len-klen_+1UL);
        }
      }

      size_t cum = threads_part_count[0];
      for (size_t i=1; i<tcount; i++) {
        cum += threads_part_count[i];
        threads_part_count[i] = cum;
      }
      

      // processor i (i > 0) 
      // start = threads_part_count[i-1], end = threads_part_count[i] //
      // threads_part_count[-1] = 0 //

      kmer_pair_buffer.resize(cum);
      num_kmers = cum;

      std::atomic<size_t> next_avail_part_idx[size] = { {0UL} };
      std::vector<size_t> partition_begin(size, 0);
      {
        size_t part_cum = 0UL;
        for (size_t i=0; i<size; i++) {
          next_avail_part_idx[i].store(part_cum);
          partition_begin[i] = part_cum;
          part_cum += partition_counts[i];
        }
      }


      // every thread extracts k-mers in parallel and fill kmer_partitions//
#pragma omp parallel
      {
        size_t my_id = omp_get_thread_num();
        const auto &my_begin_end = threads_begin_end[my_id];

        const char *beg = my_begin_end.first;
        const char *end = my_begin_end.second;
        size_t next_idx;
        int kmer_part_id;

        while (beg < end) {
          for (; (beg<end) && (*beg != '>'); ++beg) { }
          for(; (beg<end) && (*beg !='\n'); beg++) {/*noop*/}
          if ((beg+KMER_LENGTH) >= end) { break; }

          ++beg;

          kmer_t kmer=0;
          for (int i=0; (*beg != '\n') && (i<KMER_LENGTH-1); ++i, ++beg) {
            kmer = kmer_shift(kmer, char_to_el(*beg));
          }

          while ((beg < end) && (*beg != '\n')) {
            kmer = kmer_shift(kmer, char_to_el(*beg));
            kmer_part_id = retrieve_proc_id(kmer);
            next_idx = next_avail_part_idx[kmer_part_id].fetch_add(1,
                  std::memory_order_relaxed);
            if (next_idx >= kmer_pair_buffer.size()) {
              throw std::logic_error("[next_avail_part_idx]: invariant failed");
            }
            kmer_pair_buffer[ next_idx ] = KmerPairs{kmer, 1};
            ++beg;
          }
          ++beg;
        }
      }

      {
        KmerPairCompare kmer_pair_compare;

        // parallel sort each partitions rather than entire k-mer array //
#pragma omp parallel for
        for (size_t i=0; i<size; i++) {
          auto itr_beg = kmer_pair_buffer.begin() + partition_begin[i];
          auto itr_end = itr_beg + partition_counts[i];
          oneapi::tbb::parallel_sort( itr_beg, itr_end, kmer_pair_compare);
          partition_counts[i].store(0);
        }

        std::vector<size_t> uniq_kmer_counts(tcount, 0UL);
        std::vector<size_t> kmer_starts(tcount, 0UL);

        // reset the partition counts //


        // parallel duplicate removal, frequency count and partition count 
#pragma omp parallel
        {
          size_t my_id = omp_get_thread_num(); 
          size_t data_size = kmer_pair_buffer.size(); 
          size_t chunk_size = data_size/tcount;
          int curr_freq = 1, kmer_part_id;

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
          size_t uniq_idx = start_idx;

          for (size_t i=start_idx+1; i<end_idx; i++) {
            if (kmer_pair_buffer[uniq_idx].seq == kmer_pair_buffer[i].seq) {
              curr_freq++;
            } else {
              kmer_part_id = retrieve_proc_id(kmer_pair_buffer[uniq_idx].seq);
              partition_counts[kmer_part_id].fetch_add(1UL,
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
            kmer_part_id = retrieve_proc_id(kmer_pair_buffer[uniq_idx].seq);
            partition_counts[kmer_part_id].fetch_add(1UL,
                  std::memory_order_relaxed);
            kmer_pair_buffer[uniq_idx].k_count = curr_freq;
          }
          uniq_kmer_counts[ my_id ] = (uniq_idx - start_idx) + 1UL;
        }

#pragma omp parallel for
        for (size_t i=0; i<size; i++) {
          scounts_kmer[i] = int(partition_counts[i]);
        }

        // prefix sum //
        // TODO(vamsikku): 
        std::vector<size_t> prefix_sum(tcount, 0UL);
        size_t uniq_cum = uniq_kmer_counts[0UL]; 
        prefix_sum[0UL] = 0UL;
        for (size_t i=1; i<tcount; i++) {
          prefix_sum[i] = uniq_cum;
          uniq_cum += uniq_kmer_counts[i];
        }
        kmer_send_buffer.resize(uniq_cum);

        

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

        size_t total_scounts = 0;

        for (size_t i=0; i<size; i++) {
          total_scounts += scounts_kmer[i];
        }

        if (total_scounts != kmer_send_buffer.size()) {
          printf("total_soucnts = %lu kmer_send_buffer.size()=%lu\n",
              total_scounts, kmer_send_buffer.size());
          throw std::logic_error("[total_scounts != kmer_send_buffer.size()]\n");
        }

        {
          // TRANSFER and reset //
          transfer_kmers (scounts_kmer, kmer_send_buffer);
          num_kmers = 0;
          kmer_send_buffer.clear();
        }
      }
  }
  *n_kmers = 0;
}





///////////////////////////////////////////////////////////////////////////////

void Sliding_window_concurrent_v1_implicit(const char *ptr,
      size_t length, size_t *n_kmers, 
      std::vector< std::vector<kmer_t> > &partial_kmer_counts) {

  size_t p=0;
  std::vector<int> scounts_kmer (size,0);
  std::vector<KmerPairsImplicit> kmer_pair_buffer;
  std::vector<KmerPairs> kmer_send_buffer;
  size_t kpos=0;
  size_t num_kmers;

  {
      const char *symbol_begin_ = ptr;
      const char *symbol_end_ = ptr + length;
      const char skip_symbol_ = '>';
      const char delimiter_ = '\n';
      size_t read_len_ = 100UL;
      size_t klen_ = (size_t) KMER_LENGTH;
      // Find the k-mers a thread would process //
      size_t tcount = 0;
#pragma omp parallel
      {
        if (!omp_get_thread_num())
          tcount = omp_get_num_threads();
      }

      std::vector<size_t> threads_part_count(tcount, 0UL);
      std::vector<std::pair<const char *, const char *> > 
          threads_begin_end(tcount); 

      // STEP-1 find part counts of each thread in parallel //
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
        }
      }

      size_t cum = threads_part_count[0];
      for (size_t i=1; i<tcount; i++) {
        cum += threads_part_count[i];
        threads_part_count[i] = cum;
      }

      // processor i (i > 0) 
      // start = threads_part_count[i-1], end = threads_part_count[i] //
      // threads_part_count[-1] = 0 //

      kmer_pair_buffer.resize(cum);
      num_kmers = cum;

      // every thread extracts k-mers in parallel and fill kmer_partitions//
#pragma omp parallel
      {
        size_t my_id = omp_get_thread_num();
        size_t start_idx = my_id ? threads_part_count[my_id-1] : 0UL;;
        size_t end_idx = threads_part_count[my_id];
        const auto &my_begin_end = threads_begin_end[my_id];
        

        const char *beg = my_begin_end.first;
        const char *end = my_begin_end.second;
        int read_len = 0UL;

        while (beg < end) {
          for (; (beg<end) && (*beg != '>'); ++beg) { }
          for(; (beg<end) && (*beg !='\n'); beg++) {/*noop*/}
          if ((beg+KMER_LENGTH) >= end) { break; }

          ++beg;

          const char *read_begin = beg;
          // compute the readlen /
          read_len = 0UL;
          while ((beg < end) && (*beg != '\n')) { ++beg; ++read_len;}
          if (read_len > KMER_LENGTH) {

            for (int i=0; i<((read_len-KMER_LENGTH)+1); i++, start_idx++) {
              //NOTE: we need to compute partition ids after sorting //
              kmer_pair_buffer[start_idx] = 
                  KmerPairsImplicit{(read_begin+i), 
                  KmerPairsImplicitUtils::get_part_id( 
                      (read_begin+i), KMER_LENGTH) };
            }

          }
          ++beg;
        }

        if (start_idx != end_idx) {
          throw std::logic_error("Invalid partition size");
        }
      }

      {
        KmerPairsImplicitComparePartID kmer_pair_compare(KMER_LENGTH);
        {
          oneapi::tbb::parallel_sort(kmer_pair_buffer.begin(), 
              kmer_pair_buffer.begin()+num_kmers, kmer_pair_compare);
        }

        //TODO(vamsikku): parallel duplicate removal must be adapted to
        
        
        std::atomic<int> partition_counts[size] = { {0} };
        std::vector<size_t> uniq_kmer_counts(tcount, 0UL);
        std::vector<size_t> kmer_starts(tcount, 0UL);

        // parallel duplicate removal, frequency count and partition count 
#pragma omp parallel
        {
          size_t my_id = omp_get_thread_num(); 
          size_t data_size = num_kmers; 
          size_t chunk_size = data_size/tcount;
          int kmer_part_id;

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
            const char *prev_kmer = kmer_pair_buffer[start_idx-1UL].kmer_ptr_;
            while ( (start_idx < end_idx) && 
                    ( !strncmp(prev_kmer, kmer_pair_buffer[start_idx].kmer_ptr_,
                                KMER_LENGTH)) ) {
              ++start_idx;
            }
          }
          kmer_starts[my_id] = start_idx;
          int curr_freq = 1;
          size_t uniq_idx = start_idx;

          for (size_t i=start_idx+1; i<end_idx; i++) {
            if (!strncmp(kmer_pair_buffer[uniq_idx].kmer_ptr_, 
                  kmer_pair_buffer[i].kmer_ptr_, KMER_LENGTH)) {
              curr_freq++;
            } else {
              // update the partition counts //
              kmer_part_id = kmer_pair_buffer[uniq_idx].k_count;
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
                !(strncmp(kmer_pair_buffer[uniq_idx].kmer_ptr_, 
                    kmer_pair_buffer[i].kmer_ptr_, KMER_LENGTH)); i++) 
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

#pragma omp parallel for
        for (size_t i=0; i<size; i++) {
          scounts_kmer[i] = int(partition_counts[i]);
        }


        // prefix sum //
        std::vector<size_t> prefix_sum(tcount, 0UL);
        size_t cum = uniq_kmer_counts[0UL]; 
        prefix_sum[0UL] = 0UL;
        for (size_t i=1; i<tcount; i++) {
          prefix_sum[i] = cum;
          cum += uniq_kmer_counts[i];
        }
        kmer_send_buffer.resize(cum);

        // now explicitly extract the kmers //
#pragma omp parallel 
        {
          size_t my_id = omp_get_thread_num(); 
          size_t start_idx = kmer_starts[my_id];
          size_t end_idx = start_idx + uniq_kmer_counts[my_id];
          size_t kmer_send_buffer_idx = prefix_sum[my_id];
          int kmer_part_id, kmer_freq;
          const char *kmer_ptr;

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
            kmer_ptr = kmer_pair_buffer[i].kmer_ptr_;
            kmer_freq = kmer_pair_buffer[i].k_count;

            // extract the k-mer //
            kmer_t kmer = 0;
            for (int i=0; i<KMER_LENGTH; i++) {
              kmer = kmer_shift(kmer, char_to_el(kmer_ptr[i])); 
            }
            kmer_send_buffer[kmer_send_buffer_idx].seq = kmer;
            kmer_send_buffer[kmer_send_buffer_idx++].k_count = kmer_freq;
            // update the partition counts //
            kmer_part_id = retrieve_proc_id(kmer);
            partition_counts[kmer_part_id].fetch_add(1,
                std::memory_order_relaxed);
          }
        }
        kmer_pair_buffer.clear();

        ////////////////////////////////////////////////////////////////////
        {
          // TRANSFER and reset //
          transfer_kmers (scounts_kmer, kmer_send_buffer);
          num_kmers = 0;
          kmer_send_buffer.clear();
        }
      }

  }

  *n_kmers = 0;
}

// Implicit mode: we don't explicitly encode the k-mers//
void Sliding_window_concurrent_v2_old(const char *ptr,
      size_t length, size_t *n_kmers, 
      std::vector<std::vector<kmer_t>> &partial_kmer_counts) {

  size_t p=0;
  std::vector<int> scounts_kmer (size,0);
  std::vector<KmerPairsImplicit> kmer_pair_buffer;
  std::vector<KmerPairs> kmer_send_buffer, kmer_send_buffer_final;
  size_t kpos=0;
  size_t num_kmers;

  {
      const char *symbol_begin_ = ptr;
      const char *symbol_end_ = ptr + length;
      const char skip_symbol_ = '>';
      const char delimiter_ = '\n';
      size_t read_len_ = 100UL;
      size_t klen_ = (size_t) KMER_LENGTH;
      // Find the k-mers a thread would process //
      size_t tcount = 0;
#pragma omp parallel
      {
        if (!omp_get_thread_num())
          tcount = omp_get_num_threads();
      }

      std::vector<size_t> threads_part_count(tcount, 0UL);
      std::vector<std::pair<const char *, const char *> > 
          threads_begin_end(tcount); 

      // STEP-1 find part counts of each thread in parallel //
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
        }
      }

      size_t cum = threads_part_count[0];
      for (size_t i=1; i<tcount; i++) {
        cum += threads_part_count[i];
        threads_part_count[i] = cum;
      }

      // processor i (i > 0) 
      // start = threads_part_count[i-1], end = threads_part_count[i] //
      // threads_part_count[-1] = 0 //

      kmer_pair_buffer.resize(cum);
      num_kmers = cum;

      // every thread extracts k-mers in parallel and fill kmer_partitions//
#pragma omp parallel
      {
        size_t my_id = omp_get_thread_num();
        size_t start_idx = my_id ? threads_part_count[my_id-1] : 0UL;;
        size_t end_idx = threads_part_count[my_id];
        const auto &my_begin_end = threads_begin_end[my_id];
        

        const char *beg = my_begin_end.first;
        const char *end = my_begin_end.second;
        int read_len = 0UL;

        while (beg < end) {
          for (; (beg<end) && (*beg != '>'); ++beg) { }
          for(; (beg<end) && (*beg !='\n'); beg++) {/*noop*/}
          if ((beg+KMER_LENGTH) >= end) { break; }

          ++beg;

          const char *read_begin = beg;
          // compute the readlen /
          read_len = 0UL;
          while ((beg < end) && (*beg != '\n')) { ++beg; ++read_len;}
          if (read_len > KMER_LENGTH) {

            for (int i=0; i<((read_len-KMER_LENGTH)+1); i++, start_idx++) {
              //NOTE: we need to compute partition ids after sorting //
              kmer_pair_buffer[start_idx] = 
                  KmerPairsImplicit{(read_begin+i), 1UL};
            }

          }
          ++beg;
        }

        if (start_idx != end_idx) {
          throw std::logic_error("Invalid partition size");
        }
      }

      {
        KmerPairsImplicitCompare kmer_pair_compare(KMER_LENGTH);
        {
          oneapi::tbb::parallel_sort(kmer_pair_buffer.begin(), 
              kmer_pair_buffer.begin()+num_kmers, kmer_pair_compare);
        }

        //TODO(vamsikku): parallel duplicate removal must be adapted to
        
        
        std::atomic<int> partition_counts[size] = { {0} };
        std::vector<size_t> uniq_kmer_counts(tcount, 0UL);
        std::vector<size_t> kmer_starts(tcount, 0UL);

        // parallel duplicate removal, frequency count and partition count 
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
            const char *prev_kmer = kmer_pair_buffer[start_idx-1UL].kmer_ptr_;
            while ( (start_idx < end_idx) && 
                    ( !strncmp(prev_kmer, kmer_pair_buffer[start_idx].kmer_ptr_,
                                KMER_LENGTH)) ) {
              ++start_idx;
            }
          }
          kmer_starts[my_id] = start_idx;
          int curr_freq = 1, kmer_part_id;
          size_t uniq_idx = start_idx;

          for (size_t i=start_idx+1; i<end_idx; i++) {
            if (!strncmp(kmer_pair_buffer[uniq_idx].kmer_ptr_, 
                  kmer_pair_buffer[i].kmer_ptr_, KMER_LENGTH)) {
              curr_freq++;
            } else {
              // update the frequency //
              kmer_pair_buffer[uniq_idx].k_count = curr_freq;
              // copy the new unique kmer into next slot //
              kmer_pair_buffer[++uniq_idx] = kmer_pair_buffer[i];
              curr_freq = 1;
            }
          }

          // the k-mer at the boundary may have partial frequency //
          for (size_t i=end_idx; (i<data_size) && 
                !(strncmp(kmer_pair_buffer[uniq_idx].kmer_ptr_, 
                    kmer_pair_buffer[i].kmer_ptr_, KMER_LENGTH)); i++) 
          {
            curr_freq++; 
          }

          // trailing case //
          {
            // update the frequency //
            kmer_pair_buffer[uniq_idx].k_count = curr_freq;
          }
          uniq_kmer_counts[ my_id ] = (uniq_idx - start_idx) + 1UL;
        }


        // prefix sum //
        // TODO(vamsikku): 
        std::vector<size_t> prefix_sum(tcount, 0UL);
        size_t cum = uniq_kmer_counts[0UL]; 
        prefix_sum[0UL] = 0UL;
        for (size_t i=1; i<tcount; i++) {
          prefix_sum[i] = cum;
          cum += uniq_kmer_counts[i];
        }
        kmer_send_buffer.resize(cum);

        // now explicitly extract the kmers and fill the partition counts //
#pragma omp parallel 
        {
          size_t my_id = omp_get_thread_num(); 
          size_t start_idx = kmer_starts[my_id];
          size_t end_idx = start_idx + uniq_kmer_counts[my_id];
          size_t kmer_send_buffer_idx = prefix_sum[my_id];
          int kmer_part_id, kmer_freq;
          const char *kmer_ptr;

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
            kmer_ptr = kmer_pair_buffer[i].kmer_ptr_;
            kmer_freq = kmer_pair_buffer[i].k_count;

            // extract the k-mer //
            kmer_t kmer = 0;
            for (int i=0; i<KMER_LENGTH; i++) {
              kmer = kmer_shift(kmer, char_to_el(kmer_ptr[i])); 
            }
            kmer_send_buffer[kmer_send_buffer_idx].seq = kmer;
            kmer_send_buffer[kmer_send_buffer_idx++].k_count = kmer_freq;
            // update the partition counts //
            kmer_part_id = retrieve_proc_id(kmer);
            partition_counts[kmer_part_id].fetch_add(1,
                std::memory_order_relaxed);
          }
        }

        kmer_pair_buffer.clear();


#pragma omp parallel for
        for (size_t i=0; i<size; i++) {
          scounts_kmer[i] = int(partition_counts[i]);
        }

        kmer_send_buffer_final.resize(kmer_send_buffer.size());

        std::atomic<size_t> next_avail_part_idx[size] = { {0UL} };
        cum = 0UL;

        for (size_t i=0; i<size; i++) {
          next_avail_part_idx[i].store(cum);
          cum += scounts_kmer[i];
        }

        if (cum != kmer_send_buffer.size()) {
          throw std::logic_error("[Partition Invariant Failed]: uniq-kmers="
              +std::to_string(kmer_send_buffer.size())+" total-parts="
              +std::to_string(cum));
        }


        // in parallel read from kmer_send_buffer and write it to 
        // kmer_send_buffer_final at appropriate index //

#pragma omp parallel 
        {
          size_t my_id = omp_get_thread_num(); 
          size_t start_idx = prefix_sum[my_id]; //idx into the kmer_send_buffer 
          size_t end_idx = start_idx + uniq_kmer_counts[my_id];
          int kmer_part_id;
          size_t next_idx;
          kmer_t kmer;

          for (size_t i=start_idx; i<end_idx; i++) {
            if (i >= kmer_send_buffer.size()) {
              std::string message = "[kmer_send_buffer]: invariant violation "
                "i=" + std::to_string(i)+" n="+
                  std::to_string(kmer_send_buffer.size());

              throw std::logic_error(message.c_str());
            }

            // find the partition of kmer 
            // get the next index 
            // write //
            kmer = kmer_send_buffer[i].seq;
            kmer_part_id = retrieve_proc_id(kmer);
            next_idx = next_avail_part_idx[kmer_part_id].fetch_add(1,
                  std::memory_order_relaxed);
            if (next_idx >= kmer_send_buffer_final.size()) {
              std::string message = "next_idx invariant violated s="+
                std::to_string(kmer_send_buffer_final.size());
              message += " idx="+std::to_string(next_idx);
              throw std::logic_error(message.c_str());
            }
            kmer_send_buffer_final[next_idx] = kmer_send_buffer[i];
          }
        }

        kmer_send_buffer.clear();

        ////////////////////////////////////////////////////////////////////

        {
          // TRANSFER and reset //
          transfer_kmers (scounts_kmer, kmer_send_buffer_final);
          num_kmers = 0;
          kmer_send_buffer.clear();
          kmer_send_buffer_final.clear();
        }
      }

  }
  *n_kmers = 0;
}


// Parallel read and sort 
// TODO(vamsikku): will only run in 1 batch //
void Sliding_window_concurrent_v1(const char *ptr,
      size_t length, size_t *n_kmers, 
      std::vector<std::vector<kmer_t>> &partial_kmer_counts) {
  size_t p=0;
  std::vector<int> scounts_kmer (size,0);
  std::vector<KmerPairs> kmer_pair_buffer; //k_count as hash place holder//
  std::vector<KmerPairs> kmer_send_buffer;
  size_t kpos=0;
  size_t num_kmers;
  {
      const char *symbol_begin_ = ptr;
      const char *symbol_end_ = ptr + length;
      const char skip_symbol_ = '>';
      const char delimiter_ = '\n';
      size_t read_len_ = 100UL;
      size_t klen_ = (size_t) KMER_LENGTH;
      // Find the k-mers a thread would process //
      size_t tcount = 0;
#pragma omp parallel
      {
        if (!omp_get_thread_num())
          tcount = omp_get_num_threads();
      }

      std::vector<size_t> threads_part_count(tcount, 0UL);
      std::vector<std::pair<const char *, const char *> > 
          threads_begin_end(tcount); 

      // STEP-1 find part counts of each thread in parallel //
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
        }
      }

      size_t cum = threads_part_count[0];
      for (size_t i=1; i<tcount; i++) {
        cum += threads_part_count[i];
        threads_part_count[i] = cum;
      }

      // processor i (i > 0) 
      // start = threads_part_count[i-1], end = threads_part_count[i] //
      // threads_part_count[-1] = 0 //

      kmer_pair_buffer.resize(cum);
      num_kmers = cum;
      // every thread extracts k-mers in parallel and fill kmer_partitions//
#pragma omp parallel
      {
        size_t my_id = omp_get_thread_num();
        size_t start_idx = my_id ? threads_part_count[my_id-1] : 0UL;;
        size_t end_idx = threads_part_count[my_id];
        const auto &my_begin_end = threads_begin_end[my_id];

        const char *beg = my_begin_end.first;
        const char *end = my_begin_end.second;

        while (beg < end) {
          for (; (beg<end) && (*beg != '>'); ++beg) { }
          for(; (beg<end) && (*beg !='\n'); beg++) {/*noop*/}
          if ((beg+KMER_LENGTH) >= end) { break; }

          ++beg;

          kmer_t kmer=0;
          for (int i=0; (*beg != '\n') && (i<KMER_LENGTH-1); ++i, ++beg) {
            kmer = kmer_shift(kmer, char_to_el(*beg));
          }

          while ((beg < end) && (*beg != '\n')) {
            kmer = kmer_shift(kmer, char_to_el(*beg));
            kmer_pair_buffer[ start_idx++ ] = 
                KmerPairs{kmer, retrieve_proc_id(kmer)};
            ++beg;
          }
          ++beg;
        }

        if (start_idx != end_idx) {
          throw std::logic_error("Invalid partition size");
        }
      }

      {

        KmerPairCompareHashed kmer_pair_compare;
        oneapi::tbb::parallel_sort(kmer_pair_buffer.begin(), 
              kmer_pair_buffer.begin()+num_kmers, kmer_pair_compare);
        
        std::atomic<int> partition_counts[size] = { {0} };
        std::vector<size_t> uniq_kmer_counts(tcount, 0UL);
        std::vector<size_t> kmer_starts(tcount, 0UL);

        // parallel duplicate removal, frequency count and partition count 
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




#pragma omp parallel for
        for (size_t i=0; i<size; i++) {
          scounts_kmer[i] = int(partition_counts[i]);
        }

        // prefix sum //
        // TODO(vamsikku): 
        std::vector<size_t> prefix_sum(tcount, 0UL);
        size_t cum = uniq_kmer_counts[0UL]; 
        prefix_sum[0UL] = 0UL;
        for (size_t i=1; i<tcount; i++) {
          prefix_sum[i] = cum;
          cum += uniq_kmer_counts[i];
        }
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



        {
          // TRANSFER and reset //
          transfer_kmers (scounts_kmer, kmer_send_buffer);
          num_kmers = 0;
          kmer_send_buffer.clear();
        }
      }
  }
  *n_kmers = 0;
}



void Sliding_window_concurrent_v0(const char *ptr, size_t length, size_t *n_kmers, 
                     std::vector<std::vector<kmer_t>> &partial_kmer_counts) {
  size_t p=0;
  std::vector<int> scounts_kmer (size,0);
  std::vector< std::vector<KmerPairs> > kmer_partitions(size);
  std::vector<KmerPairs> kmer_send_buf;
  size_t kpos=0;

  /*find start of a read*/
  for(; ptr[p]!='>' && p<length; p++) {/*noop*/ }

  size_t num_kmers=*n_kmers;
  kmer_t kmer = 0; 
  int min_pos=0, tmp_pos=0;

  while(p<length) {
    assert(ptr[p]=='>'); /*this will be true*/

    /*skip till newline*/
    for(; p<length && ptr[p]!='\n'; p++) {/*noop*/ }
    p++; /*skip the newline*/

    if(p+KMER_LENGTH > length) break; /*too short a read*/

    kmer=0;
    int i;

    for(i=0; ptr[p]!='\n' && i<KMER_LENGTH-1; i++) {
      kmer = kmer_shift(kmer, char_to_el(ptr[p]));
      p++;
    }

    while(p<length && ptr[p]!='\n') {
      //double TM3 = MPI_Wtime();
      kmer = kmer_shift(kmer, char_to_el(ptr[p]));
   
      p++;
      //double TM4 = MPI_Wtime();
      //kmer_shift_time += (TM4-TM3);

      double T1 = MPI_Wtime();
      kmer_partitions[retrieve_proc_id(kmer)].push_back(KmerPairs{kmer, 1});
      double T2 = MPI_Wtime();
      vec_insert_time += (T2-T1);

      num_kmers++;

      if (num_kmers > MAX_KMER_COUNT){
          double t25 = MPI_Wtime ();
          double T5 = MPI_Wtime();
          {

            KmerPairCompare kmer_pair_compare;

#pragma omp parallel for
            for (int t=0; t<size; t++) {
              oneapi::tbb::parallel_sort(kmer_partitions[t].begin(),
                    kmer_partitions[t].end(), kmer_pair_compare);
              scounts_kmer[t] = 0;
              if (!kmer_partitions[t].empty()) {
                auto& sorted_kmers = kmer_partitions[t];
                int curr_freq = 1;
                size_t uniq_idx = 0UL;

                for (size_t i=1; i<sorted_kmers.size(); i++) {
                  if (sorted_kmers[i].seq == sorted_kmers[uniq_idx].seq) {
                    ++curr_freq;
                  } else {
                    sorted_kmers[uniq_idx].k_count = curr_freq;
                    sorted_kmers[++uniq_idx] = sorted_kmers[i];
                    curr_freq = 1;
                  }
                }
                sorted_kmers[uniq_idx].k_count = curr_freq;
                scounts_kmer[t] = ++uniq_idx;
              }
            }

            std::vector<size_t> kmer_starts(size, 0UL);
            kmer_starts[0] = 0UL;
            size_t cum = scounts_kmer[0];
            for (size_t i=1; i<size; i++) {
              kmer_starts[i] = cum;
              cum += scounts_kmer[i];
            }

            kmer_send_buf.resize(cum);

#pragma omp parallel for
            for (int t=0; t<size; t++) {
              size_t beg = kmer_starts[t];
#pragma omp parallel for
              for (int i=0; i<scounts_kmer[t]; i++) {
                kmer_send_buf[beg+i] = kmer_partitions[t][i];
              }
            }
          }

          double T6 = MPI_Wtime();
          tmap_insert_time += (T6-T5);
          transfer_kmers (scounts_kmer, kmer_send_buf);
          num_kmers = 0;
          for (int t=0; t<size; t++) {
            kmer_partitions[t].clear();
            scounts_kmer[t] = 0;
          }
          kmer_send_buf.clear();
          double t26 = MPI_Wtime ();
          sl_win_time_int += (t26-t25);
       } // end of if condition
    } // end of while loop
    p++; /*skip the newline*/
  }

  if (num_kmers > 0 ) {
      {
          double t25 = MPI_Wtime ();
          double T5 = MPI_Wtime();
          {

            KmerPairCompare kmer_pair_compare;

#pragma omp parallel for
            for (int t=0; t<size; t++) {
              oneapi::tbb::parallel_sort(kmer_partitions[t].begin(),
                    kmer_partitions[t].end(), kmer_pair_compare);
              if (!kmer_partitions[t].empty()) {
                auto& sorted_kmers = kmer_partitions[t];
                int curr_freq = 1;
                size_t uniq_idx = 0UL;

                for (size_t i=1; i<sorted_kmers.size(); i++) {
                  if (sorted_kmers[i].seq == sorted_kmers[uniq_idx].seq) {
                    ++curr_freq;
                  } else {
                    sorted_kmers[uniq_idx].k_count = curr_freq;
                    sorted_kmers[++uniq_idx] = sorted_kmers[i];
                    curr_freq = 1;
                  }
                }
                sorted_kmers[uniq_idx].k_count = curr_freq;
                scounts_kmer[t] = ++uniq_idx;
              }
            }

            std::vector<size_t> kmer_starts(size, 0UL);
            kmer_starts[0] = 0UL;
            size_t cum = scounts_kmer[0];
            for (size_t i=1; i<size; i++) {
              kmer_starts[i] = cum;
              cum += scounts_kmer[i];
            }

            kmer_send_buf.resize(cum);

#pragma omp parallel for
            for (int t=0; t<size; t++) {
              size_t beg = kmer_starts[t];
#pragma omp parallel for
              for (int i=0; i<scounts_kmer[t]; i++) {
                kmer_send_buf[beg+i] = kmer_partitions[t][i];
              }
            }
          }

          double T6 = MPI_Wtime();
          tmap_insert_time += (T6-T5);
          transfer_kmers (scounts_kmer, kmer_send_buf);
          num_kmers = 0;
          scounts_kmer.clear(); 
          double t26 = MPI_Wtime ();
          sl_win_time_int += (t26-t25);
       } // end of if condition
  }

  *n_kmers = num_kmers;
}

// Algo: 
// 1. since the MAX_KMER_COUNT is fixed used it as a fixed buffer to fill 
// the kmers (current not done in parallel) //
// 2. sort first on parition and then on actual kmer
// 3. parallel remove duplicates inplace
// 4. parallel write to kmer_send_buffer to be consumed by AlltoAllV 
void Sliding_window_concurrent_v1_old(const char *ptr, size_t length, 
      size_t *n_kmers, std::vector<std::vector<kmer_t>> &partial_kmer_counts) {
  size_t p=0;
  std::vector<int> scounts_kmer(size,0);
  std::vector<KmerPairs> kmer_pair_buffer; //k_count as hash place holder//
  std::vector<KmerPairs> kmer_send_buffer;

  //kmer_pair_buffer.resize(MAX_KMER_COUNT);

  for(; ptr[p]!='>' && p<length; p++) {/*noop*/ }

  kmer_t kmer = 0; 
  size_t num_kmers = 0;
  bool final_sweep = false;
  while(p<length) {
    assert(ptr[p]=='>'); /*this will be true*/

    /*skip till newline*/
    for(; p<length && ptr[p]!='\n'; p++) {/*noop*/ }
    p++; /*skip the newline*/

    if(p+KMER_LENGTH > length) break; /*too short a read*/

    kmer=0;
    int i;

    for(i=0; ptr[p]!='\n' && i<KMER_LENGTH-1; i++) {
      kmer = kmer_shift(kmer, char_to_el(ptr[p]));
      p++;
    }

    while(p<length && ptr[p]!='\n') {
      kmer = kmer_shift(kmer, char_to_el(ptr[p]));
      p++;

      if (num_kmers >= MAX_KMER_COUNT) {
        throw std::logic_error("num_kmers invariant violation");
      }

      {
        double T1 = MPI_Wtime();
        kmer_pair_buffer.push_back(KmerPairs{kmer, retrieve_proc_id(kmer)});
        ++num_kmers;
      //  kmer_pair_buffer[num_kmers++] = KmerPairs{kmer, retrieve_proc_id(kmer)};
        double T2 = MPI_Wtime();
        vec_insert_time += (T2-T1);
      }

      if (num_kmers >= MAX_KMER_COUNT) {
CONCURRENT_PROCESS_KMERS:
        double t25 = MPI_Wtime ();
        double T5 = MPI_Wtime();
        {
          KmerPairCompareHashed kmer_pair_compare;
          oneapi::tbb::parallel_sort(kmer_pair_buffer.begin(), 
              kmer_pair_buffer.begin()+num_kmers, kmer_pair_compare);

          ////////////////////////////////////////////////////////////////////
          size_t tcount = 1;
#pragma omp parallel
          { //TODO(vamsikku): why do we need a pragma to get # of threads // 
            if (!omp_get_thread_num()) {
              tcount = omp_get_num_threads();
            }
          }
          
          std::atomic<size_t> partition_counts[size] = { {0} };
          std::vector<size_t> uniq_kmer_counts(tcount, 0UL);
          std::vector<size_t> kmer_starts(tcount, 0UL);

          // parallel duplicate removal, frequency count and partition count 
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

#pragma omp parallel for
          for (size_t i=0; i<size; i++) {
            scounts_kmer[i] = int(partition_counts[i]);
          }

          // prefix sum //
          // TODO(vamsikku): 
          std::vector<size_t> prefix_sum(tcount, 0UL);
          size_t cum = uniq_kmer_counts[0UL]; 
          prefix_sum[0UL] = 0UL;
          for (size_t i=1; i<tcount; i++) {
            prefix_sum[i] = cum;
            cum += uniq_kmer_counts[i];
          }
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
        }

        {
          // TRANSFER and reset //
          double T6 = MPI_Wtime();
          tmap_insert_time += (T6-T5);
          transfer_kmers (scounts_kmer, kmer_send_buffer);
          num_kmers = 0;
          kmer_send_buffer.clear();
          double t26 = MPI_Wtime ();
          sl_win_time_int += (t26-t25);
        }

        if (final_sweep) { goto DONE; }
      } // end of if condition

    } // end of while loop
    p++; /*skip the newline*/
  }

  if (num_kmers) {
    final_sweep = true;
    goto CONCURRENT_PROCESS_KMERS;
  } 

DONE:
  *n_kmers = num_kmers;
}



void Sliding_window_concurrent_old(
      const char *ptr, size_t length, size_t *n_kmers, 
      std::vector<std::vector<kmer_t>> &partial_kmer_counts) {
  double t25 = MPI_Wtime ();
  size_t tcount = 0;
#pragma omp parallel
  {
    if (!omp_get_thread_num())
      tcount = omp_get_num_threads();
  }

  //NOTE: find unique k-mers 
  std::vector<int> scounts_kmer (size,0);
  std::vector<KmerPairs> kmer_send_buf;

  const char *symbol_begin_ = ptr;
  const char *symbol_end_ = ptr + length;
  size_t klen_ = 32UL; // TODO(vamsikku): make this a parameter //
  size_t read_len_ = 100UL; //TOD0(vamsikku): make this a parameter //
  char skip_symbol_ = '>';
  char delimiter_ = '\n';
 
  /////////////////////////////////// NEW ALGO /////////////////////////////////
  kmer_strncmp_t kmer_strncmp;
  std::vector< std::pair<const char *, size_t> > ordered_permutation;
  std::vector<size_t> threads_part_count(tcount, 0UL);
  std::vector<std::pair<const char *, const char *> > threads_begin_end(tcount);

  std::atomic<size_t> read_count = {0UL};
  std::atomic<size_t> ordered_permutation_length = {0UL};
  

  //STEP-1: determine the pre-allocation size // 
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
      }

      //TODO(vamsikku): filling the kmer-send buffer can also be done
      //concurrently
      size_t total_uniq_kmer_count = 0;
      for (size_t i=0; i<tcount; i++) {
        total_uniq_kmer_count += uniq_kmer_counts[i]; 
      }
      printf("unique kmers=%lu\n", total_uniq_kmer_count);


#if 0
      std::vector<std::pair<KmerPairs, size_t> > temp_buffer;
      kmer_send_buf.resize(total_uniq_kmer_count);
      temp_buffer.resize(total_uniq_kmer_count);

      size_t temp_buffer_idx=0, dest_id;
      kmer_t kmer;
      KmerPairs kmer_pair;

      for (auto itr=ordered_permutation.begin(); 
            itr!=ordered_permutation.end(); ++itr) {
        // extract the k-mer //
        kmer=0; 
        char const * ptr = itr->first;
        for (size_t i=0; i<KMER_LENGTH-1; i++, ++ptr) {
          kmer = kmer_shift(kmer, char_to_el(*ptr));
        }
        kmer_pair.seq = kmer;
        kmer_pair.k_count = (int) itr->second;
        dest_id = retrieve_proc_id(kmer);
        temp_buffer[ temp_buffer_idx++ ] = std::make_pair(kmer_pair, dest_id);
        scounts_kmer[ dest_id ]++;
      }

      std::vector<int> scounts_idx(size, 0);
      int cum_scounts = scounts_kmer[0];
      for (size_t i=1; i<size; i++) { cum_scounts += scounts_kmer[i]; }

      for (size_t i=0; i<temp_buffer.size(); i++) {
        dest_id = temp_buffer[i].second;
        kmer_send_buf[ scounts_idx[dest_id]++ ] = temp_buffer[i].first;
      }

      transfer_kmers (scounts_kmer, kmer_send_buf);
      scounts_kmer.clear(); 
      double t26 = MPI_Wtime ();
      sl_win_time_int += (t26-t25);
#endif
}




void Sliding_window_stripped (const char *ptr, size_t length, size_t *n_kmers, 
                     std::vector<std::vector<kmer_t>> &partial_kmer_counts) {
  size_t p=0;
  std::vector<int> scounts_kmer (size,0);
  std::vector<KmerPairs> kmer_send_buf;
  size_t kpos=0;

  /*find start of a read*/
  for(; ptr[p]!='>' && p<length; p++) {/*noop*/ }

  size_t num_kmers=*n_kmers;
  kmer_t kmer = 0; 
  lmer_t lmer_out = 0, min_lmer=0, tmp_lmer=0;
  uint64_t min_lmer_freq=0, tmp_lmer_freq=0;
  int min_pos=0, tmp_pos=0;

  while(p<length) {
    assert(ptr[p]=='>'); /*this will be true*/

    /*skip till newline*/
    for(; p<length && ptr[p]!='\n'; p++) {/*noop*/ }
    p++; /*skip the newline*/

    if(p+KMER_LENGTH > length) break; /*too short a read*/

    kmer=0, lmer_out=0;
    min_lmer=0, min_lmer_freq=0;
    tmp_lmer=0, tmp_lmer_freq=0;
    min_pos=0, tmp_pos=0;
    int i;

    for(i=0; ptr[p]!='\n' && i<KMER_LENGTH-1; i++) {
      kmer = kmer_shift(kmer, char_to_el(ptr[p]));

      if (i<LMER_LENGTH-1) 
          lmer_out = lmer_shift(lmer_out, char_to_el(ptr[p]));
      else {
           lmer_out = lmer_shift(lmer_out, char_to_el(ptr[p]));

           tmp_lmer = lmer_out;
           tmp_lmer_freq = global_lmer_frequency[tmp_lmer];
           tmp_pos = i-(LMER_LENGTH-1);

           if (i == LMER_LENGTH-1) {
               min_lmer = tmp_lmer;
               min_lmer_freq = tmp_lmer_freq;
           }
           else {
                if (tmp_lmer_freq < min_lmer_freq) {
                    min_lmer = tmp_lmer;
                    min_lmer_freq = tmp_lmer_freq;
                    min_pos = tmp_pos;
                }
           }
      }
      p++;
    }

    while(p<length && ptr[p]!='\n') {
      //double TM3 = MPI_Wtime();
      kmer = kmer_shift(kmer, char_to_el(ptr[p]));
      lmer_out = lmer_shift(lmer_out, char_to_el(ptr[p]));
      uint64_t lmer_out_freq = global_lmer_frequency[lmer_out];
   
      //double TM5 = MPI_Wtime();
      if (min_pos < 0) {
          recalculate_min_lmer(kmer, &min_lmer, &min_lmer_freq, &min_pos);
          num_recalculate_lmer++;
      }
      //double TM6 = MPI_Wtime();
      //kmer_recal_time += (TM6-TM5);

      if (lmer_out_freq < min_lmer_freq) {
          min_lmer = lmer_out;
          min_lmer_freq = lmer_out_freq;
          min_pos = KMER_LENGTH-LMER_LENGTH;
          //min_lmer_freq = global_lmer_frequency[lmer_out];
      }
      p++;
      min_pos--;
      //double TM4 = MPI_Wtime();
      //kmer_shift_time += (TM4-TM3);

      double T1 = MPI_Wtime();
      partial_kmer_counts[retrieve_proc_id(min_lmer)].push_back(kmer);
      double T2 = MPI_Wtime();
      vec_insert_time += (T2-T1);

      num_kmers++;

      if (num_kmers > MAX_KMER_COUNT){
          double t25 = MPI_Wtime ();
          double T5 = MPI_Wtime();
          {
            for (int t=0; t<size; t++)
            {
              int counter=1; kpos=0;
              sort(partial_kmer_counts[t].begin(), partial_kmer_counts[t].end());
              kmer_t prev=partial_kmer_counts[t][0];
              for(size_t i = 1; i < (partial_kmer_counts[t].size()); i++)
              {
                if (partial_kmer_counts[t][i] == prev) { counter++; }
                else {
                  kmer_send_buf.push_back(KmerPairs{prev, counter});
                  kpos++;
                  counter=1;
                  prev=partial_kmer_counts[t][i];
                }
              }
              kmer_send_buf.push_back(KmerPairs{prev, counter});
              kpos++;

              scounts_kmer[t] = kpos;
              partial_kmer_counts[t].clear();
              partial_kmer_counts[t].shrink_to_fit();
            }
          }

          double T6 = MPI_Wtime();
          tmap_insert_time += (T6-T5);
          transfer_kmers (scounts_kmer, kmer_send_buf);
          num_kmers = 0;
          scounts_kmer.clear(); 
          double t26 = MPI_Wtime ();
          sl_win_time_int += (t26-t25);
       } // end of if condition
    } // end of while loop
    p++; /*skip the newline*/
  }
  *n_kmers = num_kmers;
  scounts_kmer.shrink_to_fit();
}



void Sliding_window (const char *ptr, size_t length, size_t *n_kmers, 
                     std::vector<std::vector<kmer_t>> &partial_kmer_counts)
{
#ifdef LMER_DEBUG
    FILE *fp_d;
    char debug_file_name[25];
    char proc_id[3];

    sprintf(proc_id, "%d", rank); 
    strcpy(debug_file_name,"debug_p");
    strcpy(&debug_file_name[strlen(debug_file_name)],proc_id);
    strcpy(&debug_file_name[strlen(debug_file_name)],".log");
    fp_d = fopen (debug_file_name, "w");
    
    /*check to see if it opened okay */
    if (fp_d == NULL)
    {
		printf ("Error opening proc %d 's dump file \n", rank);
		exit (0);
    }
#endif

  size_t p=0;
  std::vector<int> scounts_kmer (size,0);
  std::vector<KmerPairs> kmer_send_buf;
  size_t kpos=0;

  //char kmer_out[WINDW_SIZE_PLUS];

  /*find start of a read*/
  for(; ptr[p]!='>' && p<length; p++) {/*noop*/ }

  size_t num_kmers=*n_kmers;
  kmer_t kmer = 0; 
  lmer_t lmer_out = 0, min_lmer=0, tmp_lmer=0;
  uint64_t min_lmer_freq=0, tmp_lmer_freq=0;
  int min_pos=0, tmp_pos=0;

  while(p<length) {
    assert(ptr[p]=='>'); /*this will be true*/

    /*skip till newline*/
    for(; p<length && ptr[p]!='\n'; p++) {/*noop*/ }
    p++; /*skip the newline*/

    if(p+KMER_LENGTH > length) break; /*too short a read*/

    kmer=0, lmer_out=0;
    min_lmer=0, min_lmer_freq=0;
    tmp_lmer=0, tmp_lmer_freq=0;
    min_pos=0, tmp_pos=0;
    int i;

    //double TM1 = MPI_Wtime();
    for(i=0; ptr[p]!='\n' && i<KMER_LENGTH-1; i++) {
      kmer = kmer_shift(kmer, char_to_el(ptr[p]));

      if (i<LMER_LENGTH-1) 
          lmer_out = lmer_shift(lmer_out, char_to_el(ptr[p]));
      else {
           lmer_out = lmer_shift(lmer_out, char_to_el(ptr[p]));

           tmp_lmer = lmer_out;
           tmp_lmer_freq = global_lmer_frequency[tmp_lmer];
           tmp_pos = i-(LMER_LENGTH-1);

           if (i == LMER_LENGTH-1) {
               min_lmer = tmp_lmer;
               min_lmer_freq = tmp_lmer_freq;
           }
           else {
                if (tmp_lmer_freq < min_lmer_freq) {
                    min_lmer = tmp_lmer;
                    min_lmer_freq = tmp_lmer_freq;
                    min_pos = tmp_pos;
                }
           }
      }
      p++;
    }
    //double TM2 = MPI_Wtime();
    //kmer_shift_time += (TM2-TM1);

    while(p<length && ptr[p]!='\n') {
      //double TM3 = MPI_Wtime();
      kmer = kmer_shift(kmer, char_to_el(ptr[p]));
      lmer_out = lmer_shift(lmer_out, char_to_el(ptr[p]));
      uint64_t lmer_out_freq = global_lmer_frequency[lmer_out];
   
      //double TM5 = MPI_Wtime();
      if (min_pos < 0) {
          recalculate_min_lmer(kmer, &min_lmer, &min_lmer_freq, &min_pos);
          num_recalculate_lmer++;
      }
      //double TM6 = MPI_Wtime();
      //kmer_recal_time += (TM6-TM5);

      if (lmer_out_freq < min_lmer_freq) {
          min_lmer = lmer_out;
          min_lmer_freq = lmer_out_freq;
          min_pos = KMER_LENGTH-LMER_LENGTH;
          //min_lmer_freq = global_lmer_frequency[lmer_out];
      }
      p++;
      min_pos--;
      //double TM4 = MPI_Wtime();
      //kmer_shift_time += (TM4-TM3);

#ifdef LMER_DEBUG
      fprintf(fp_d, "kmer: %lu, lmer: %lu, lmer_freq: %lu, min_lmer: %lu, min_freq: %lu\n", kmer,lmer_out,lmer_out_freq,min_lmer,min_lmer_freq);
#endif
#if 0
      /* Murmerhash method */
      /*
      calculate owner using MurmurHash and update corresponding scounts_kmer
      got = tmp_k_map.find(kmer);
      if ( got == tmp_k_map.end() ) 
          scounts_kmer[(int) MurmurHash64A ((char*)&kmer, sizeof(kmer), SEED, size)]++;

      tmp_k_map[kmer]++;
      */

#endif
      double T1 = MPI_Wtime();
#ifdef NOOPT
      kmers_per_proc[retrieve_proc_id(min_lmer)].push_back(kmer);
#endif
      partial_kmer_counts[retrieve_proc_id(min_lmer)].push_back(kmer);
      double T2 = MPI_Wtime();
      vec_insert_time += (T2-T1);

      num_kmers++;

      if (num_kmers > MAX_KMER_COUNT){
          double t25 = MPI_Wtime ();
          // initiate collective communication to pass k-mers and their respective counts to rightful owners
          // calculate global owner of each k-mer and populate the k-mer and count to their respective 'p'th local buffers
          // if global position of a k-mer in my k_map != me, delete k-mer from my k_map
          // iterate through k-mers recieved after collective communication ends, and add k-mers to my k_map 
          // reset num_kmers count to 0.
          
          double T5 = MPI_Wtime();
#ifdef NOOPT
          for (int t=0; t<size; t++)
          {
		  int counter=1;
                  sort(kmers_per_proc[t].begin(), kmers_per_proc[t].end());
                  kmer_t prev=kmers_per_proc[t][0];
		  for(size_t i = 1; i < (kmers_per_proc[t].size()); i++)
		  {
			 if (kmers_per_proc[t][i] == prev) {
			     counter++;
                         } else {
			     kmer_cnt_tmp_buf[t].push_back(counter);
			     counter=1;
                             prev=kmers_per_proc[t][i];
			 } 
		     
		  }
                  kmer_cnt_tmp_buf[t].push_back(counter);

                  kmers_per_proc[t].erase( unique( kmers_per_proc[t].begin(), kmers_per_proc[t].end() ), kmers_per_proc[t].end() );
                  assert(kmers_per_proc[t].size() == kmer_cnt_tmp_buf[t].size());
                  scounts_kmer[t] = kmer_cnt_tmp_buf[t].size(); 
                  //scounts_kmer[t] = std::accumulate(kmer_cnt_tmp_buf[t].begin(), kmer_cnt_tmp_buf[t].end(), 0); 
          }
#else
          for (int t=0; t<size; t++)
          {
		  int counter=1;
                  kpos=0;
                  sort(partial_kmer_counts[t].begin(), partial_kmer_counts[t].end());
                  kmer_t prev=partial_kmer_counts[t][0];
		  for(size_t i = 1; i < (partial_kmer_counts[t].size()); i++)
		  {
			 if (partial_kmer_counts[t][i] == prev) {
			     counter++;
                         } else {
                             kmer_send_buf.push_back(KmerPairs{prev, counter});
                             //kmers_per_proc[t].push_back(KmerPairs{prev, counter});
			     kpos++;
                             counter=1;
                             prev=partial_kmer_counts[t][i];
			 }
		     
		  }
                  //kmers_per_proc[t].push_back(KmerPairs{prev, counter});
                  kmer_send_buf.push_back(KmerPairs{prev, counter});
                  kpos++;

                  //kmers_per_proc[t].erase( unique( kmers_per_proc[t].begin(), kmers_per_proc[t].end() ), kmers_per_proc[t].end() );
                  //assert(kmers_per_proc[t].size() == kmer_cnt_tmp_buf[t].size());
                  //scounts_kmer[t] = kmers_per_proc[t].size();
                  scounts_kmer[t] = kpos;
                  partial_kmer_counts[t].clear();
                  partial_kmer_counts[t].shrink_to_fit();
          }
          /*
          for (int t=0; t<size; t++)
          {
                  sort(partial_kmer_counts[t].begin(), partial_kmer_counts[t].end());
                  std::vector<kmer_t>::iterator low,up;
		  for(auto it=partial_kmer_counts[t].begin(); it != partial_kmer_counts[t].end(); )
		  {
                      low=std::lower_bound (partial_kmer_counts[t].begin(), partial_kmer_counts[t].end(), *it);
                      up= std::upper_bound (partial_kmer_counts[t].begin(), partial_kmer_counts[t].end(), *it);
                      int ncount = std::count (low, up, *it);
                      kmers_per_proc[t].push_back(KmerPairs{*it, ncount});
                      it = up;
		  }
                  scounts_kmer[t] = kmers_per_proc[t].size(); 
          }
          */
#endif
          double T6 = MPI_Wtime();
          tmap_insert_time += (T6-T5);

          /*double T5 = MPI_Wtime();
          for(std::vector<kmer_t>::iterator it = kmer_tmp_buf.begin(); it != kmer_tmp_buf.end(); ++it) {
              if (tmp_k_map[*it]++ == 0)
                  scounts_kmer[(int) MurmurHash64A ((char*)&(*it), sizeof(*it), SEED, size)]++;
          }
          double T6 = MPI_Wtime();
          tmap_insert_time += (T6-T5);
          */

          //clear the partial counts
          /*
          for(int k=0; k<size; k++)
              partial_kmer_counts[k].clear();         
          */

#ifdef NOOPT
          transfer_kmers (scounts_kmer, kmers_per_proc, kmer_cnt_tmp_buf);
#else
          transfer_kmers (scounts_kmer, kmer_send_buf);
#endif
          num_kmers = 0;

#ifdef NOOPT
          for (int t=0; t<size; t++) {
               kmers_per_proc[t].clear();
               kmer_cnt_tmp_buf[t].clear();
          }
#endif
          //memset(scounts_kmer, 0, size*sizeof(*scounts_kmer));
          scounts_kmer.clear(); 
          double t26 = MPI_Wtime ();
          sl_win_time_int += (t26-t25);
       } // end of if condition
    } // end of while loop
#ifdef LMER_DEBUG
    fprintf(fp_d, "----------------------------\n");
#endif
    p++; /*skip the newline*/
  }

#ifdef LMER_DEBUG
          fclose (fp_d);
#endif

  *n_kmers = num_kmers;
  scounts_kmer.shrink_to_fit();

}

void process_remaining_kmers(
                     std::vector<std::vector<kmer_t>> &partial_kmer_counts) 
{
    if (partial_kmer_counts.empty()) { return; }

     std::vector<int> scounts_kmer (size,0);
     std::vector<KmerPairs> kmer_send_buf;
     size_t kpos=0;

      double T7 = MPI_Wtime();
          for (int t=0; t<size; t++)
          {
             if(partial_kmer_counts[t].size())
             {
		  int counter=1;
                  kpos=0;
                  sort(partial_kmer_counts[t].begin(), partial_kmer_counts[t].end());
                  kmer_t prev=partial_kmer_counts[t][0];
		  for(size_t i = 1; i < (partial_kmer_counts[t].size()); i++)
		  {
			 if (partial_kmer_counts[t][i] == prev) {
			     counter++;
                         } else {
                             kmer_send_buf.push_back(KmerPairs{prev, counter});
			     kpos++;
			     counter=1;
                             prev=partial_kmer_counts[t][i];
			 }
		     
		  }
                  kmer_send_buf.push_back(KmerPairs{prev, counter});
                  kpos++;

                  //kmers_per_proc[t].erase( unique( kmers_per_proc[t].begin(), kmers_per_proc[t].end() ), kmers_per_proc[t].end() );
                  //assert(kmers_per_proc[t].size() == kmer_cnt_tmp_buf[t].size());
                  //scounts_kmer[t] = kmers_per_proc[t].size(); 
                  scounts_kmer[t] = kpos;
                  partial_kmer_counts[t].clear();
                  partial_kmer_counts[t].shrink_to_fit();
              }
           }
      /*
      for(int t=0; t<size; t++)
      {
           if(partial_kmer_counts[t].size())
           {
                  sort(partial_kmer_counts[t].begin(), partial_kmer_counts[t].end());
                  std::vector<kmer_t>::iterator low,up;
		  for(auto it=partial_kmer_counts[t].begin(); it != partial_kmer_counts[t].end(); )
		  {
                      low=std::lower_bound (partial_kmer_counts[t].begin(), partial_kmer_counts[t].end(), *it);
                      up= std::upper_bound (partial_kmer_counts[t].begin(), partial_kmer_counts[t].end(), *it);
                      int ncount = std::count (low, up, *it);
                      kmers_per_proc[t].push_back(KmerPairs{*it, ncount});
                      it = up;
		  }
                  scounts_kmer[t] = kmers_per_proc[t].size(); 

             
             // if (kmers_per_proc[t].size() != kmer_cnt_tmp_buf[t].size())
             // {
             //     fprintf(stderr, "counters not matching, for rank: %d, t: %d!! kmers_per_proc size: %lu, kmer_cnt_tmp_buf: %lu\n",
             //             rank, t, kmers_per_proc[t].size(), kmer_cnt_tmp_buf[t].size());

             //     
             //     char pros_id[3];
             //     char outfile_name[25];

             //     sprintf(pros_id, "%d", rank); 
             //     strcpy(outfile_name,"tmp_kmers_");
             //     strcpy(&outfile_name[strlen(outfile_name)],pros_id);
             //     strcpy(&outfile_name[strlen(outfile_name)],".log");
             //     FILE *ft = fopen(outfile_name, "w");
             //     if (ft == NULL)
             //     {
             //        printf("Error opening tmp file!\n");
             //        exit(1);
             //     }
             //     
             //     for (size_t m=0; m<kmers_per_proc[t].size(); m++)
             //          fprintf(ft, "m: %lu, kmer: %lu\n", m, kmers_per_proc[t][m]);
             //     fprintf(ft, "----------- end of kmers ------------\n");
             //     for (size_t m=0; m<kmer_cnt_tmp_buf[t].size(); m++)
             //          fprintf(ft, "m: %lu, cnt: %d\n", m, kmer_cnt_tmp_buf[t][m]);

             //     fclose(ft);

             // }
             // assert(kmers_per_proc[t].size() == kmer_cnt_tmp_buf[t].size());
             // scounts_kmer[t] = kmer_cnt_tmp_buf[t].size();
               
           }
          }*/
          double T8 = MPI_Wtime();
          tmap_insert_time += (T8-T7);

          //clear the partial counts
          for(int k=0; k<size; k++)
              partial_kmer_counts[k].clear();         
 
          transfer_kmers (scounts_kmer, kmer_send_buf);
         
          scounts_kmer.clear();
          scounts_kmer.shrink_to_fit();

}

void print_kmer_count_timers()
{

    MPI_Reduce(&alltoall_time, &global_alltoall_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) printf ("Average time for Alltoall across all procs (secs): %f \n", 
                            (double)global_alltoall_time/(double)size);

    MPI_Reduce(&pack_sbuf_time, &global_pack_sbuf_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) printf ("Average time for pack_sbuf_time across all procs (secs): %f \n", 
                            (double)global_pack_sbuf_time/(double)size);

    MPI_Reduce(&alltoallv_time, &global_alltoallv_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) printf ("Average time for AlltoallV across all procs (secs): %f \n", 
                            (double)global_alltoallv_time/(double)size);

    MPI_Reduce(&unpack_rbuf_time, &global_unpack_rbuf_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) printf ("Average time for unpack_rbuf_time across all procs (secs): %f \n", 
                            (double)global_unpack_rbuf_time/(double)size);

    MPI_Reduce(&unpack_rbuf_sort, &global_unpack_rbuf_sort, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) printf ("Average time for unpack_rbuf_time:sort across all procs (secs): %f \n", 
                            (double)global_unpack_rbuf_sort/(double)size);

    MPI_Reduce(&unpack_rbuf_insert, &global_unpack_rbuf_insert, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) printf ("Average time for unpack_rbuf_time:insert across all procs (secs): %f \n", 
                            (double)global_unpack_rbuf_insert/(double)size);

    MPI_Reduce(&unpack_rbuf_acc, &global_unpack_rbuf_acc, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) printf ("Average time for unpack_rbuf_time:acc across all procs (secs): %f \n", 
                            (double)global_unpack_rbuf_acc/(double)size);

    MPI_Reduce(&sl_win_time_int, &global_sl_win_time_int, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) printf ("Average time for partial Sliding Window across all procs (secs): %f \n", 
                            (double)global_sl_win_time_int/(double)size);

    MPI_Reduce(&sl_win_time, &global_sl_win_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) printf ("Average time for Sliding Window across all procs (secs): %f \n", 
                            (double)global_sl_win_time/(double)size);

    MPI_Reduce(&sl_lmer_freq, &global_sl_lmer_freq, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) printf ("Average time for lmer Sliding Window across all procs (secs): %f \n", 
                            (double)global_sl_lmer_freq/(double)size);

    MPI_Reduce(&vec_insert_time, &global_vec_insert_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) printf ("Average time for vector insert across all procs (secs): %f \n", 
                            (double)global_vec_insert_time/(double)size);

    MPI_Reduce(&tmap_insert_time, &global_tmap_insert_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) printf ("Average time for temp_map insert across all procs (secs): %f \n", 
                            (double)global_tmap_insert_time/(double)size);

    MPI_Reduce(&num_recalculate_lmer, &global_num_recalculate_lmer, 1, MPI_UINT64_T, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) printf ("Average number of re-calculations of min l-mer across all procs (secs): %f \n", 
                            (double)global_num_recalculate_lmer/(double)size);

    //MPI_Reduce(&kmer_recal_time, &global_kmer_recal_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    //if (rank == 0) printf ("Average time for recalculating min-mers across all procs (secs): %f \n", 
    //                        (double)global_kmer_recal_time/(double)size);

    //MPI_Reduce(&kmer_shift_time, &global_kmer_shift_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    //if (rank == 0) printf ("Average time for k-mer construction across all procs (secs): %f \n", 
    //                        (double)global_kmer_shift_time/(double)size);

    MPI_Reduce(&kmer_count_time, &global_kmer_count_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0) printf ("Average time for performing k-mer counting across all procs (secs): %f \n", 
                            (double)global_kmer_count_time/(double)size);

    uint64_t all_proc_kmer_count = 0;
    uint64_t tmp_kmer_count = kmer_proc_buf.size();
    //assert(kmer_proc_buf.size() == kmer_cnt_proc_buf.size());

    MPI_Reduce(&tmp_kmer_count, &all_proc_kmer_count, 1, MPI_UINT64_T, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0) { 
        printf("Total distinct k-mer entries across all proc's: %lu\n", all_proc_kmer_count); 
        printf("Number of batch iterations: %d \n", num_batch_transfers);
    }

}

void free_kmer_count_buffers()
{
   
    kmer_proc_buf.clear();
    kmer_proc_buf.shrink_to_fit();

}

void perform_kmer_counting (const char *read_data, size_t length)
{

    double start_t = MPI_Wtime ();

    //calculate frequencies of all l-mers in the Read dataset
    double time_l1 = MPI_Wtime ();
    Sliding_window_l(read_data, length);
    double time_l2 = MPI_Wtime ();
    sl_lmer_freq = time_l2 - time_l1;

    // Perform Allreduce to obtain global lmer counts across all proc's
    int num_lmers = pow(4, LMER_LENGTH);
    //global_lmer_frequency.reserve(num_lmers);

    MPI_Allreduce (lmer_frequency.data(), global_lmer_frequency.data(), num_lmers, MPI_UINT64_T, MPI_SUM, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

#ifdef NOOPT
    std::vector< std::vector<kmer_t> > kmers_per_proc; //, std::vector<kmer_t>);
    std::vector< std::vector<int> > kmer_cnt_tmp_buf; //, std::vector<int>);
    for (int i=0; i<size; i++) {
         kmers_per_proc.push_back(std::vector<kmer_t> ());
         kmer_cnt_tmp_buf.push_back(std::vector<int> ());
    }
#else
    std::vector< std::vector<kmer_t> > partial_kmer_counts(size);
#endif
    size_t num_kmers = 0;

    double t23 = MPI_Wtime ();

    if (rank==0)
        fprintf(stderr, "Starting sliding window\n");

    if (USE_CONCURRENT_KMER_ALGO) {
      if (USE_CONCURRENT_KMER_ALGO_VERSION == 1) {
        Sliding_window_concurrent_v1(read_data,
            length, &num_kmers, partial_kmer_counts);
      } else if (USE_CONCURRENT_KMER_ALGO_VERSION == 2) {
        Sliding_window_concurrent_v1_partition_and_sort(read_data, length, 
            &num_kmers, partial_kmer_counts);
      } else if (USE_CONCURRENT_KMER_ALGO_VERSION == 3) {
        Sliding_window_concurrent_v1_partition_and_sort_extra_space(read_data, 
            length, &num_kmers, partial_kmer_counts);
      } else if (USE_CONCURRENT_KMER_ALGO_VERSION == 4 ) {
        Sliding_window_concurrent_v1_implicit(read_data, length, &num_kmers,
            partial_kmer_counts);
      } else if (USE_CONCURRENT_KMER_ALGO_VERSION == 5) {
        intel::esc::Async_Distributed_Kmer_Counter<kmer_t, KmerPairs>
            async_dist_kmer_counter(size, rank, (size_t)KMER_LENGTH,
                  MAX_KMER_COUNT);


        num_batch_transfers += async_dist_kmer_counter.run(read_data, 
              read_data+length, kmer_proc_buf);
        num_kmers = 0;
      } else {
        Sliding_window_concurrent_v0(read_data, length, &num_kmers,
            partial_kmer_counts);
      }
    } else {
      Sliding_window(read_data, length, &num_kmers, partial_kmer_counts);
    }

    double t24 = MPI_Wtime ();
    sl_win_time = t24 - t23;

    if (rank==0)
        fprintf(stderr, "Completed sliding window\n");
        
    // initiate communication for the residual kmers
    //printf("rank: %d, num_kmers: %d\n", rank, num_kmers);
    if (num_kmers) {
         
          if (rank==0) fprintf(stderr, "remaining sliding window, num_kmers: %d\n", num_kmers);
#ifdef NOOPT
          process_remaining_kmers (kmers_per_proc, kmer_cnt_tmp_buf);
#else
          process_remaining_kmers (partial_kmer_counts);
#endif
    }

#ifdef NOOPT
    kmers_per_proc.clear();
    kmers_per_proc.shrink_to_fit();
    kmer_cnt_tmp_buf.clear();
    kmer_cnt_tmp_buf.shrink_to_fit();
#else
    for (int i=0; i<size; i++)
         partial_kmer_counts[i].clear();
    partial_kmer_counts.shrink_to_fit(); 
#endif

    lmer_frequency.clear();
    global_lmer_frequency.clear();

    double end_t = MPI_Wtime ();

    kmer_count_time = end_t - start_t;

    MPI_Barrier(MPI_COMM_WORLD);

    if (rank==0)
        fprintf(stderr, "Completed kmer counting\n");

    print_kmer_count_timers();
}


