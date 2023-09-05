#include <cstdio>
#include <vector>
#include <string>
#include <errno.h>

#include <oneapi/tbb/global_control.h>
#include <concurrent_kmer_counter.hpp>
#include <generate_kmers.hpp>

#include "stop_clock.hpp"


typedef intel::esc::kmer_utils_t<size_t> kmer_utils_t;
typedef intel::esc::Concurrent_Kmer_Partial_Counter<size_t>
    concurrent_kmer_counter_t;

struct count_unique_kmers_t {
  count_unique_kmers_t() : count_(0UL), kmer_freq_(0UL) {}

  template<typename iterator>
  void operator()(iterator beg, iterator end) const {
    kmer_utils_t kutils(32UL);

    for (; beg!=end; ++beg) { ++count_; kmer_freq_+=(size_t)beg->second;
      //std::string kmer_str = kutils.get_kmer_string(beg->first);
      //printf("%s %lu\n", kmer_str.c_str(), beg->second);
    }
  }
  mutable std::atomic<size_t> count_;
  mutable std::atomic<size_t> kmer_freq_;
}; // struct count_unique_kmers_t //

int main(int argc, char **argv) {
  size_t tcount=8;
  bool use_explicit_mode = false;
  bool use_implicit_mode_v0 = false;

  if (argc < 2) {
    printf("Usage: ./concurrent-kmer-table-builder {reads.fasta} [tcount] "
          "[--use_implicit|--use_explicit}]");
    return -1;
  }
  if (argc >= 3) {
    if (sscanf(argv[2], "%lu", &tcount) != 1) {
      tcount = 8;
      printf("WARNING: unable to scan thread count %s defaulting to %lu\n",
          argv[2], tcount);
    }
  }


  if (argc >= 4) {
    use_explicit_mode = (strcmp(argv[3], "--use_explicit") == 0);
    use_implicit_mode_v0 = (strcmp(argv[3], "--use_implicit_v0") == 0);
  }


  omp_set_num_threads(tcount);

  printf("reading %s\n", argv[1]);
  FILE *fptr = fopen(argv[1], "r");
  if (!fptr) {
    throw std::logic_error("Unable to open file: "+std::string(argv[1]));
  }
  fseek(fptr, 0, SEEK_END);
  long ret_file_size = ftell(fptr);
  if (ret_file_size < 0) {
    throw std::logic_error("Unable to determine file size");
  }
  size_t file_size = (size_t) ret_file_size;
  fseek(fptr, 0, SEEK_SET);

  std::vector<char> read_buffer;
  read_buffer.resize(file_size, '\0');

  size_t ret = fread(read_buffer.data(),  1, file_size, fptr);
  printf("file_size=%lu ret=%lu\n", file_size, ret);
  if (ret != file_size) {
    throw std::logic_error("fread() failed\n");
  }
  fclose(fptr);
  const char *reads_beg = read_buffer.data();
  const char *reads_end = reads_beg + read_buffer.size();

#pragma omp parallel 
  {
    if (!omp_get_thread_num()) {
      tcount = omp_get_num_threads();
    }
  }
  printf("tcount=%lu\n", tcount);
  oneapi::tbb::global_control global_control(
      oneapi::tbb::global_control::max_allowed_parallelism, tcount);



  concurrent_kmer_counter_t kmer_counter(32, reads_beg, reads_end, 100, 
      std::numeric_limits<size_t>::max());

  count_unique_kmers_t counter;
  { 
    intel::esc::stop_clock_t 
          clock("[concurrent build elapsed time]", stdout);
    size_t rcount ;
    if (use_explicit_mode) {
      printf("mode=explicit\n");
      rcount = kmer_counter.concurrent_count_and_process_omp_explicit(counter);
    } else if (use_implicit_mode_v0) {
      printf("mode=implicit-v0\n");
      rcount =
        kmer_counter.concurrent_count_and_process_omp_implicit_v0(counter);
    } else {
      printf("mode=implicit-v1\n");
      rcount =
        kmer_counter.concurrent_count_and_process_omp_implicit_v1(counter);
    }
    printf("Reads = %lu , Unique Kmers=%lu freq=%lu\n", rcount, 
        (size_t)counter.count_, (size_t) counter.kmer_freq_);
  }
}
