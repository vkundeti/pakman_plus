#include <cstdio>
#include <vector>
#include <string>
#include <errno.h>

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
    for (; beg!=end; ++beg) { ++count_; kmer_freq_ += beg->second; 
      //std::string kmer_str = kutils.get_kmer_string(beg->first);
      //printf("%s %lu\n", kmer_str.c_str(), beg->second);
    }
  }
  mutable std::atomic<size_t> count_;
  mutable std::atomic<size_t> kmer_freq_;
}; // struct count_unique_kmers_t //

int main(int argc, char **argv) {
  bool use_sort=false;

  if (argc < 2) {
    printf("Usage: ./concurrent-kmer-table-builder [--use_sort]");
    return -1;
  }

  if (argc >= 3) {
    use_sort = (strcmp(argv[2], "--use_sort") == 0);
  }

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

  concurrent_kmer_counter_t kmer_counter(32, reads_beg, reads_end, 100, 
      std::numeric_limits<size_t>::max());

  count_unique_kmers_t counter;
  { 
    intel::esc::stop_clock_t 
          clock("[default build elapsed time]", stdout);
    size_t rcount ;
    if (use_sort) {
      printf("mode=sort\n");
      rcount = kmer_counter.default_count_and_process(counter);
    } else {
      printf("mode=hash_map\n");
      rcount = kmer_counter.default_count_and_process(counter);
    }
    printf("Reads = %lu , Unique Kmers=%lu freq=%lu\n", rcount, 
        (size_t)counter.count_, (size_t)counter.kmer_freq_);
  }
}
