#ifndef GENERATE_KMERS_HPP
#define GENERATE_KMERS_HPP
#include <algorithm>
#include <cstring>
#include <unordered_map>
#include <stdexcept>
#include <iostream>

namespace intel {
namespace esc {
#if 0
template<typename kmer_t>
class Concurrent_Kmer_Frequency_Table {
  public:

    Concurrent_Kmer_Frequency_Table(size_t n) : n_(n) { init(); }

  private:
    size_t n_;
    std::vector<partial_table_t> frequency_tables_;
}; // class Concurrent_Kmer_Freqency_Table //
#endif

template<typename kmer_t>
struct partial_kmer_frequency_compare_t {
  typedef std::pair<kmer_t, size_t> key_t;
  bool operator()(const key_t& a, const key_t& b) const {
    return a.first < b.first;
  }
}; // struct partial_kmer_frequency_compare_t //

template<typename kmer_t>
struct kmer_frequency_t {
  kmer_frequency_t(const kmer_t& kmer=kmer_t(0), size_t frequency=0)
    : kmer_(kmer), frequency_(frequency) {}

  inline bool operator<(const kmer_frequency_t& o) const {
    return kmer_ < o.kmer_;
  }

  inline size_t encoded_length() const {
    return sizeof(kmer_t)+sizeof(size_t);
  }

  char* encode(char *byte_stream) const {
    memcpy(byte_stream, &kmer_, sizeof(kmer_t));
    byte_stream += sizeof(kmer_t);
    memcpy(byte_stream, &frequency_, sizeof(size_t));
    byte_stream += sizeof(size_t);
    return byte_stream;
  }

  char* decode(char *byte_stream) {
    memcpy(&kmer_, byte_stream, sizeof(kmer_t));
    byte_stream += sizeof(kmer_t);
    memcpy(&frequency_, byte_stream, sizeof(size_t));
    byte_stream += sizeof(size_t);
    return byte_stream;
  }

  bool operator==(const kmer_frequency_t& o) const { 
    return (kmer_ == o.kmer_) && (frequency_ == o.frequency_) ;
  }

  kmer_t kmer_;
  size_t frequency_;
}; // struct kmer_frequency_t //

template<typename kmer_t=uint64_t>
struct kmer_utils_t {
  kmer_utils_t(size_t klen) : klen_(klen) {}

  char get_last_symbol_of_kmer(const kmer_t& kmer) const {
    kmer_t mask = (kmer_t)(3);
    kmer_t v = mask & kmer;
    switch (v) {
      case 0:
        return 'A';
      case 1:
        return 'C';
      case 2:
        return 'T';
      default:
        return 'G';
    }
  }

  inline std::string append_kmer_to_string(const std::string& in_bp_str,
      const kmer_t& kmer, size_t in_len=std::numeric_limits<size_t>::max()) {
    std::string out_bp_str = in_bp_str + get_kmer_string(kmer, in_len);
    return out_bp_str;
  }

  inline std::string prepend_kmer_to_string(const std::string& in_bp_str,
      const kmer_t& kmer, size_t in_len=std::numeric_limits<size_t>::max()) {
    std::string out_bp_str = get_kmer_string(kmer, in_len) + in_bp_str;
    return out_bp_str;
  }

  kmer_t symbol_to_value(char c) const {
    switch (c) {
      case 'A':
        return kmer_t(0);
      case 'C':
        return kmer_t(1);
      case 'T':
        return kmer_t(2);
      case 'G':
        return kmer_t(3);
      default:
        throw std::logic_error("Invalid symbol");
    }
  }

  std::string get_kmer_string(const kmer_t& kmer,
      size_t in_len=std::numeric_limits<size_t>::max()) const {
    kmer_t shift_kmer = kmer;
    std::string kmer_string;
  
    size_t len = std::min(klen_, in_len);
    char c = get_last_symbol_of_kmer(shift_kmer);
    kmer_string.push_back(c);
    for (size_t i=1; i<len; ++i) {
      shift_kmer >>= 2;
      c = get_last_symbol_of_kmer(shift_kmer);
      kmer_string.push_back(c);
    }
    std::reverse(kmer_string.begin(), kmer_string.end());
    return kmer_string;
  }

  kmer_t get_kmer(const std::string& str) const {
    if (str.length() > klen_) {
      throw std::logic_error("[Kvalue Invariant Failed]: str.len()="
            +std::to_string(str.length())+" <=" + std::to_string(klen_));
    }
    kmer_t kmer = kmer_t(0);

    for (size_t i=0; i<str.length(); i++) {
      kmer <<= 2UL;
      kmer += symbol_to_value(str[i]);
    }
    return kmer;
  }

  size_t klen_;
}; // struct kmer_utils_t //

template<typename in_kmer_t=uint64_t>
struct kmer_value_partition_t {
  typedef in_kmer_t kmer_t;
  kmer_value_partition_t(const size_t& klen=32UL, const size_t& n=1UL,
      size_t min_size=4UL)
    : klen_(klen), n_(n), chunk_size_() , min_size_(min_size)
  {
    if ((n < 1) || !klen_ || (klen_ > 32)) {
      throw std::logic_error("Invalid partition state ");
    }

    if (klen_ == 32) {
      max_value_ = (kmer_t)(-1);
    } else {
      max_value_ = (kmer_t)((kmer_t)(1)<<(2*klen_)) - (kmer_t)(1);
    }
    chunk_size_ = std::max(min_size, max_value_/n);
    printf("[Partition Info]: kval=%lu nval=%lu max_val_=%lu chunk_size_=%lu\n",
        klen, n_, max_value_, chunk_size_);

  }

  void reset(size_t klen, size_t n) {
    klen_ = klen;
    n_ = n;
    if ((n < 1) || !klen_ || (klen_ > 32)) {
      throw std::logic_error("Invalid partition state ");
    }
    if (klen_ == 32) {
      max_value_ = (kmer_t)(-1);
    } else {
      max_value_ = (kmer_t)((kmer_t)(1)<<(2*klen_)) - (kmer_t)(1);
    }
    chunk_size_ = std::max(min_size_, max_value_/n_);
    printf("[Partition Info]: kval=%lu nval=%lu max_val_=%lu chunk_size_=%lu\n",
        klen, n_, max_value_, chunk_size_);
  }

  size_t operator()(const kmer_t& kmer) const {
    return std::min(n_-1, kmer/chunk_size_); 
  }

  kmer_t min_value_of_partition(size_t thread_id) const {
    kmer_t part_min = ((kmer_t)thread_id)*((kmer_t)chunk_size_);
    return (part_min > max_value_) ? kmer_t(max_value_) : part_min;
  }

  kmer_t max_value_of_partition(size_t thread_id) const {
    kmer_t part_min = ((kmer_t)thread_id)*((kmer_t)chunk_size_);
    if (part_min > max_value_) { return kmer_t(0); }
    if (thread_id == n_-1) { return max_value_; } // avoid overflow in kmer_t//

    kmer_t part_max =
      (kmer_t)((kmer_t)thread_id+kmer_t(1))*(kmer_t)chunk_size_ - (kmer_t)(1);
    return std::min(part_max, max_value_);
  }

  bool has_valid_partition(size_t thread_id) const {
    return min_value_of_partition(thread_id) <= 
        max_value_of_partition(thread_id);
  }

  kmer_t max_kmer_val() const {
    return max_value_; 
  }

  kmer_t prefix(const kmer_t& kmer) const {
    kmer_t ret = kmer;
    return (ret >> 2);
  }

  kmer_t suffix(const kmer_t& kmer) const {
    kmer_t mask = (kmer_t(1) << (2UL*(klen_-1UL)))-kmer_t(1) ;
    return (kmer & mask);
  }

  kmer_t get_kmer(const std::string& str) const {
    if (str.length() > klen_) {
      throw std::logic_error("[Kvalue Invariant Failed]: str.len()="
            +std::to_string(str.length())+" <=" + std::to_string(klen_));
    }
    kmer_t kmer = kmer_t(0);

    for (size_t i=0; i<str.length(); i++) {
      kmer <<= 2UL;
      kmer += symbol_to_value(str[i]);
    }
    return kmer;
  }

  kmer_t symbol_to_value(char c) const {
    switch (c) {
      case 'A':
        return kmer_t(0);
      case 'C':
        return kmer_t(1);
      case 'T':
        return kmer_t(2);
      case 'G':
        return kmer_t(3);
      default:
        throw std::logic_error("Invalid symbol");
    }
  }

  char get_last_symbol_of_kmer(const kmer_t& kmer) const {
    kmer_t mask = (kmer_t)(3);
    kmer_t v = mask & kmer;
    switch (v) {
      case 0:
        return 'A';
      case 1:
        return 'C';
      case 2:
        return 'T';
      default:
        return 'G';
    }
  }

  std::string get_kmer_string(const kmer_t& kmer,
      size_t in_len=std::numeric_limits<size_t>::max()) const {
    kmer_t shift_kmer = kmer;
    std::string kmer_string;
  
    size_t len = std::min(klen_, in_len);
    char c = get_last_symbol_of_kmer(shift_kmer);
    kmer_string.push_back(c);
    for (size_t i=1; i<len; ++i) {
      shift_kmer >>= 2;
      c = get_last_symbol_of_kmer(shift_kmer);
      kmer_string.push_back(c);
    }
    std::reverse(kmer_string.begin(), kmer_string.end());
    return kmer_string;
  }


  void dump() const {
    printf("\n\nchunk_size=%lu\n", chunk_size_);
    for (kmer_t k=0; k<=max_kmer_val(); k++) {
      std::cout << "kmer=" << k << " partition=" << (*this)(k) << std::endl;
    }
  }

  size_t max_value_;
  size_t klen_;
  size_t n_;
  size_t chunk_size_;
  size_t min_size_;
}; // struct kmer_value_partition_t //

template<typename KmerNativeType=uint64_t, typename KmerNativeTypeExt=uint64_t>
class Generate_Kmers_From_Read {
  public:
    ////////////////////////////////////////////////////////////////////////////
    typedef KmerNativeType kmer_t;
    typedef std::unordered_map<char, kmer_t> symbol_table_t;
    ////////////////////////////////////////////////////////////////////////////
    Generate_Kmers_From_Read(size_t klen=1UL,
        const char*rbeg=NULL, const char *rend=NULL) : symbol_table_(), mask_(),
      klen_(klen), rbeg_(), rend_(), rlen_(0UL), curr_idx_(), kmer_() {
        if (!klen || (kmer_bits_needed(klen) >= bits_for_kmer_type()) ) {
          throw std::logic_error("[Klen invariant violation]: klen >=1 and "
              "kmer should fit in native type");
        }
        init_once();
        init(rbeg_, rend_);
    }

    size_t bits_for_kmer_type() const {
      size_t bit_count = sizeof(kmer_t)*8UL;
      return bit_count;
    }


    size_t kmer_bits_needed(size_t klen) const {
      size_t len = 1, bit_count = 1;
      while (len < klen) {
        ++bit_count;
        len <<= 1;
      }
      return (2*bit_count);
    }

    char get_last_symbol_of_kmer(const kmer_t& kmer) const {
      kmer_t mask = (kmer_t)(3);
      kmer_t v = mask & kmer;
      switch (v) {
        case 0:
          return 'A';
        case 1:
          return 'C';
        case 2:
          return 'T';
        default:
          return 'G';
      }
    }

    char get_first_symbol_of_kmer(const kmer_t& kmer) const {
      kmer_t shifted_kmer = (kmer) >> (2*(klen_-1)); 
      return get_last_symbol_of_kmer(shifted_kmer);
    }

    std::string get_kmer_string(const kmer_t& kmer) const {
      kmer_t shift_kmer = kmer;
      std::string kmer_string;
     
      char c = get_last_symbol_of_kmer(shift_kmer);
      kmer_string.push_back(c);
      for (size_t i=1; i<klen_; ++i) {
        shift_kmer >>= 2;
        c = get_last_symbol_of_kmer(shift_kmer);
        kmer_string.push_back(c);
      }
      std::reverse(kmer_string.begin(), kmer_string.end());
      return kmer_string;
    }

    bool is_valid() const {
      return (rend_ > rbeg_) && ((rend_ - rbeg_) >= klen_);
    }

    void reset(const char *beg, const char *end) {
      init(beg, end);
    }

    // precondition: reached_end() is false //
    Generate_Kmers_From_Read& operator++() {
      ++curr_idx_;
      if (reached_end()) { return *this; }
      kmer_ <<= 2;
      char c = rbeg_[(klen_-1UL)+curr_idx_];
      check_symbol_validity(c);
      kmer_ += symbol_table_[c];
      kmer_ = mask_ & kmer_;
      return *this;
    }

    const kmer_t& operator*() const { return kmer_; }

    bool operator==(const Generate_Kmers_From_Read& other) const {
      return (!is_valid() || reached_end()) && 
             (!other.is_valid() || other.reached_end());

    }

    bool operator!=(const Generate_Kmers_From_Read& other) const {
      return !(*this == other);
    }

  private:

    // precondition: is_valid() must be true //
    bool reached_end() const {
      return (curr_idx_ > rlen_) || 
        (rlen_ - curr_idx_) < klen_;
    }

    void init_once() {
      // 2^{2*klen_} - 1 // 111111111111 (sequence of 2*klen_ ones)
      mask_ = (klen_==32) ? ((kmer_t)-1) : 
        (kmer_t) (kmer_t(1) << (2*klen_)) - (kmer_t) 1;
      symbol_table_['A'] = kmer_t(0);
      symbol_table_['C'] = kmer_t(1);
      symbol_table_['T'] = kmer_t(2);
      symbol_table_['G'] = kmer_t(3);
    }

    void init(const char *rbeg, const char *rend) {
      if (klen_ < 1) { 
        throw std::logic_error("K-value invariant violation");
      }

      rbeg_ = rbeg; rend_ = rend;
      if (!is_valid()) { return; }
      rlen_ = (rend_ - rbeg_);
      curr_idx_ = 0UL;
      init_first_kmer();
    }

    // precondition: is_valid() must be true //
    void init_first_kmer() {
      char c = *rbeg_;
      check_symbol_validity(c);
      kmer_ = (kmer_t)(symbol_table_[c]);
      for (size_t i=1; i<klen_; i++) {
        kmer_ <<= 2;
        c = rbeg_[i];
        check_symbol_validity(c);
        kmer_ += symbol_table_[c];
      }
    }

    void check_symbol_validity(char c) const {
      if ( !((c=='A') || (c=='T') || (c=='G') || (c=='C')) ) {
        fprintf(stderr, "[WARNING] Invalid symbol %c in the read: ");
        for (const char *ptr=rbeg_; ptr!=rend_; ++ptr) {
          fprintf(stderr, "%c", *ptr);
        }
        fprintf(stderr, "\n");
      }
    }

    symbol_table_t symbol_table_;
    kmer_t mask_;
    size_t klen_;
    const char *rbeg_;
    const char *rend_;
    size_t rlen_;
    size_t curr_idx_;
    kmer_t kmer_;
}; // class Generate_Kmers_From_Read //



} // namespace esc //
} // namespace intel //
#endif
