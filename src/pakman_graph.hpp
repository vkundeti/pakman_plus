#ifndef PAKGRAPH_HPP
#define PAKGRAPH_HPP
#include <algorithm>
#include <cassert>
#include <cmath>
#include <generate_kmers.hpp>
#include <set>
#include <vector>

#include <kmer_traits.hpp>
#include <data_exchange_protocol.hpp>
#include <pakman_graph_dumper.hpp>

namespace intel {
namespace esc {

// Pakman_Graph: is a collection of sub-graphs across the processors //
template< typename KmerType, typename KmerTraits=kmer_traits_t<KmerType>,
          typename in_data_exchange_t=intel::esc::Data_Exchange_Protocol<>
        >
class Pakman_Subgraph {
  public:
  //////////////////////////////////////////////////////////////////////////////
  typedef KmerType kmer_t;
  typedef kmer_t k1mer_t; // 1 less symbol than a kmer //
  typedef kmer_utils_t<kmer_t> kutils_t;
  typedef KmerTraits kmer_traits;
  typedef in_data_exchange_t data_exchange_t;

  struct suffix_t {
    suffix_t(const std::string& str="", int count=0, int vcount=0,
          bool term=true)
      : str_(str), count_(count), visit_count_(vcount), has_term_flag_(term) {}

    bool is_terminal() const { return has_term_flag_; }
    void reset(int count=1, int visit_count=0, bool has_term=true) {
      count_ = count; visit_count_ = visit_count; has_term_flag_=has_term;
    }
    const std::string& str() const { return str_; }
    void set_string(const std::string& s) { str_=s; }
    void append_string(const std::string s) { str_ += s; }
    void set_term_flag(bool flag) { has_term_flag_ = flag; }
    int count_;
    int visit_count_;

    private: // enforce the invariant that //
    std::string str_;
    bool has_term_flag_;
  }; // struct suffix_t//
  typedef std::vector<suffix_t> suffixes_t;

  struct wire_t { // prefix to suffix connection //
    wire_t(const size_t suffix_idx, int count, int offset)
      : suffix_idx_(suffix_idx), count_(count), offset_(offset) {}

    const int& visit_count() const { return count_; }

    size_t suffix_idx_;
    int count_;
    int offset_;
  }; // struct wire_t //
  typedef std::vector<wire_t> wires_t;

  struct prefix_t : public suffix_t { // same as 
    prefix_t(const std::string& str="", int count=1, int vcount=0,
        bool term=true)
      :  suffix_t(str, count, vcount, term), wires_() {} 

    void clear_connections() { wires_.clear(); }
    void connect_to_suffix(size_t suffix_idx, int count, int offset) {
      wires_.push_back(wire_t(suffix_idx, count, offset) );
    }

    size_t wire_count() const { return wires_.size(); }
    wires_t wires_;
  }; // struct prefix_t //

  typedef std::vector<prefix_t> prefixes_t;

  struct visit_count_sorter_t {
    visit_count_sorter_t() : kutils_(kmer_traits::kmer_length()) {}

    inline bool operator()(const suffix_t& a, const suffix_t& b) const {
      return 
       (a.visit_count_ != b.visit_count_) ?  (a.visit_count_ > b.visit_count_) :
        (a.count_ != b.count_) ? (a.count_ > b.count_) :
          (a.is_terminal() != b.is_terminal()) ?
            (!a.is_terminal() > !b.is_terminal()) 
            : suffix_ordering(a.str(), b.str()) ; 
    }

    inline bool suffix_ordering(const std::string& suf_a,
          const std::string& suf_b) const {
      if (suf_a.length() != suf_b.length()) { 
        return suf_a.length() < suf_b.length();
      }

      for (size_t i=0; i<suf_a.length(); i++) {
        if (suf_a[i] != suf_b[i]) { 
          return kutils_.symbol_to_value(suf_a[i]) <
              kutils_.symbol_to_value(suf_b[i]);
        }
      }
      return false;
    }

    kutils_t kutils_;
  }; // struct visit_count_sorter_t//


  struct kmer_freq_t {
    kmer_freq_t(const kmer_t& kmer=kmer_t(), const size_t& freq=0UL) :
      kmer_(kmer), frequency_(freq) {}

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

    bool operator==(const kmer_freq_t& o) const { 
      return (kmer_ == o.kmer_) && (frequency_ == o.frequency_) ;
    }


    kmer_t kmer_;
    size_t frequency_;
  }; // struct kmer_freq_t //

  struct node_t {
    ////////////////////////////////////////////////////////////////////////////
    struct pred_iterator_t {
      pred_iterator_t(const node_t& nd, size_t pidx) : nd_(nd), pidx_(pidx) {
        if (!reached_end() && nd_.is_prefix_terminal(pidx_)) { ++(*this); }
      }
      bool operator==(const pred_iterator_t& o) const {
        return reached_end() && o.reached_end(); 
      }
      bool operator!=(const pred_iterator_t& o) const {
        return !(*this == o);
      }
      bool reached_end() const { return (pidx_ >= nd_.prefixes_.size()); }
      pred_iterator_t& operator++() {
        do { ++pidx_; } while(!reached_end() && nd_.is_prefix_terminal(pidx_));
        return *this;
      }
      kmer_t operator*() const { return nd_.pred_suffix_pair(pidx_).first; }
      const node_t& nd_;
      size_t pidx_;
    }; // pred_iterator_t //

    struct succ_iterator_t {
      succ_iterator_t(const node_t& nd, size_t sidx) : nd_(nd), sidx_(sidx) {
        if (!reached_end() && nd_.is_suffix_terminal(sidx_)) { ++(*this); }
      }
      bool operator==(const succ_iterator_t& o) const {
        return reached_end() && o.reached_end(); 
      }
      bool reached_end() const { return (sidx_ >= nd_.suffixes_.size()); }
      bool operator!=(const succ_iterator_t& o) const {
        return !(*this == o);
      }
      succ_iterator_t& operator++() {
        do { ++sidx_; } while(!reached_end() && nd_.is_suffix_terminal(sidx_));
        return *this;
      }
      kmer_t operator*() const { return nd_.prefix_succ_pair(sidx_).second; }
      const node_t& nd_;
      size_t sidx_;
    }; // succ_iterator_t //
    ////////////////////////////////////////////////////////////////////////////

    node_t(const k1mer_t& k1mer=kmer_t())
      :  k1mer_(k1mer), prefixes_(), suffixes_(), wires_(), 
        prefix_count_(), suffix_count_() {
        for (size_t i=0; i<4UL; i++) {
          prefix_count_[i] = 0UL;
          suffix_count_[i] = 0UL;
        }
        prefixes_.push_back(prefix_t("", 0UL, 0UL, true));
        suffixes_.push_back(suffix_t("", 0UL, 0UL, true));
    }

    void clear() { prefixes_.clear(); suffixes_.clear(); }

    void add_prefix_incremental(char c, size_t delta) {
      switch(c) {
        case 'A':
          prefix_count_[0]+=delta; return;

        case 'C':
          prefix_count_[1]+=delta; return;

        case 'T':
          prefix_count_[2]+=delta; return;

        default:
          prefix_count_[3]+=delta; 
      }
    }

    void add_suffix_incremental(char c, size_t delta) {
      switch(c) {
        case 'A':
          suffix_count_[0]+=delta; return;

        case 'C':
          suffix_count_[1]+=delta; return;

        case 'T':
          suffix_count_[2]+=delta; return;

        default:
          suffix_count_[3]+=delta; 
      }
    }

    bool is_balanced() const {
      int prefix_vc = 0, suffix_vc = 0;
      for (const prefix_t& p : prefixes_) {
        prefix_vc += p.visit_count_;
      }

      for (const suffix_t& s : suffixes_) {
        suffix_vc += s.visit_count_;
      }
      return (prefix_vc == suffix_vc);
    }

    // predecessor iterator //
    pred_iterator_t begin_pred() const { return pred_iterator_t(*this, 0UL); }
    pred_iterator_t end_pred() const {
      return pred_iterator_t(*this, prefixes_.size());
    }

    // successor iterator //
    succ_iterator_t begin_succ() const { return succ_iterator_t(*this, 0UL); }
    succ_iterator_t end_succ() const {
      return succ_iterator_t(*this, suffixes_.size());
    }

    bool is_prefix_terminal(size_t pidx) const { 
      return prefixes_[pidx].is_terminal(); 
    }
    bool is_suffix_terminal(size_t sidx) const {
      return suffixes_[sidx].is_terminal();
    }

    std::string concat_prefix_node(size_t prefix_idx) const {
      size_t k1len = kmer_traits::kmer_length()-1UL;
      kutils_t kutil(k1len);
      return kutil.append_kmer_to_string(prefixes_[prefix_idx].str(), k1mer_);
    }

    std::string k1mer_str() const {
      size_t k1len = kmer_traits::kmer_length()-1UL;
      kutils_t kutil(k1len);
      return kutil.get_kmer_string(k1mer_);
    }

    const wires_t& prefix_wires(size_t prefix_idx) const {
      return prefixes_[prefix_idx].wires_;
    }

    const std::string& prefix(size_t prefix_idx) const {
      return prefixes_[prefix_idx].str();
    }
    const size_t prefix_count(size_t prefix_idx) const {
      return prefixes_[prefix_idx].count_;
    }


    const std::string& suffix(size_t suffix_idx) const {
      return suffixes_[suffix_idx].str();
    }
    const size_t suffix_count(size_t suffix_idx) const {
      return suffixes_[suffix_idx].count_;
    }


    //TODO(vamsikku): replace k1len with a static variable state
    kmer_t predecessor_of_prefix(size_t prefix_idx) const {
      return pred_suffix_pair(prefix_idx).first;
    }

    kmer_t successor_of_suffix(size_t suffix_idx) const {
      return prefix_succ_pair(suffix_idx).second;
    }


    // Precondition: prefix_idx is not terminal//
    std::pair<kmer_t, std::string> pred_suffix_pair(size_t prefix_idx) const {
      size_t k1len = kmer_traits::kmer_length()-1UL;
      kutils_t kutil(k1len);
      std::string str =
          kutil.append_kmer_to_string(prefixes_[prefix_idx].str(), k1mer_);
      if (str.length() < k1len) {
        throw std::logic_error("[Prefix+k1mer invariant]: failed");
      }

      std::string prefix = str.substr(0, k1len);
      std::string suffix = str.substr(k1len);
      return std::make_pair(kutil.get_kmer(prefix), suffix);
    }

    // Precondition: suffix_idx is not terminal //
    std::pair<std::string, kmer_t> prefix_succ_pair(size_t suffix_idx) const {
      size_t k1len = kmer_traits::kmer_length()-1UL;
      kutils_t kutil(k1len);
      std::string str =
          kutil.prepend_kmer_to_string(suffixes_[suffix_idx].str(), k1mer_);
      if (str.length() < k1len) {
        throw std::logic_error("[Prefix+k1mer invariant]: failed");
      }

      std::string prefix = str.substr(0, (str.length()-k1len));
      std::string suffix = str.substr((str.length()-k1len));
      return std::make_pair(prefix, kutil.get_kmer(suffix));
    }

    bool does_node_meet_threshold(size_t threshold) const {
      for (size_t i=0; i<4; i++) {
        if ( (prefix_count_[i] >= threshold) || 
              (suffix_count_[i] >= threshold)){
          return true;
        }
      }
      return false;
    }

    void reset_terminals_and_clear_connections() {
      for (auto itr=prefixes_.begin(); itr!=prefixes_.end(); ++itr) {
        itr->clear_connections();
      }
    }

    void clear_wires() {
      for (prefix_t& prefix : prefixes_) {
        prefix.wires_.clear();
      }
    }

    void rebuild_prefix_suffix_connections() {
      visit_count_sorter_t sorter;

      reset_terminals_and_clear_connections();
      std::sort(prefixes_.begin(), prefixes_.end(), sorter);
      std::sort(suffixes_.begin(), suffixes_.end(), sorter);

      int total_visit_count = balance_suffix_prefix_visit_counts(), delta=0;
      int suf_used_visit_count = 0, pre_used_visit_count = 0;

      auto pitr = prefixes_.begin(), pitr_end = prefixes_.end();
      auto sitr = suffixes_.begin(), sitr_end = suffixes_.end();

      // Loop-Invariant1: suffix_used_visit_count <= sitr->visit_count_ // 
      // Loop-Invariant2: prefix_used_visit_count <= pitr->visit_count_ // 
      while ((pitr != pitr_end) && (sitr != sitr_end)) {
        delta = std::min(pitr->visit_count_ - pre_used_visit_count,
                         sitr->visit_count_ - suf_used_visit_count);
        pre_used_visit_count += delta; suf_used_visit_count += delta;
        // connection prefix to suffix //
        pitr->connect_to_suffix((sitr - suffixes_.begin()),
              delta, suf_used_visit_count);
        if (pre_used_visit_count == pitr->visit_count_) { 
          ++pitr; pre_used_visit_count = 0; //ensure Loop-Invariant1 holds //
        }
        if (suf_used_visit_count == sitr->visit_count_) {
          ++sitr; suf_used_visit_count = 0; //ensure Loop-Invariant2 holds //
        }
      }
    }

    //Precondition: suffixes and prefixes are sorted based on visit count //
    int balance_suffix_prefix_visit_counts() {
      int total_suf_count = 0UL;
      int total_pre_count = 0UL;

      for (const prefix_t& p : prefixes_) {
        total_pre_count += p.visit_count_;
      }
      for (const suffix_t& s : suffixes_) {
        total_suf_count += s.visit_count_;
      }

      if (total_suf_count < total_pre_count) {
        assert(!suffixes_.empty() && !(suffixes_.back().visit_count_));
        suffixes_.back().visit_count_=(total_pre_count - total_suf_count); 
        suffixes_.back().count_=1;
      } else if (total_suf_count > total_pre_count) {
        assert(!prefixes_.empty() && !(prefixes_.back().visit_count_));
        prefixes_.back().visit_count_=(total_suf_count - total_pre_count); 
        prefixes_.back().count_=1;
      }
      return std::max(total_suf_count, total_pre_count);
    }

    k1mer_t k1mer_;
    prefixes_t prefixes_;
    suffixes_t suffixes_;
    wires_t wires_;
    size_t prefix_count_[4];
    size_t suffix_count_[4];
  }; // struct node_t //
  ////////////////////////// message types [beg]////////////////////////////////
  typedef std::pair<size_t, kmer_freq_t> kmer_freq_mesg_t;
  typedef std::vector<kmer_freq_mesg_t> kmer_freq_messgs_t;
  struct mnode_compact_t {
    mnode_compact_t() : concat_string_(), append_terminal_(),
      is_k1mer_at_suffix_(), count_(), visit_count_(), flanking_length_() {}

    mnode_compact_t(const std::string& str, bool term, bool k1suf,
        int count, int visit_count, size_t flen) 
      : concat_string_(str), append_terminal_(term), is_k1mer_at_suffix_(k1suf),
        count_(count), visit_count_(visit_count), flanking_length_(flen) {}

    inline size_t encoded_length() const {
      // encode string ending with \0 //
      return (concat_string_.length()+1UL) 
             + (1UL) // pack append terminal and is_k1mer_at_suffix into byte
             + sizeof(int)
             + sizeof(int)
             + sizeof(size_t);
    }

    char* encode(char *byte_stream) const {
      if (concat_string_.length()) {
        memcpy(byte_stream, concat_string_.c_str(), concat_string_.length());
        byte_stream += concat_string_.length();
      }
      // add a null char //
      *byte_stream = '\0'; 
      byte_stream++;


      unsigned char val = ((unsigned char) append_terminal_) + 
          (2 * (unsigned char) is_k1mer_at_suffix_);

      *byte_stream = val;
      byte_stream++;

      memcpy(byte_stream, &count_, sizeof(int));
      byte_stream += sizeof(int);

      memcpy(byte_stream, &visit_count_, sizeof(int));
      byte_stream += sizeof(int);

      memcpy(byte_stream, &flanking_length_, sizeof(int));
      byte_stream += sizeof(size_t);
      return byte_stream;
    }

    char* decode(char *byte_stream) {
      concat_string_.clear();

      // decode concat_string //
      while(*byte_stream) {
        concat_string_.push_back(*byte_stream);
        ++byte_stream;
      }
      ++byte_stream; // move past null stream //

      append_terminal_ = (*byte_stream)&1;
      is_k1mer_at_suffix_ = (*byte_stream)&2;
      byte_stream++;

      memcpy(&count_, byte_stream, sizeof(int));
      byte_stream += sizeof(int);

      memcpy(&visit_count_, byte_stream, sizeof(int));
      byte_stream += sizeof(int);

      memcpy(&flanking_length_, byte_stream, sizeof(int));
      byte_stream += sizeof(size_t);
      return byte_stream;
    }

    bool operator==(const mnode_compact_t& o) const {
      return (concat_string_ == o.concat_string_) &&
        (append_terminal_ == o.append_terminal_) &&
        (is_k1mer_at_suffix_ == o.is_k1mer_at_suffix_) &&
        (count_ == o.count_) && (visit_count_ == o.visit_count_) &&
        (flanking_length_ == o.flanking_length_);
    }

    std::string concat_string_;
    bool append_terminal_;
    bool is_k1mer_at_suffix_;
    int count_;
    int visit_count_;
    size_t flanking_length_; // invariant: (k-1)+flanking_length_
  }; // struct mnode_compact_t//

  typedef std::pair<size_t, mnode_compact_t> mnode_compact_mesg_t;
  typedef std::vector<mnode_compact_mesg_t> mnode_compact_messgs_t;
  typedef std::vector<mnode_compact_t> mnode_compacts_t;
  ////////////////////////// message types [end]////////////////////////////////

  typedef std::unordered_map<kmer_t, node_t> local_nodes_t;
  typedef std::vector<kmer_freq_t> kmer_freqs_t;
  typedef std::vector< kmer_freq_t > remote_kmers_t; 
  typedef typename local_nodes_t::const_iterator const_local_node_iterator_t;
  typedef typename local_nodes_t::iterator local_node_iterator_t;
  typedef typename local_nodes_t::iterator node_handle_t;
  typedef std::vector<node_handle_t> node_handles_t;
  typedef std::vector< local_nodes_t > subgraph_t;
  typedef typename subgraph_t::const_iterator const_subgraph_iterator_t;
  typedef std::set<kmer_t> rewire_set_t;
  //////////////////////////////////////////////////////////////////////////////

  private:

  inline static size_t get_thread_count() {
    size_t n=0;
#pragma omp parallel
    {
      if (!omp_get_thread_num()) {
        n = omp_get_num_threads();
      }
    }
    return n;
  }

  public:

  Pakman_Subgraph(size_t klen=kmer_traits::kmer_length())
    : k1mer_partition_(), kmer_partition_(), klen_(klen), n_() 
  {
    n_ = get_thread_count();
    kmer_partition_.reset(klen_, n_);
    k1mer_partition_.reset(klen_-1UL, n_);
    local_nodes_.resize(n_);
    rewire_sets_.resize(n_);
  }

  Pakman_Subgraph(size_t klen, size_t n)
    : k1mer_partition_(), kmer_partition_(), klen_(klen), n_(n) 
  {
    kmer_partition_.reset(klen_, n_);
    k1mer_partition_.reset(klen_-1UL, n_);
    local_nodes_.resize(n_);
    rewire_sets_.resize(n_);
  }

  protected:


  void process_node_suffix(size_t tid, const kmer_t& k1mer,
      const std::string& suffix, size_t kmer_count, size_t coverage) {
    auto &nodes = local_nodes_[tid];
    auto itr = nodes.find(k1mer);

    if (itr == nodes.end()) {
      itr = nodes.insert(std::make_pair(k1mer, node_t(k1mer))).first;
    }
    node_t& node_info = itr->second;
    node_info.suffixes_.push_back(
        suffix_t(suffix, kmer_count, 
              (int) ::ceil((double)(kmer_count)/(double)coverage), false));
  }

  void process_node_prefix(size_t tid, const kmer_t& k1mer,
      const std::string& prefix, size_t kmer_count, size_t coverage) {
    auto &nodes = local_nodes_[tid];
    auto itr = nodes.find(k1mer);

    if (itr == nodes.end()) {
      itr = nodes.insert(std::make_pair(k1mer, node_t(k1mer))).first;
    }
    node_t& node_info = itr->second;
    node_info.prefixes_.push_back(
        prefix_t(prefix, kmer_count, 
              (int) ::ceil((double)(kmer_count)/(double)coverage), false));
  }

  public:

  template<typename partitioned_kmer_freq_table_t>
  void build_graph(const partitioned_kmer_freq_table_t& freq_table,
      size_t threshold, size_t cov=100) {
    size_t n = n_;

    // PAR-STEP-1.1[generate & xchg]: identify remote k-mers and process k-mers
    std::vector<kmer_freq_messgs_t> remote_kmer_xchg_in(n); // global-access //
    std::vector<kmer_freqs_t> remote_kmer_xchg_out(n); // local-access //
#pragma omp parallel 
    { // begin omp parallel //
      size_t tid = omp_get_thread_num();
      const auto &my_table = freq_table[tid];
      auto kbegin=my_table.begin();
      auto kend=my_table.end();
      auto &my_exchange_queue = remote_kmer_xchg_in[tid];
      for (auto itr=kbegin; itr!=kend; ++itr) {
        size_t kmer_frequency = itr->second;
        if (itr->second <= threshold) { continue; }

        kmer_t kmer = itr->first;
        size_t kmer_count = itr->second;
        kmer_freq_t kmer_freq(kmer, kmer_count);
        std::string kmer_str = kmer_partition_.get_kmer_string(kmer);
        std::string kmer_str_first = kmer_str.substr(0,1); 
        std::string kmer_str_last = kmer_str.substr(klen_-1UL, klen_);
        kmer_t prefix = kmer_partition_.prefix(kmer);
        kmer_t suffix = kmer_partition_.suffix(kmer);

        //prefix_tid < n_ , suffix_tid < n_ //
        size_t prefix_tid = k1mer_partition_(prefix);
        size_t suffix_tid = k1mer_partition_(suffix);

        if (prefix_tid == tid) {
          process_node_suffix(tid, prefix, kmer_str_last, kmer_count, cov);
        } else {
          my_exchange_queue.push_back(std::make_pair(prefix_tid, kmer_freq));
        }

        if (suffix_tid == tid) {
          process_node_prefix(tid, suffix, kmer_str_first, kmer_count, cov);
        } else {
          if (suffix_tid == prefix_tid) { continue; } //k-mer already sent //
          my_exchange_queue.push_back(std::make_pair(suffix_tid, kmer_freq));
        }
      } // for //
    } // end omp parallel //
    // PAR-STEP-1.2[exchange]: do the exchange using data_exchange collective // 
    data_exchange_t::exchange_data(
        remote_kmer_xchg_in, remote_kmer_xchg_out, kmer_freq_t());


    // PAR-STEP-2: process received kmers purely locally //
#pragma omp parallel
    {
      size_t tid = omp_get_thread_num();
      for (const kmer_freq_t& kmer_freq : remote_kmer_xchg_out[tid]) {
        const kmer_t& kmer = kmer_freq.kmer_;
        const size_t& kmer_count = kmer_freq.frequency_;
        std::string kmer_str =kmer_partition_.get_kmer_string(kmer);
        std::string kmer_str_first = kmer_str.substr(0,1); 
        std::string kmer_str_last = kmer_str.substr(klen_-1UL, klen_);
        kmer_t prefix = kmer_partition_.prefix(kmer);
        kmer_t suffix = kmer_partition_.suffix(kmer);

        size_t prefix_tid = k1mer_partition_(prefix);
        size_t suffix_tid = k1mer_partition_(suffix);
        if (prefix_tid == tid) {
          process_node_suffix(tid, prefix, kmer_str_last, kmer_count, cov);
        }

        if (suffix_tid == tid) {
          process_node_prefix(tid, suffix, kmer_str_first, kmer_count, cov);
        }
      }

      // finally build wire info //
      auto &my_nodes = local_nodes_[tid];
      for (auto itr=my_nodes.begin(); itr!=my_nodes.end(); ++itr) {
        (itr->second).rebuild_prefix_suffix_connections();
      }
    }
  } // build_graph() //

  template<typename ContigOutput>
  void compact_graph(ContigOutput contig_out) {
    size_t n = n_, k1len = klen_-1UL;

    //print_pakgraph_adj(*this);
    //print_pakgraph_suf_pre(*this);
    // PAR-STEP-1.1: find id_set (locally) //
    std::vector<node_handles_t> independent_sets(n);
    compute_local_independent_set(independent_sets);

    //dump_independent_set(*this, n, independent_sets);

    // PAR-STEP-2.1: [generate & xchg] messages //
    std::vector<mnode_compact_messgs_t> mn_compact_xchg_in(n);
    std::vector<mnode_compacts_t> mn_compact_xchg_out(n);
#pragma omp parallel
    {
      size_t tid = omp_get_thread_num();
      // clear rewire sets //
      rewire_sets_[tid].clear();
      auto &my_exchange_queue = mn_compact_xchg_in[tid];
      auto &my_nodes = local_nodes_[tid];
      for (node_handle_t& nitr : independent_sets[tid]) {
        const node_t& mnode = nitr->second;
        auto k1mer_str = mnode.k1mer_str();

        // foreach (prefix-suffix) connection//
        for (size_t pidx=0; pidx < mnode.prefixes_.size(); ++pidx) { 
          auto pred = mnode.predecessor_of_prefix(pidx);
          const auto& prefix_str = mnode.prefix(pidx);
          bool terminal_prefix = mnode.is_prefix_terminal(pidx);

          //TODO(vamsikku): create a new function //
          // mnode, pidx, contig_out, terminal_prefix, my_exchange_queue //
          for (const wire_t& w : mnode.prefix_wires(pidx)) {
            auto sidx = w.suffix_idx_;
            auto succ = mnode.successor_of_suffix(sidx);
            const auto& suffix_str = mnode.suffix(sidx);
            bool terminal_suffix = mnode.is_suffix_terminal(sidx);

            auto concat_str = prefix_str + k1mer_str + suffix_str;
            if (terminal_prefix && terminal_suffix) { // report contig //
              *(contig_out[tid]) = concat_str; 
              continue;
            }
           
            int count =
              std::min(mnode.prefix_count(pidx), mnode.suffix_count(sidx));
            if (!terminal_prefix && (pred != mnode.k1mer_)) {
              // has a valid predecessor //
              size_t pred_tid = k1mer_partition_(pred);
              bool append_term = terminal_suffix || (succ == mnode.k1mer_); 
              size_t flanking_len = suffix_str.length();
              mnode_compact_t mn_compact(concat_str, append_term, false, count,
                    w.visit_count(), flanking_len);

              if (pred_tid == tid) {
                compact_macro_node(tid, mn_compact);
              } else {
                my_exchange_queue.push_back(
                    std::make_pair(pred_tid, mn_compact));
              }
            }

            if (!terminal_suffix && (succ != mnode.k1mer_)) {
              // has a valid successor //
              size_t succ_tid = k1mer_partition_(succ);
              bool prepend_term = terminal_prefix || (pred == mnode.k1mer_); 
              size_t flanking_len = prefix_str.length();
              mnode_compact_t mn_compact(concat_str, prepend_term, true, count,
                    w.visit_count(), flanking_len) ;


              if (succ_tid == tid) {
                compact_macro_node(tid, mn_compact);
              } else {
                my_exchange_queue.push_back(
                    std::make_pair(succ_tid, mn_compact));
              }
            }
          } // foreach (prefix-suffix) connection //
        }

        // erase the node //
        my_nodes.erase(nitr);
      }
    }
    // PAR-STEP-2.2: [exchange]
    data_exchange_t::exchange_data(mn_compact_xchg_in, mn_compact_xchg_out,
        mnode_compact_t());

#pragma omp parallel
    {
      // PAR-STEP-3: update macro node suffixes (fully locally) //
      size_t tid = omp_get_thread_num();
      for (const mnode_compact_t& mn_compact : mn_compact_xchg_out[tid]) {
        compact_macro_node(tid, mn_compact);
      }

      // PAR-STEP=4: //
      for (const kmer_t& mnode_kmer : rewire_sets_[tid]) {
        rewire_mnode(tid, mnode_kmer);
      }
    }
  }

  void rewire_mnode(size_t tid, kmer_t mnode_kmer) {
    auto &my_nodes = local_nodes_[tid];
    auto nitr = my_nodes.find(mnode_kmer);
    if (nitr == my_nodes.end()) {
      throw std::logic_error("[rewire_mnode]: invariant failed missing kmer");
    }
    node_t& node = nitr->second;
    node.rebuild_prefix_suffix_connections();
  }

  void compact_macro_node(size_t tid, const mnode_compact_t& mn_compact) {
    size_t k1len = kmer_traits::kmer_length()-1UL;
    kutils_t kutil(k1len);
    const auto &concat = mn_compact.concat_string_;
    size_t flanking_len = mn_compact.flanking_length_;
    int count = mn_compact.count_, visit_count = mn_compact.visit_count_;
    bool append_terminal = mn_compact.append_terminal_;

    auto &my_nodes = local_nodes_[tid];
    auto &my_rewire_set = rewire_sets_[tid];
    // STEP-1: find the node //
    const std::string k1mer_str = mn_compact.is_k1mer_at_suffix_ ?
        concat.substr(concat.length()-k1len, k1len) : concat.substr(0, k1len);
    kmer_t search_k1mer = kutil.get_kmer(k1mer_str);

    auto nitr = my_nodes.find(search_k1mer);
    if (nitr == my_nodes.end()) {
      throw std::logic_error("[mnode_compact_t]: invariant failed missing kmer="
            + k1mer_str);
    }

    node_t& node = nitr->second;
    my_rewire_set.insert(search_k1mer);

    if (mn_compact.is_k1mer_at_suffix_) { 
      //        [ flanking len    ]                          [k1len]
      //concat= <.....new_str.....><...prefix_search_key....><k1mer> //
      auto prefix_search_key =
          concat.substr(flanking_len, concat.length()-(k1len+flanking_len));
      auto new_str = concat.substr(0, (concat.length()-k1len));
      find_and_update(node.prefixes_, prefix_search_key, new_str,
          count, visit_count, append_terminal);

    } else { 
      //        [k1mer]                           [ flanking len    ]                          
      //concat= <k1mer><....suffix_search_key....><...new_str.......> //
      auto suffix_search_key =
          concat.substr(k1len, concat.length()-(k1len+flanking_len));
      auto new_str = concat.substr(k1len, (concat.length()-k1len));
      find_and_update(node.suffixes_, suffix_search_key, new_str,
          count, visit_count, append_terminal);
    }
  }

  template<typename Container>
  void find_and_update(Container &input, 
      const std::string& search_key, const std::string& new_str,
      int count, int visit_count, bool append_terminal) {
    bool append_to_input = true;
    for (auto itr=input.begin(); itr!=input.end(); ++itr) {
      if (itr->str() == search_key) {
        if (!itr->is_terminal()) {
          // replace //
          itr->set_string(new_str);
          itr->reset(count, visit_count, append_terminal);
          append_to_input = false;
        } 
        break;
      }
    }

    if (append_to_input) {
      // append it //
      typedef typename Container::value_type fix_t;
      fix_t new_fix(new_str, count, visit_count, append_terminal);
      input.push_back(new_fix);
    }
  }

  template<typename Output>
  void compute_local_independent_set(Output& output) {
#pragma omp parallel 
    { // begin parallel //
      int tid = omp_get_thread_num();
      local_nodes_t &my_nodes = local_nodes_[tid];
      auto out_itr = std::back_inserter(output[tid]);
      for (auto itr=my_nodes.begin(); itr!=my_nodes.end(); ++itr) {
        kmer_t k1mer = itr->first;
        const node_t &node = itr->second;
        auto pred_itr = node.begin_pred();
        for (; pred_itr!=node.end_pred(); ++pred_itr) {
          if ( *pred_itr > k1mer) { break; }
        }
        if (pred_itr != node.end_pred()) { continue; }

        auto succ_itr = node.begin_succ();
        for (; succ_itr!=node.end_succ(); ++succ_itr) {
          if ( *succ_itr > k1mer) { break; }
        }
        if (succ_itr != node.end_succ()) { continue; }

        *out_itr = itr;
        ++out_itr;
      }
    } // end parallel //
  }

  const_local_node_iterator_t begin_nodes(size_t tid=0UL) const {
    if ( tid > local_nodes_.size()) {
      throw std::logic_error("[Thread count invariant failed]: ");
    }
    return (local_nodes_[tid]).begin();
  }

  const_local_node_iterator_t end_nodes(size_t tid=0UL) const {
    if ( tid > local_nodes_.size()) {
      throw std::logic_error("[Thread count invariant failed]: ");
    }
    return (local_nodes_[tid]).begin();
  }

  const_subgraph_iterator_t begin_subgraphs() const {
    return local_nodes_.begin();
  }

  const_subgraph_iterator_t end_subgraphs() const {
    return local_nodes_.end();
  }

  size_t kmer_length() const { return klen_; }

  size_t total_node_count() const {
    //TODO(vamsikku): use prefix sum
    size_t count=0UL;
    for (const auto & subgraph: local_nodes_) {
      count += subgraph.size();
    }
    return count;
  }
  
  protected:
  intel::esc::kmer_value_partition_t<kmer_t> k1mer_partition_;
  intel::esc::kmer_value_partition_t<kmer_t> kmer_partition_;
  size_t klen_;
  size_t n_;
  subgraph_t local_nodes_;
  std::vector<rewire_set_t> rewire_sets_;
}; // struct Pakman_Subgraph //


} // namesapce esc //
} // namespace intel //
#endif
