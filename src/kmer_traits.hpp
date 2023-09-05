#ifndef KMER_TRAITS_HPP
#define KMER_TRAITS_HPP

namespace intel {
namespace esc {

template<typename KmerType>
struct kmer_traits_t {
static inline size_t kmer_length(void) { return 32UL; }
}; // struct kmer_traits //

} // namespace esc //
} // namespace intel // 
#endif
