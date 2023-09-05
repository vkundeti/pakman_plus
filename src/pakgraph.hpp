#ifndef BUILD_PAKGRAPH_HPP
#define BUILD_PAKGRAPH_HPP

#include <type_traits>
#include "distributed_integer_sort.hpp"
#include "record_traits.hpp"
#include "kmer.h"

namespace intel {
namespace esc {

struct suffix_prefix_node_t {
  std::static_assert(sizeof(kmer_t)==sizeof(size_t),
      "[Kmer Size Invariant]: kmer type and size_t type mismatch");
}; // struct suffix_prefix_node_t //

class Pakman_Graph {

  public:

  private:

}; // class Pakman_Graph //

} // namespace esc //
} // namespace intel //
#endif
