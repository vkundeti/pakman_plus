#ifndef MPI_RANK_TRAITS_HPP
#define MPI_RANK_TRAITS_HPP

#include <mpi.h>
#include <stdexcept>

namespace intel {
namespace esc {

struct mpi_node_traits_t {
  static inline int node_rank(MPI_Comm comm=MPI_COMM_WORLD) {
    int rank;
    int ret;
    ret = MPI_Comm_rank(comm, &rank);
    if (ret != MPI_SUCCESS) {
      throw std::logic_error("MPI_Comm_rank failed with code: "+
            std::to_string(ret));
    }
    return rank;
  }

  static inline int node_count(MPI_Comm comm=MPI_COMM_WORLD) {
    int node_count;
    int ret;
    ret = MPI_Comm_size(comm, &node_count);
    if (ret != MPI_SUCCESS) {
      throw std::logic_error("MPI_Comm_size failed with code: "+
            std::to_string(ret));
    }
    return node_count;
  }
}; // struct mpi_node_traits_t //

} // namespace esc //
} // namespace intel //
#endif
