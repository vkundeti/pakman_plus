#ifndef GRAPH_TRAITS_HPP
#define GRAPH_TRAITS_HPP
#include <mpi.h>
#include <iostream>

namespace intel {
namespace esc {


template<typename T>
struct graph_traits {
  typedef int edge_iterator_t;
  // vertex_handle_t should be unique for every vertex in the graph//
  typedef unsigned long vertex_t;
  typedef void vertex_record_traits;
  typedef bool weight_t;

  // weak ordering functor on vertices //
  // bool operator()(const vertex_t& a, const vertex_t& b) //
  typedef int vertex_ordering_t;

  // hash for vertices //
  // size_t operator()(const vertex_t& a) //
  typedef int vertex_hash_t;

  static vertex_t invalid_vertex();

  template<typename EdgeHandle>
  static vertex_t begin_vertex(const EdgeHandle&);

  template<typename EdgeHandle>
  static vertex_t end_vertex(const EdgeHandle&);

  template<typename T1>
  static bool mpi_data_type(const T1, MPI_Datatype& mpi_dtype);

  template<typename EdgeHandle>
  static weight_t weight(const EdgeHandle&);

}; // struct graph_traits //

// edge is (u,v) //
template<typename V, typename W>
struct Edge {
  typedef Edge edge_t;

  Edge(const V& u, const V& v, const W& w=W(0)) : u_(u), v_(v), w_() {} 

  bool operator<(const edge_t& o) const {
    return (v_ == o.v_) ? (u_ < o.u_) : (v_ < o.v_);
  }

  bool is_canonical() const {
    return (u_ < v_);
  }

  bool operator==(const edge_t& o) const {
    return (u_==o.u_) && (v_ == o.v_) && (w_ == o.w_);
  }

  void print() const {
    std::cout << "(u=" << u_ << " v=" << v_ << " w=" << w_ << ")" << std::endl;
  }

  V u_;
  V v_;
  W w_;
}; // struct edge_t //


} // namespace esc  //
} // namespace intel //
#endif
