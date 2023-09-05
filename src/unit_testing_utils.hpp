#ifndef UNIT_TESTING_UTILS_HPP
#define UNIT_TESTING_UTILS_HPP
#include <limits>

#include <record_traits.hpp>
#include <graph_traits.hpp>

struct unit_testing_record_traits {
  typedef size_t record_t;
  typedef double weight_t;

  static size_t record_length_in_bytes() { return sizeof(record_t); }
  static void create_record(record_t& record, void *byte_begin,
        void *byte_end) {
    void *res = memcpy(&record, byte_begin, record_length_in_bytes()); 
    ((void)res);
    assert(res == &record);
  }

  static void print_byte_string(const record_t& v, FILE *fptr=stdout) {
    size_t rlen = record_length_in_bytes(); 
    const unsigned char * rbytes = (const unsigned char *)(&v);
    fprintf(fptr, " ");
    for (size_t i=0; i<rlen; i++) {
      size_t val =  (size_t) rbytes[i];
      char bit=(1 << 7);
      fprintf(fptr, " [");
      for (size_t k=0; k<8; k++, bit >>= 1) {
        fprintf(fptr, "%d", bit&val ? 1 : 0);
      }
      fprintf(fptr, "] (%lu) ", val);
    }
    fprintf(fptr, " ");
  }

  //NOTE: sometimes when bucket_bytes is not a multiple of record length
  //we treat as if we had padded zeros on the right to make it a multiple
  static size_t bytes_to_index(const record_t& record,
        size_t start_byte, size_t end_byte) {
    size_t rlen = record_length_in_bytes();
    if ((start_byte >= end_byte) || (start_byte >= rlen)) {
      throw std::logic_error("[bytes_to_index]: invariant violoation "
            "start_byte="+std::to_string(start_byte)+
            " end_byte="+std::to_string(end_byte));
    }
    const unsigned char * rbytes = (const unsigned char *)(&record);

    size_t index = 0;
    for (--end_byte; end_byte > start_byte; --end_byte) {
      index = index*256 + 
        (size_t) ( (end_byte < rlen) ? rbytes[end_byte] : 0UL);
    }

    if (end_byte == start_byte) {
      index = index*256UL + (size_t) rbytes[start_byte];
    }
    return index;
  }

  static bool equivalent_record(const record_t& a, const record_t& b) {
    return (a == b);
  }

}; // struct unit_testing_record_traits //

struct unit_testing_graph_traits {
  typedef typename intel::esc::int_record_traits::record_t vertex_t;
  typedef intel::esc::int_record_traits vertex_record_traits;
  typedef double weight_t;
  struct edge_t { /* (u,v) */
    bool operator<(const edge_t& e) { return w_ < e.w_; }
    weight_t w_;
    vertex_t u_;
    vertex_t v_;
  }; // struct edge_t // 

  template<typename EdgeHandle>
  static vertex_t begin_vertex(const EdgeHandle&);

  template<typename EdgeHandle>
  static vertex_t end_vertex(const EdgeHandle&);

  template<typename EdgeHandle>
  static weight_t weight(const EdgeHandle&);
  
  inline static vertex_t invalid_vertex() { 
    return std::numeric_limits<vertex_t>::max();
  }

  template<typename T1>
  static bool mpi_data_type(const T1, MPI_Datatype& mpi_dtype);
}; // struct unit_testing_graph_traits //


template<>
inline typename unit_testing_graph_traits::vertex_t 
  unit_testing_graph_traits::begin_vertex(
    const typename unit_testing_graph_traits::edge_t& e) {
  return e.u_;
}

template<>
inline typename unit_testing_graph_traits::vertex_t 
  unit_testing_graph_traits::end_vertex(
    const typename unit_testing_graph_traits::edge_t& e) {
  return e.v_;
}

template<>
inline typename unit_testing_graph_traits::weight_t 
  unit_testing_graph_traits::weight(
      const typename unit_testing_graph_traits::edge_t& e) {
  return e.w_;
}

#endif
