#ifndef RECORD_TRAITS_HPP
#define RECORD_TRAITS_HPP
#include <cassert>
#include <cstring>
#include <stdexcept>

namespace intel {
namespace esc {

template<typename T>
struct record_traits {
  typedef T record_t;
  static size_t record_length_in_bytes();
  // [start_byte_index, end_byte) //
  static size_t bytes_to_index(const record_t&,
        size_t start_byte, size_t end_byte);
  // [byte_begin, byte_end) //
  static void create_record(record_t&, void *byte_begin, void *byte_end);
  bool equivalent_record(const record_t&, const record_t&);
}; // struct record_traits //

// Some useful traits //
struct size_t_record_traits {
  typedef size_t record_t;

  static size_t record_length_in_bytes() { return sizeof(record_t); }
  static void create_record(record_t& record, void *byte_begin,
        void *byte_end) {
    void *res = memcpy(&record, byte_begin, record_length_in_bytes()); 
    ((void)res);
    assert(res == &record);
  }

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
}; // struct size_t_record_traits //

struct int_record_traits {
  typedef int record_t;

  static size_t record_length_in_bytes() { return sizeof(record_t); }
  static void create_record(record_t& record, void *byte_begin,
        void *byte_end) {
    void *res = memcpy(&record, byte_begin, record_length_in_bytes()); 
    ((void)res);
    assert(res == &record);
  }

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
}; // struct int_record_traits //

} // namespace esc
} // namespace intel
#endif
