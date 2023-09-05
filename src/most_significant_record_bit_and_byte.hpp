#ifndef MSB_RECORD_BIT_AND_BYTE_HPP
#define MSB_RECORD_BIT_AND_BYTE_HPP

#include <cstddef>
#include <mpi.h>
#include <record_traits.hpp>

namespace intel {
namespace esc {

template<typename RecordTraits>
class Most_Significant_Record_Bit_Byte {
  public:
    ////////////////////////////////////////////////////////////////////////////
    typedef RecordTraits record_traits;
    typedef typename record_traits::record_t record_t;

    struct msbb_t {
      msbb_t(int byte_idx=-1, unsigned char bit_idx=8)
        : byte_index_(byte_idx), bit_index_(bit_idx) {}

      bool is_valid() const {
        return (byte_index_ >= 0) && (bit_index_ < 8);
      }

      bool operator>(const msbb_t& o) const {
        return (byte_index_ == o.byte_index_) ? 
            (bit_index_ > o.bit_index_) : (byte_index_ > o.byte_index_);
      }

      bool operator==(const msbb_t& o) const {
        return (byte_index_ == o.byte_index_) && (bit_index_ == o.bit_index_);
      }

      void print(char c1='[', char c2=']', char c3='\n') const {
        printf("%cbyte_idx=%d bit_idx=%hhu%c%c",
              c1, byte_index_, bit_index_, c2, c3);
      }

      int byte_index_; // [0, record_lenght)
      unsigned char bit_index_; // [0,8)
    }; // struct msbb_t // 
    ////////////////////////////////////////////////////////////////////////////



    static size_t extract_bits_into_size_t(const record_t& record,
        const msbb_t& msbb, size_t bit_count=8) {
      size_t curr_byte_idx = msbb.byte_index_;
      size_t curr_bit_idx = msbb.bit_index_;
      size_t curr_bit_count = 0;
      size_t rlen = record_traits::record_length_in_bytes();
      size_t val =0;

      //first set of bits //
      {
        size_t shifted_bits = curr_bit_idx+1;
        size_t mask = size_t(1<<shifted_bits)-1;
        size_t rbyte = record_traits::bytes_to_index(record,
              curr_byte_idx, curr_byte_idx+1);
        rbyte &= mask;

        if (bit_count < shifted_bits) {
          rbyte >>= (shifted_bits - bit_count);
          shifted_bits = bit_count;
        }

        val <<= shifted_bits;
        val += rbyte;
        curr_bit_count += shifted_bits;
        // next byte and bit index to the left //
        curr_byte_idx--;
        curr_bit_idx=7UL;
      }

      // middle set of bytes //
      while ((curr_byte_idx < rlen) && (curr_bit_count < bit_count) 
            && ((bit_count - curr_bit_count) >= 8)) {
        val <<= 8;
        size_t rbyte = record_traits::bytes_to_index(record,
              curr_byte_idx, curr_byte_idx+1);
        val += rbyte;

        curr_bit_count +=8;
        curr_byte_idx--;
      }

      // last set of bits //
      if ((curr_byte_idx < rlen) && (curr_bit_count < bit_count)) {
        // this means (bit_count - curr_bit_count) < 8 //
        size_t delta = (bit_count - curr_bit_count);
        size_t rbyte = record_traits::bytes_to_index(record,
              curr_byte_idx, curr_byte_idx+1);
        rbyte >>= (8-delta);
        val <<= delta;
        val += rbyte;
      }
      return val;
    }

    static void apply_msbb_reduce(void *input_buffer, void *inout_buffer,
        int *len, MPI_Datatype *dtype ) {
      msbb_t *u = (msbb_t *)(input_buffer);
      msbb_t *v = (msbb_t *)(inout_buffer);
      msbb_t *output_array = (msbb_t *)(inout_buffer);
      for (int i=0; i<*len; i++) {
        //u[i].print('[',']','>'); v[i].print('[',']');
        output_array[i] = (u[i] > v[i]) ? u[i] : v[i];
      }
    }


    void create_msbb_type(MPI_Datatype &msbb_type) const {
      //[ (int, 0, 1)   , (unsigned char, offsetof(msbb_t, bit_index_, 1)) ] //
      MPI_Aint displacements[2] =
        { offsetof(msbb_t, byte_index_), offsetof(msbb_t, bit_index_) };
      int block_lengths[2] = {1, 1};
      MPI_Datatype native_types[2] = {MPI_INTEGER, MPI_UNSIGNED_CHAR};

      int ret = MPI_Type_create_struct(2, block_lengths, displacements,
            native_types, &msbb_type );
      if (ret != MPI_SUCCESS) {
        throw std::logic_error("[MPI_Type_create_struct]: failed");
      }
      ret = MPI_Type_commit(&msbb_type);
    }

    msbb_t extract_ms_bit_and_byte(const record_t& record) const {
      msbb_t ret_value;
      size_t rlen = record_traits::record_length_in_bytes();
      if (!rlen) { return ret_value; }

      // (byte:n-1) (byte:n-2) ..... (byte:0) //
      unsigned char val=0;
      for (size_t r=rlen-1; r<rlen; --r) {
        if ((val=(unsigned char)record_traits::bytes_to_index(record,r,r+1))) {
          ret_value.byte_index_ = (int) r; 
          break;
        }
      }

      if (val) {
        // find the highest bit in this byte //
        unsigned char b=(1<<7);
        unsigned char bit_index = 7;
        while (!(b&val)) {
          b >>= 1;
          bit_index--;
        }
        ret_value.bit_index_ = bit_index;
      }
      return ret_value;
    }

    // (byte_index (left-to-right), bit-index) //
    template<typename RecordPartitionIterator>
    msbb_t find(RecordPartitionIterator pbegin,
          RecordPartitionIterator pend) const {
      //TODO(vamsikku): static assert the iterator on record type //

      // largest byte and highest bit //


      msbb_t curr_max;
      for (; pbegin != pend; ++pbegin) {
        msbb_t cvalue = extract_ms_bit_and_byte(*pbegin);
        if (cvalue > curr_max) {
          curr_max = cvalue;
        }
      }

      bool is_commutable=true;
      MPI_Op operation;
      int ret = MPI_Op_create(
          Most_Significant_Record_Bit_Byte::apply_msbb_reduce,
          (int)(is_commutable), &operation);
      if (ret != MPI_SUCCESS) {
        throw std::logic_error("MPI_OP_create failed");
      }

      MPI_Datatype mpi_msbb_type; 
      create_msbb_type(mpi_msbb_type);

      msbb_t global_max;
      ret = MPI_Allreduce(&curr_max, &global_max, 1,
            mpi_msbb_type, operation, MPI_COMM_WORLD);
      if (ret != MPI_SUCCESS) {
        throw std::logic_error("MPI_Reduce failed");
      }
      return global_max;
    }


}; // class Most_Significant_Record_Bit_Byte //


} // namespace esc //
} // namespace intel //
#endif
