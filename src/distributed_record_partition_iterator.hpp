#ifndef DISTRIBUTED_RECORD_PARTITION_ITERATOR_HPP 
#define DISTRIBUTED_RECORD_PARTITION_ITERATOR_HPP

#include <cstdio>
#include <list>
#include <stdexcept>
#include <string>
#include <unistd.h>
#include <vector>
#include <cerrno>
#include <cstring>


namespace intel {
namespace esc {

// Distributed_Record_Partition_Iterator: iterate over distributed records
// in parallel over a file system which can seek in O(1) time //
class Distributed_Record_Partition_Iterator {
  public:
    Distributed_Record_Partition_Iterator(const std::string& delim="",
          FILE *fptr=NULL, size_t n=1, size_t my_rank=0, size_t min_size=4096,
          bool iterator_mode=true)
      : delim_(delim), fptr_(fptr), n_(n), my_rank_(my_rank),
        begin_offset_(), end_offset_(), min_size_(min_size), process_buffer_(),
        circular_buffer_(), close_file_ptr_(false)
    {
      if (fptr_) {
        if (!delim.length()) {
          throw std::logic_error("Delimiter cannot be empty!");
        }
        init_offsets();
        if (iterator_mode && !stream_reached_end()) {
          init_circular_buffer();
          has_next_record();
        }
      }

    }

  private:

#if 0
    Distributed_Record_Partition_Iterator (
        const Distributed_Record_Partition_Iterator& o) {
      if (!o.fptr_) {
        if (fptr_ && close_file_ptr_) {
          fclose(fptr_);
          close_file_ptr_ = false;
        }
        fptr_ = NULL;
        begin_offset_ = o.begin_offset_;
        end_offset_ = o.end_offset_;
        return;
      }
      // reopen the file stream //
      int fd = fileno(o.fptr_);
      printf("fd returned = %d \n", fd);
      fptr_ = fdopen(fd, "r");
      printf("strerror=%s\n", strerror(errno));
      if (!fptr_) {
        throw std::logic_error("Failed to get file descriptor"
              " in copy constructor");
      }
      long seek_offset = ftell(o.fptr_);
      if (seek_offset < 0) {
        throw std::logic_error("Unable to determine the seek offset");
      }

      int ret = fseek(fptr_, seek_offset, SEEK_SET);
      if (ret) {
        throw std::logic_error("Unable to determine the seek offset");
      }

      close_file_ptr_ = true;
      n_ = o.n_;
      my_rank_ = o.my_rank_;
      begin_offset_ = o.begin_offset_;
      end_offset_ = o.end_offset_;
      min_size_ = o.min_size_;
      process_buffer_ = o.process_buffer_;
      circular_buffer_ = o.circular_buffer_;
    }

    Distributed_Record_Partition_Iterator& operator=(
        const Distributed_Record_Partition_Iterator& o) {

      if (!o.fptr_) {
        if (fptr_ && close_file_ptr_) {
          fclose(fptr_);
          close_file_ptr_ = false;
        }
        fptr_ = NULL;
        begin_offset_ = o.begin_offset_;
        end_offset_ = o.end_offset_;
        return *this;
      }
      //NOTE: user owns 

      FILE *orig_fptr_ = fptr_;
      int fd = fileno(o.fptr_);
      fptr_ = fdopen(fd, "r");
      if (!fptr_) {
        throw std::logic_error("Failed to get file descriptor"
              " in copy constructor");
      }

      if (close_file_ptr_) {
        fclose(orig_fptr_);
        close_file_ptr_ = false;
      }
      long seek_offset = ftell(o.fptr_ );
      if (seek_offset < 0) {
        throw std::logic_error("Unable to determine the seek offset");
      }

      int ret = fseek(fptr_, seek_offset, SEEK_SET);
      if (ret) {
        throw std::logic_error("Unable to determine the seek offset");
      }

      close_file_ptr_ = true;
      n_ = o.n_;
      my_rank_ = o.my_rank_;
      begin_offset_ = o.begin_offset_;
      end_offset_ = o.end_offset_;
      min_size_ = o.min_size_;
      process_buffer_ = o.process_buffer_;
      circular_buffer_ = o.circular_buffer_;
    }
#endif

  public:

    inline bool has_next_record() {
      if (circular_buffer_.empty()) {
        throw std::logic_error("Invariant violation: empty circular buffer");
      }
      // append data until next match or EOF or empty record //
      do {
        while (!is_match(delim_.begin(), delim_.end(),
                circular_buffer_.begin(), circular_buffer_.end())) {
          process_buffer_.push_back(circular_buffer_.front());

          if(!read_next_symbol_into_circular_buffer()) { break; }
        }
      } while (process_buffer_.empty() &&  // make progress on empty-records //
               read_next_symbol_into_circular_buffer());
      return !process_buffer_.empty();
    }

    const std::string& operator*() const { return process_buffer_; }

    inline bool next() {
      process_buffer_.clear();
      return read_next_symbol_into_circular_buffer();
    }
    inline void operator++() {
      if (next()) {
        has_next_record();
      }
    }

    // Only invalid iterators are equivalent //
    inline bool operator==(
        const Distributed_Record_Partition_Iterator& o) const {
      return reached_end() && o.reached_end();
    }

    inline bool operator!=(
        const Distributed_Record_Partition_Iterator& o) const {
      return !(*this == o);
    }

    std::pair<size_t, size_t> get_offsets() const {
      return std::make_pair(begin_offset_, end_offset_);
    }

    size_t curr_offset() const { return begin_offset_; }

    void print_status() const {
      printf("processor-%lu [%lu,%lu] %0.5lf \n", my_rank_,
          begin_offset_, end_offset_, double(begin_offset_/end_offset_));
    }

  private:

    bool reached_end() const {
      return process_buffer_.empty() && stream_reached_end();
    }


    bool read_next_symbol_into_circular_buffer() {
      if (!circular_buffer_.empty()) {
        circular_buffer_.pop_front();
      }
      if (!stream_reached_end()) {
        int c = fgetc(fptr_);
        if (c == EOF) { return false; }
        circular_buffer_.push_back((char) c);
        ++begin_offset_;
      }
      return !circular_buffer_.empty();
    }

    bool stream_reached_end() const { return (begin_offset_ >= end_offset_); }

    void init_circular_buffer() {
      for (size_t i=0; (i<delim_.size()) && (begin_offset_ < end_offset_);
          i++, begin_offset_++) {
        int byte = fgetc(fptr_);
        if (byte == EOF) { return; }
        circular_buffer_.push_back((char)byte);
      }
    }

    void init_offsets() {
      if (!fptr_) { return; }
      // determine the file-size in bytes //

      fseek(fptr_, 0, SEEK_END);
      long ret_file_size = ftell(fptr_);
      if (ret_file_size < 0) {
        throw std::logic_error("Unable to determine file size");
      }
      size_t file_size = (size_t) ret_file_size;
      size_t chunk_size = std::max(min_size_,
            file_size/std::min(file_size, n_));

      begin_offset_ = chunk_size*my_rank_;
      end_offset_ = file_size;
      if (begin_offset_> file_size) {
        begin_offset_ = file_size;
        end_offset_ = begin_offset_;
        return;
      }

      end_offset_ = std::min( file_size, begin_offset_ + chunk_size);

#if 0
      printf("file_size=%lu chunk_size = %lu b=%lu e=%lu\n",
          file_size, chunk_size, begin_offset_, end_offset_);
#endif


      if (my_rank_) {
        begin_offset_ = adjust_offset_until_delimiter(begin_offset_);
      }

      if (my_rank_ +1 < n_) {
        end_offset_ = adjust_offset_until_delimiter(end_offset_);
      } else {
        end_offset_ = file_size;
      }

      if (!stream_reached_end() && 
          fseek(fptr_, (long) begin_offset_, SEEK_SET) < 0) {
        throw std::logic_error("Bad seek");
      }
    }

    size_t adjust_offset_until_delimiter(size_t current_offset) {
      if (fseek(fptr_, (long) current_offset, SEEK_SET) < 0) {
        throw std::logic_error("Bad seek");
      }

      std::list<char> circular_buffer;
      int byte;
      for (size_t i=0; (i<delim_.size() && ((byte=fgetc(fptr_))!=EOF)); i++) {
          circular_buffer.push_back((char)byte);
          ++current_offset;
      }

      while (!is_match(delim_.begin(), delim_.end(),
            circular_buffer.begin(), circular_buffer.end()) &&
            ((byte=fgetc(fptr_))!=EOF)) {
        circular_buffer.pop_front();
        circular_buffer.push_back((char)(byte));
        ++current_offset;
      }
      return current_offset;
    }

    template<typename IteratorA, typename IteratorB>
    bool is_match(IteratorA abegin, IteratorA aend,
          IteratorB bbegin, IteratorB bend) const {
      while(!(abegin == aend) && !(bbegin == bend) && (*abegin == *bbegin)) {
        ++abegin; ++bbegin;
      }
      return (abegin == aend) && (bbegin == bend);
    }


    ////////////////////////////////////////////////////////////////////////////
    std::string delim_;
    FILE *fptr_;
    size_t n_; // number of processors //
    size_t my_rank_; //my rank //
    size_t begin_offset_;
    size_t end_offset_;
    size_t min_size_;
    std::string process_buffer_;
    std::list<char> circular_buffer_;
    bool close_file_ptr_;
    ////////////////////////////////////////////////////////////////////////////
}; // class Distributed_Record_Partition_Iterator//


} // namespace esc //
} //namespace intel //
#endif
