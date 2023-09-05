#ifndef HRES_TIMERS_HPP
#define HRES_TIMERS_HPP

#include <chrono>
#include <string>

namespace intel {
namespace esc {

class stop_clock_t {
  public:
    ////////////////////////////////////////////////////////////////////////////
    typedef typename std::chrono::time_point<std::chrono::high_resolution_clock>
      time_sample_t;
    typedef typename std::chrono::high_resolution_clock hres_timer_t;
    ////////////////////////////////////////////////////////////////////////////

    stop_clock_t(const std::string& label, FILE *fptr=stdout,
                 bool dont_report_at_end=false)
      : label_(label), beg_(), end_(), fptr_(fptr), 
        dont_report_at_end_(dont_report_at_end) { beg_ = hres_timer_t::now(); }

    void reset_label(const std::string& new_label) {
      label_ = new_label;
    }

    void report_elapsed_time() {
      end_ = hres_timer_t::now();
      std::chrono::duration<double> time_span =
        std::chrono::duration_cast<std::chrono::duration<double>>(end_ - beg_);
      fprintf(fptr_, "[StopClock: %s]: %g sec\n\n",
            label_.c_str(), time_span.count());
      fflush(fptr_);
    }

    double get_elapsed_time_in_sec() const {
      time_sample_t end = hres_timer_t::now();
      std::chrono::duration<double> time_span =
        std::chrono::duration_cast<std::chrono::duration<double>>(end - beg_);
      return time_span.count();
    }

    const char* get_label() const { return label_.c_str(); }

    ~stop_clock_t() {
      if (!dont_report_at_end_) report_elapsed_time();
    }

  private:
    std::string label_;
    time_sample_t beg_;
    time_sample_t end_;
    FILE *fptr_;
    bool dont_report_at_end_;
}; // class stop_clock_t //

} // namespace esc //
} // namespace intel //
#endif
