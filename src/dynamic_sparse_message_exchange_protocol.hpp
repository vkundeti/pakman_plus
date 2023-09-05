#ifndef DYNAMIC_SPARSE_MESSAGE_EXCHANGE_PROTOCOL_HPP
#define DYNAMIC_SPARSE_MESSAGE_EXCHANGE_PROTOCOL_HPP
#include <mpi.h>
#include <list>
#include <unordered_map>

namespace intel {
namespace esc {

struct message_queue_traits {

template<typename MessageIterator>
static size_t destination_id(const MessageIterator& itr);

// Currently there is only 1 datatype the protocol works // 
static void mpi_data_type(MPI_Datatype&);

template<typename MessageIterator>
static void* mpi_data_buffer(const MessageIterator&);

template<typename MessageIterator>
static size_t element_count(const MessageIterator&);

template<typename MessageIterator>
static int channel(const MessageIterator&);

}; // struct message_traits //

struct async_message_processor_t {
  bool operator()(int source, int tag, const MPI_Datatype&,
        void *message, size_t elem_count);
}; // struct async_message_processor_t //


//TODO(vamsikku): parameterize by MPI_COMM_WORLD //
template<typename MessageQueueTraits>
class Dynamic_Sparse_Message_Exchange_Protocol {
  public:
    typedef MessageQueueTraits mq_traits;

    struct default_channel_selector_t {

      template<typename T>
      int operator()(T itr) {
        return mq_traits::channel(itr);
      }
    }; // default_channel_selector_t //


  template<typename MessageQueueIterator, typename AsyncMessageProcessor>
  static void async_send_receive_and_process(MessageQueueIterator mbegin, 
      MessageQueueIterator mend, AsyncMessageProcessor& message_processor,
      bool stateful_destination=false, bool skip_self_messages=false,
      size_t my_rank=0, bool cleanup_call_back=false, bool debug=false) {
    default_channel_selector_t channel_selector;

    async_send_receive_and_process_gen(mbegin, mend, message_processor,
        stateful_destination, skip_self_messages, my_rank, cleanup_call_back,
        channel_selector, debug);
  }


  // MessageProcessor 
  template<typename MessageQueueIterator, typename AsyncMessageProcessor,
    typename ChannelSelector=default_channel_selector_t>
  static void async_send_receive_and_process_gen(MessageQueueIterator mbegin, 
      MessageQueueIterator mend, AsyncMessageProcessor& message_processor,
      bool stateful_destination, bool skip_self_messages,
      size_t my_rank, bool cleanup_call_back, 
      ChannelSelector &channel_selector, bool debug=false) {
    typedef typename std::list<MPI_Request>::iterator send_request_iterator_t;
    std::list<MPI_Request> send_requests;
    std::unordered_map<MPI_Request *, void *> cleanup_buffers;
    MPI_Datatype data_type;
    mq_traits::mpi_data_type(data_type);
    std::vector<unsigned char> receive_buffer; 
    int dtype_size = -1;
    int channel = channel_selector(mbegin); //mq_traits::channel(mbegin);
    {
      
      int ret = MPI_Type_size(data_type, &dtype_size);
      if ((ret != MPI_SUCCESS) || (dtype_size <= 0)) {
        throw std::logic_error("[MPI_Type_size]: invalid datatype size");
      }

      MPI_Aint lb, extent;
      ret = MPI_Type_get_extent(data_type, &lb, &extent);
      if (ret != MPI_SUCCESS) {
        throw std::logic_error("[MPI_Type_get_extent] failed");
      }
      dtype_size = std::max(dtype_size, (int) extent);
    }

    // STEP-0: send all my messages using MPI_Issend //
    for (; mbegin!=mend; ++mbegin) {

      size_t destination = !stateful_destination ? 
          mq_traits::destination_id(mbegin) : 
          mq_traits::destination_id(mbegin, message_processor);

      if (skip_self_messages && (destination == my_rank) ) { continue; }

      if (debug) {
        printf("p-%lu sending message to p-%lu\n", my_rank, destination);
        fflush(stdout);
      }

      void *mbuffer = !cleanup_call_back ?
            mq_traits::mpi_data_buffer(mbegin)
          : mq_traits::mpi_data_buffer(mbegin, message_processor);
      size_t elem_count = mq_traits::element_count(mbegin);
      send_requests.push_back(MPI_REQUEST_NULL);
      int ret = MPI_Issend(mbuffer, (int)elem_count, data_type,
            (int) destination, channel, MPI_COMM_WORLD, &send_requests.back());
      if (ret != MPI_SUCCESS) {
        send_requests.pop_back();
        throw std::logic_error("[MPI_Issend] failed");
      }

      if (cleanup_call_back) {
        cleanup_buffers[ &send_requests.back() ] = mbuffer;
      }
    }

    //use sparse data exchange to receive data and distribute keys//
    MPI_Request ibarrier_request;
    bool entered_barrier_phase = false;
    do { // communication protocol: sparse data exchange  //
      //////////////////////////////////////////////////////////////////////
      // 1. Probe for a message and receive it if there is one            //
      //////////////////////////////////////////////////////////////////////
      MPI_Status status;
      int flag = 0;
      int ret = MPI_Iprobe(MPI_ANY_SOURCE, channel,
            MPI_COMM_WORLD, &flag, &status);
      if ((ret == MPI_SUCCESS) && flag) {
        //  1.1 If there is some message then receive it //
        int recv_count=0;

        //  1.2 Get the count of bytes //
        ret = MPI_Get_count(&status, data_type, &recv_count);
        if ((ret != MPI_SUCCESS) || (recv_count <0) ) {
          throw std::logic_error("MPI_Get_count: failed " 
              "or invalid recv_byte_count " +
                std::to_string(recv_count));
        }
        
        
        // 1.3 Receive //
        size_t buffer_size = size_t(dtype_size)*size_t(recv_count);
        receive_buffer.resize(buffer_size, (unsigned char)0);

        ret = MPI_Recv(receive_buffer.data(), recv_count, data_type,
            status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD,
              MPI_STATUS_IGNORE);

        if (ret != MPI_SUCCESS) {
          throw std::logic_error("MPI_Recv: failed ");
        }

        if (debug) {
         printf("p-%lu received %d edges from p-%d\n",
             my_rank, recv_count, status.MPI_SOURCE);
        }
        ///////////////////// PROCESS //////////////////////////////////////////
        bool pstatus = message_processor(status.MPI_SOURCE, status.MPI_TAG,
            data_type, receive_buffer.data(), (size_t) recv_count );
        if (!pstatus) {
          throw std::logic_error("Processing of message failed");
        }
        ////////////////////////////////////////////////////////////////////////
      } else if (ret != MPI_SUCCESS) {
        throw std::logic_error("MPI_Iprobe: failed");
      }
      //////////////////////////////////////////////////////////////////////

      if (!entered_barrier_phase) {
        //////////////////////////////////////////////////////////////////////
        // 2. Check if all my sent messages are received by some one.       //
        //////////////////////////////////////////////////////////////////////
        if (!send_requests.empty()) {
          std::vector<send_request_iterator_t> received_requests;
          for (send_request_iterator_t sitr=send_requests.begin();
                sitr!=send_requests.end(); ++sitr) {
            int ret = MPI_Test( &(*sitr), &flag, MPI_STATUS_IGNORE);

            if ((ret == MPI_SUCCESS) && flag) {

              if (!cleanup_buffers.empty()) {
                auto citr = cleanup_buffers.find( &(*sitr) );
                if (citr == cleanup_buffers.end()) {
                  throw std::logic_error("[Invalid cleanup buffer state]");
                }
                mq_traits::clean_buffer(message_processor, citr->second);
                cleanup_buffers.erase(citr);
              }


              received_requests.push_back(sitr);
            } else if (ret != MPI_SUCCESS) {
              throw std::logic_error("MPI_Test: failed ");
            }
          }
          for (auto itr=received_requests.begin();
                itr!=received_requests.end(); ++itr) {
            send_requests.erase(*itr);
          }
        }

        ////////////////////////////////////////////////////////////////////
        // 3. Go into barrier phase                                       //
        ////////////////////////////////////////////////////////////////////
        if (send_requests.empty()) {
          int ret = MPI_Ibarrier(MPI_COMM_WORLD, &ibarrier_request);
          if (ret != MPI_SUCCESS) {
            throw std::logic_error("MPI_Ibarrier: failed");
          }
          entered_barrier_phase = true;
#ifdef DEBUG_MPI
          printf("processor-%lu touched barrier\n", my_rank);
#endif

        } 
      }

      if (entered_barrier_phase) {
        // 3. Test for completion of the barrier call //
        int flag = 0;
        int ret = MPI_Test( &ibarrier_request, &flag, MPI_STATUS_IGNORE);
        if ((ret == MPI_SUCCESS) && (flag==1)) {
          break; // terminate communication //
        } else if (ret != MPI_SUCCESS) {
          throw std::logic_error("MPI_Test: failed");
        }
      }
    } while(1); // sparse data exchange protocol //

  }

}; // class Dynamic_Sparse_Message_Exchange_Protocol //



} // namespace esc //
} // namespace intel //

#endif
