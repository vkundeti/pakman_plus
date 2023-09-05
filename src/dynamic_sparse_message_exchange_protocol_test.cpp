#include <gtest/gtest.h>

#include <unit_testing_utils.hpp>
#include <dynamic_sparse_message_exchange_protocol.hpp>


struct unit_testing_message_queue_traits {

template<typename Iterator>
static size_t destination_id(const Iterator& itr) { return itr->first; }

template<typename Iterator, typename Processor>
static size_t destination_id(const Iterator& itr, Processor &) 
  { return itr->first; }

template<typename Iterator>
static size_t element_count(const Iterator& itr) {
  return (itr->second).size();
}

template<typename Iterator>
static int channel(const Iterator&) { return (int)0; }

template<typename Iterator>
static void* mpi_data_buffer(const Iterator& itr){
  return (itr->second).data();
}

template<typename Iterator, typename Processor>
static void* mpi_data_buffer(const Iterator& itr, Processor& processor){
  return (itr->second).data();
}

static void mpi_data_type(MPI_Datatype& data_type) {
  data_type = MPI_INTEGER;
}

template<typename Processor>
static void clean_buffer(Processor& processor, void *) {}

}; // struct unit_testing_message_queue_traits //

struct unit_testing_message_processor_t {

  bool operator()(int source, int tag, const MPI_Datatype& data_type,
      void *message, size_t elem_count) {
    int const *data = (const int *)message;

    for (size_t i=0; i<elem_count; i++) {
      recv_data_[source].push_back( data[i] );
    }

    return true;
  }

  std::map<size_t, std::vector<int> > recv_data_;
};

class Dynamic_Sparse_Message_Exchange_Protocol_Test : public ::testing::Test {
  
  
  virtual void SetUp() {
    MPI_Comm_size(MPI_COMM_WORLD, &N_);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank_);
  }

  virtual void TearDown() {
    MPI_Barrier(MPI_COMM_WORLD);
  }

  protected:
    int my_rank_;
    int N_;

    std::map<size_t, std::map<size_t, std::vector<int> > > send_data;

}; // class Dynamic_Sparse_Message_Exchange_Protocol_Test //

typedef typename intel::esc::Dynamic_Sparse_Message_Exchange_Protocol<
  unit_testing_message_queue_traits > protocol_t; 



TEST_F(Dynamic_Sparse_Message_Exchange_Protocol_Test, basic_test) {

  if (N_ < 3) {
    printf("Test needs atleast 3 processors\n");
    return;
  }

  std::map<size_t, std::map<size_t, std::vector<int> > > send_data;

  send_data[0] = { {1, {5,8,9}} , {2, {9,3,8} } }; // p-0 //
  send_data[1] = { {2, {1}} }; // p-1//
  send_data[2] = { {0, {11}} }; // p-2 //

  std::map<size_t, std::map<size_t, std::vector<int> > > expected_recv_data;
  expected_recv_data[0] = { {2, {11}} };
  expected_recv_data[1] = { {0, {5,8,9}} };
  expected_recv_data[2] = { {0, {9,3,8}}, {1, {1}}   };

  unit_testing_message_processor_t message_processor;



  {
    auto itr = send_data.find(my_rank_); 

    if (itr != send_data.end()) {
      protocol_t::async_send_receive_and_process(
        (itr->second).begin(), (itr->second).end(), message_processor);
    } else {
      std::map<size_t, std::vector<int> > empty;
      protocol_t::async_send_receive_and_process(empty.begin(), empty.end(),
            message_processor);
    }

    {
      auto itr = expected_recv_data.find(my_rank_); 
      if (itr != expected_recv_data.end()) {
        EXPECT_EQ(message_processor.recv_data_ , itr->second);
      } else {
        EXPECT_TRUE(message_processor.recv_data_.empty());
      }
    }
  }


}




