#include <gtest/gtest.h>
#include <mpi.h>

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  ::testing::InitGoogleTest(&argc, argv);
  int result=0;
  result = RUN_ALL_TESTS();

  MPI_Finalize();
  return result;
}
