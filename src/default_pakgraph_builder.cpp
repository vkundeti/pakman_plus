#include <mpi.h>

#include <incremental_pakgraph.hpp>
#include <mpi_rank_traits.hpp>

typedef uint64_t kmer_t;

struct unit_tests_32mer_traits_t {
 static inline size_t kmer_length() { return 32UL; }
}; // struct unit_tests_kmer_traits_t //

struct unit_tests_5mer_traits_t {
 static inline size_t kmer_length() { return 5UL; }
}; // struct unit_tests_kmer_traits_t //

typedef typename intel::esc::Distributed_Pak_Graph<kmer_t,
    unit_tests_32mer_traits_t, intel::esc::mpi_node_traits_t>
      distributed_pak_graph_t;


int main(int argc, char **argv) {

  //MPI_Init(&argc, &argv);

  int provided=-1, required=MPI_THREAD_MULTIPLE;

  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);

  if (provided != required) { 
    throw std::logic_error("[Thread Level Invalid]");
  }

  distributed_pak_graph_t graph_builder;
  if (argc < 2) {
    throw std::logic_error("USAGE: "
          "./default_pakgraph_builder {fasta file} [threshold=21]");
  }

  size_t threshold = 21;

  if ((argc >= 3) && (sscanf(argv[2], "%lu", &threshold)!=1)) {
    printf("[unable to scan threshold default threshold: 21]\n");
    threshold = 21;
  }

  threshold = std::max(1UL, threshold);
  
  graph_builder.default_build(argv[1], threshold-1UL);

  MPI_Finalize();
}
