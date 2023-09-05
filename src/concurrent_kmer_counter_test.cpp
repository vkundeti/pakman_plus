#include <gtest/gtest.h>
#include <concurrent_kmer_counter.hpp>
#include <generate_kmers.hpp>

static const char *READ_TEST1=
">CACVBX020000001.1-NC_000913.3520\n"
"TACCATTCACCCGCTGGCGGGGCAGGGGGTAAATCTCGGCTTTATGGATGCTGCAGAGCTGATTGCCGAACTGAAACGGTTGCATCGTCAGGGGAAAGAC\n"
">CACVBX020000001.1-NC_000913.3519\n"
"TGGTTTGGCCGTCGCTTAGGCCAGGATGCCGCACCAGAGAAACTGGGCGTTGACCGCGCACTGGCATTCAGAGAAATCAATCAGCAACTTTATCAACTAC\n"
">CACVBX020000001.1-NC_000913.3518\n"
"GTTTGAAACTGGCCGACACGCTTCCTGGCGTTAAGCCGCAACTTATCCGCCAGGCAATGGGATTAAACGATTTGCCTGAATGGCTGCGTTAAAAATTTCT\n"
">CACVBX020000001.1-NC_000913.3517\n"
"TCGATTTCTTCCGGCTTTTTCACCACGCTGGCGGTGAGGTTTTCGTTACCGGTTTCTGTAATCACAATGTCGTCTTAAATACGAATGCCGATACCGCGAT\n"
">CACVBX020000001.1-NC_000913.3516\n"
"CCTTAGCCACTGGTTAGGACTGGATGTCCATGACGTGGGTGTTTATGGTCAGGATCGCGCGCGCATTCTGGAACCGGGCATGGTACTGACCGTAGAGCCA\n"
">CACVBX020000001.1-NC_000913.3515\n"
"GGTTAAGTCACGGGGAGCTGCCGGTACATTTGATTGAAGCGACTCCGCCAGAGTCACATGCTCATCCGGGCTTTGATGGACGAGCGATAGCGCTGGCGGC\n"
">CACVBX020000001.1-NC_000913.3514\n"
"CAACGTTGCTACTTCCGTTGCGCATGAAGGGCGCGCTTTTGAACGCTTTACGCAACATGGCCCGCTGGCGATGTTGCCGATGTCTGACGGACGCTGTTCG\n"
">CACVBX020000001.1-NC_000913.3513\n"
"GCGCTATCGCTCGTCCATCAAAGCCCGGATGAGCATGTGACTCTGGCGCAGTCGCTTCAATCAAATGTACCGGCAGCGCCCCGTGACTGAACCGGGAAAT\n"
">CACVBX020000001.1-NC_000913.3512\n"
"TTCATTGGGCAGCGCAATTTCATAGCCCGCTTCACCGGTATAACCAGTGGTGGCAATAAACAGATCGCCCGCCTGCACGCCAAAGAACGGTTTCATCCCT\n"
">CACVBX020000001.1-NC_000913.3511\n"
"CGCTGAAACTCTTGCCGGGATATCTCACTCATAACACTCTCCTTACGTTTTTTGTTTTTACTGTAGAGTCGGTTTTTGTACTTCTGGCGCGGTCGGTTGC\n"
">CACVBX020000001.1-NC_000913.3510\n"
"GATTGCCAGACGCCGATGCGCGCCAATTGCTGACAGGTACCCGCCGCCAGCGCTATCGCTCGTCCATCAAAGCCCGGATGAGCATGTGACTCTGGCGCAG\n"
">CACVBX020000001.1-NC_000913.3509\n"
"TTAGCGGCCTGGTAAAACTCGGCATCCTGAAAGGTGATGTTGATGAACTGATCGCTCAGAACGCCCATCGTCCTGTCTTTATGCATGGCCTTAGCCACTG\n"
">CACVBX020000001.1-NC_000913.3508\n"
"GGTCAACGCCCAGTTTCTCTGGCGCGGCATCCTGGCCTAAGCGACGGCCAAACCAGATCTCCGCCGTCAGGTCGCGAACGCGGTTGAACAGAACGCTGTG\n"
">CACVBX020000001.1-NC_000913.3507\n"
"CAAGAGTTTCAGCGTCGCCGTCAGGCCCTGGGGGAGCAAATGCAACCCGGCAGCGCCGCGCTGATTTTTGCTGCACCAGAAGTAACACGCAGCGCCGACA\n"
">CACVBX020000001.1-NC_000913.3506\n"
"TTTCATTGGGCAGCGCAATTTCATAGCCCGCTTCACCGGTATAACCAGTGGTGGCAATAAACAGATCGCCCGCCTGCACGCCAAAGAACGGTTTCATCCC\n"
">CACVBX020000001.1-NC_000913.3505\n"
"TCTGGTGCAGCAAAAATCAGCGCGGCGCTGCCGGGTTGCATTTGCTCCACCAGGGCCTGACGGCGACGCTGAAACTCTTGCCGGGATATCTCACTCATAA\n"
">CACVBX020000001.1-NC_000913.3504\n"
"ACGCCGCACTGTGCTTGCGGCTACGCTCATAGCGACGCAGATAAATGTACTGCCCGATGTCTTTCCCCTGACGATGCAACCGTTTCAGTTCGGCAATCAG\n"
">CACVBX020000001.1-NC_000913.3503\n"
"ACCTGGCGTAACGCTGCATTGCCCTGATCGCGTGGCTAACGTTGCCCGTACTCAGAGTCACGTTGAAGTGACGCTGGAGAGTGGCGAGACGCTGACGGGC\n"
">CACVBX020000001.1-NC_000913.3502\n"
"ACGCCAAGACGGGTGAGTAATTTTTCGCTGGCGGCATTGATAGCCGAAACGCGCAGTTGTGGTGGTGCATTCGCCGCAAGAGGTTCCTGTACGCGCTGCT\n"
"\n";

class Concurrent_Kmer_Counter_Test:  // READ, KLEN, THREAD_COUNT //  
    public testing::TestWithParam <std::tuple<std::string, size_t, size_t> > {
  protected:
    ////////////////////////////////////////////////////////////////////////////
    struct concurrent_kmer_processor_t {
      typedef intel::esc::kmer_utils_t<size_t> kmer_utils_t;

      concurrent_kmer_processor_t(size_t klen) : klen_(klen), kmer_count_(0),
        total_freq_count_(0) {}

      template<typename Iterator>
      void operator() ( Iterator begin, Iterator end) const {
        kmer_utils_t kmer_utils(klen_);
        for (; begin!=end; ++begin) {
          //std::string kmer_string = kmer_utils.get_kmer_string(begin->first);
          //printf("%s %lu\n", kmer_string.c_str(), begin->second);
          ASSERT_TRUE(begin->second >= 1UL);
          ++kmer_count_;
          total_freq_count_ += (size_t) begin->second;
        }
      }

      size_t klen_;
      mutable std::atomic<size_t> kmer_count_;
      mutable std::atomic<size_t> total_freq_count_;
    }; // struct concurrent_kmer_processor_t //

    ////////////////////////////////////////////////////////////////////////////


    typedef intel::esc::Concurrent_Kmer_Partial_Counter<size_t>
        concurrent_kmer_counter_t;
    virtual void SetUp() { }

}; // class Generate_Kmers_From_Read_Test //

TEST_P(Concurrent_Kmer_Counter_Test, concurrent_test) {
  std::string input = std::get<0>(GetParam());
  size_t klen = std::get<1>(GetParam());
 
  const char *reads_beg = input.c_str();
  const char *reads_end = reads_beg + input.length();
  concurrent_kmer_counter_t kmer_counter(klen, reads_beg, reads_end, 100,
      std::numeric_limits<size_t>::max());

  concurrent_kmer_processor_t kmer_processor(klen);

  size_t rcount = 
      kmer_counter.concurrent_count_and_process_omp_v3(kmer_processor);
  printf("reads = %lu\n", rcount);
  EXPECT_EQ( (size_t)kmer_processor.total_freq_count_, rcount*((100-klen+1)));
}


TEST_P(Concurrent_Kmer_Counter_Test, default_test) {
  std::string input = std::get<0>(GetParam());
  size_t klen = std::get<1>(GetParam());
 
  const char *reads_beg = input.c_str();
  const char *reads_end = reads_beg + input.length();
  concurrent_kmer_counter_t kmer_counter(klen, reads_beg, reads_end, 100,
      std::numeric_limits<size_t>::max());

  concurrent_kmer_processor_t kmer_processor(klen);

  size_t rcount = kmer_counter.default_count_and_process(kmer_processor);
  printf("reads = %lu\n", rcount);
  EXPECT_EQ( (size_t)kmer_processor.total_freq_count_, rcount*((100-klen+1)));
}

//TODO: fix readlen hardcoded from 100 //
INSTANTIATE_TEST_SUITE_P(Concurrent_Kmer_Counter_Test_Suite,
     Concurrent_Kmer_Counter_Test,
    ::testing::Values( 
      
      std::make_tuple(
    ">CACVBX020000001.1-NC_000913.3520\n"
    "TACCATTCACCCGCTGGCGGGGCAGGGGGTAAATCTCGGCTTTATGGATGCTGCAGAGCTGATTGCCGAACTGAAACGGTTGCATCGTCAGGGGAAAGAC\n"
    ">CACVBX020000001.1-NC_000913.3519\n"
    "TGGTTTGGCCGTCGCTTAGGCCAGGATGCCGCACCAGAGAAACTGGGCGTTGACCGCGCACTGGCATTCAGAGAAATCAATCAGCAACTTTATCAACTAC\n" , 32, 4)
                       
      )
);
