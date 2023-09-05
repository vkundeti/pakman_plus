#ifndef KMER_H
#define KMER_H

#include <inttypes.h>

#define KMER_LENGTH                    (WINDW_SIZE+1)
#define LMER_SIZE                      (pow(4, LMER_LENGTH))
#define MN_LENGTH                      (KMER_LENGTH-1)

#define MOD 2147483647
#define HL 31


#ifdef EXTEND_KMER
//static_assert(std::is_same_v<__uint128_t, unsigned __int128>);
//typedef unsigned __int128 uint128_t;
//typedef uint128_t kmer_t;
typedef __uint128_t kmer_t;
#else
typedef uint64_t kmer_t;
#endif



#endif
