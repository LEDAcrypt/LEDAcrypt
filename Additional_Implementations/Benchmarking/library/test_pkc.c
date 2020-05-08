#include "test_pkc.h"

#include "qc_ldpc_parameters.h"
#include "gf2x_limbs.h"
#include "gf2x_arith_mod_xPplusOne.h"

#include "mceliece.h"
#include "mceliece_keygen.h"
#include "mceliece_cca2_encrypt.h"
#include "mceliece_cca2_decrypt.h"
#include "api.h"
#include "rng.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>    // C99 sqrtl(...) for sample variance computation
#include <assert.h>

#include <time.h> // struct timespec; clock_gettime(...); CLOCK_REALTIME
#include "timing_and_stat.h"

/*----------------------------------------------------------------------------*/

void print_PKC_mceliece_cca2_parameters(int machine_readable)
{
if(machine_readable == 0) {
   fprintf(stderr,
           "\n -------------------------------------------------------------");
   fprintf(stderr,
           "\n LEDAcrypt-PKC  IND-CCA2 Kobara-Imai gamma (KI-g) Scheme       ");
   fprintf(stderr,
           "\n                                                             ");
   fprintf(stderr,"\n   N0:............................................%5d      ",
           N0);
   fprintf(stderr,"\n   P:.............................................%5d(b)  ",
           P);
   fprintf(stderr,
           "\n   K/N:...............................%5d/%5d = %2.1lf    ", K, N,
           K/((double)N));
   fprintf(stderr,
           "\n                                                             ");
   fprintf(stderr,"\n   H circ. block weight  V:.......................%5d      ",
           V);
   fprintf(stderr,"\n   numb. of errors. T:............................%5d      ",
           NUM_ERRORS_T);
   fprintf(stderr,
           "\n                                                             ");
   fprintf(stderr,"\n   True-RNG material byte length:.................%5d (B)  ",
           TRNG_BYTE_LENGTH);
   fprintf(stderr,"\n   Hash digest byte length:.......................%5d (B)  ",
           HASH_BYTE_LENGTH);
   fprintf(stderr,"\n   Min. Ciphertext Overhead.......................%5lu (B)  ",
           N0*NUM_DIGITS_GF2X_ELEMENT*DIGIT_SIZE_B-MAX_BYTES_IN_IWORD);
    fprintf(stderr,"\n   Max. Ciphertext Overhead.......................%5d (B)  ",
           N0*NUM_DIGITS_GF2X_ELEMENT*DIGIT_SIZE_B);   fprintf(stderr,
           "\n                                                             ");
   fprintf(stderr,
           "\n   private key size:....................%5lu(B) = %2.1lf(KiB)  ",
           (unsigned long) sizeof(privateKeyMcEliece_t),
           sizeof(privateKeyMcEliece_t) / ((double) 1024));
   fprintf(stderr,
           "\n   public  key size:....................%5lu(B) = %2.1lf(KiB) ",
           (unsigned long) sizeof(publicKeyMcEliece_t),
           (sizeof(publicKeyMcEliece_t)) / ((double) 1024));
   fprintf(stderr,
           "\n                                                             ");
   fprintf(stderr,
           "\n -------------------------------------------------------------\n");
} else {
   fprintf(stderr,
           "\n SIZE: %d & %s & %d & %5lu & %5lu & %5lu & %5lu \\\\ \n",
           CATEGORY,
         (DFR_SL_LEVEL==0) ? "$2^{-64}$ " :
                             ( (CATEGORY==1) ? "$2^{-128}$" :
                                               ((CATEGORY==3) ?"$2^{-192}$" : "$2^{-256}$") ),
           N0,
           (unsigned long) CRYPTO_SECRETKEYBYTES,
           (unsigned long) CRYPTO_PUBLICKEYBYTES, /* pubkey */
           (unsigned long) (N0*NUM_DIGITS_GF2X_ELEMENT*DIGIT_SIZE_B-MAX_BYTES_IN_IWORD),/* min overhead */
           (unsigned long) (CRYPTO_BYTES)); /*max overhead*/
}
} // end print_PKC_mceliece_cca2_parameters

/*----------------------------------------------------------------------------*/

static
double compute_time_interval(struct timespec *start, struct timespec *end)
{

   double total_time=0.0;
   long int delta_ns, delta_s;
   if (start->tv_nsec==end->tv_nsec) {
      if (start->tv_sec==end->tv_sec) {
         total_time=(double) (end->tv_sec - start->tv_sec);
      } else {
         return total_time;
      }
   }
   if (start->tv_nsec < end->tv_nsec) {
      delta_ns=end->tv_nsec - start->tv_nsec;
      delta_s=end->tv_sec - start->tv_sec;
      return ((double) delta_s)+ ((double) delta_ns)/1000000000.0;
   }

   if (start->tv_nsec >= end->tv_nsec) {
      delta_ns= (1000000000-start->tv_nsec) + end->tv_nsec;
      delta_s = end->tv_sec - (start->tv_sec + 1);
      return ((double) delta_s)+ ((double) delta_ns)/1000000000.0;
   }
   return -1.0f;
} // end compute_time_interval


/*----------------------------------------------------------------------------*/

void print_results_machine_readable(char tag[],
                                    welford_t keygen_stats,
                                    welford_t encrypt_stats,
                                    welford_t decrypt_stats){
 fprintf(stderr,
         "\n %s: %d & %s & %d ",
         tag,
         CATEGORY,
         (DFR_SL_LEVEL==0) ? "$2^{-64}$ " :
                             ( (CATEGORY==1) ? "$2^{-128}$" :
                                               ((CATEGORY==3) ?"$2^{-192}$" : "$2^{-256}$") ),
         N0
        );
 fprintf(stderr,
         "& %4.3Lf ($\\pm$ %.3Lf) & %4.3Lf ($\\pm$ %.3Lf) & %4.3Lf ($\\pm$ %.3Lf) & %4.3Lf \\\\ \n",
         welford_mean(keygen_stats),
         welford_stddev(keygen_stats),
         welford_mean(encrypt_stats),
         welford_stddev(encrypt_stats),
         welford_mean(decrypt_stats),
         welford_stddev(decrypt_stats),
         (welford_mean(encrypt_stats) + welford_mean(decrypt_stats)));
}

/*----------------------------------------------------------------------------*/

void print_results_human_readable(welford_t keygen_stats,
                                    welford_t encrypt_stats,
                                    welford_t decrypt_stats,
                                    long int NumTests,
                                    long int decodeOk,
                                    long int memcmpOk){
      fprintf(stderr,"\n\nPerformance tests.");
      fprintf(stderr, "\n\n%12sNumber of tests: %10lu", " ",
              (long unsigned int)NumTests);
      fprintf(stderr, "\nAverage key generation time: %14.3Lf (+,- %.3Lf) millisec",
              welford_mean(keygen_stats),
              welford_stddev(keygen_stats));

      fprintf(stderr,"\nNumber of correct decodings: %10ld \
                     \n  ---  ---      decryptions: %10ld",
              decodeOk, memcmpOk);

      fprintf(stderr,"\n%4sAverage encryption time: %14.3Lf (+,- %.3Lf) millisec",
              " ",
              welford_mean(encrypt_stats),
              welford_stddev(encrypt_stats)
             );
      if (memcmpOk > 1)
         fprintf(stderr,
                 "\n%4sAverage decryption time: %14.3Lf (+,- %.3Lf) millisec ... with successful decoding, and ptx and decrypted ctx match",
                 " ",
                welford_mean(decrypt_stats),
                welford_stddev(decrypt_stats)
                );
      else
         fprintf(stderr,"\nNo successful enc/dec match !!!");
}

/*----------------------------------------------------------------------------*/
#if (defined HIGH_PERFORMANCE_X86_64)
#define METRIC_ALPHA_FACTOR 1000
void print_cycles_metrics(char tag[],
                                    welford_t keygen_stats,
                                    welford_t encrypt_stats,
                                    welford_t decrypt_stats){
 fprintf(stderr,
         "\n %s: %d & %s & %d ",
         tag,
         CATEGORY,
         (DFR_SL_LEVEL==0) ? "$2^{-64}$ " :
                             ( (CATEGORY==1) ? "$2^{-128}$" :
                                               ((CATEGORY==3) ?"$2^{-192}$" : "$2^{-256}$") ),
         N0
        );
 fprintf(stderr,
         "& %4.1Lf ($\\pm$ %.1Lf) & %4.1Lf ($\\pm$ %.1Lf) & %4.1Lf ($\\pm$ %.1Lf) & %4.1Lf & %4.1Lf\\\\ \n",
         welford_mean(keygen_stats)/1E+3,
         welford_stddev(keygen_stats)/1E+3,
         welford_mean(encrypt_stats)/1E+3,
         welford_stddev(encrypt_stats)/1E+3,
         welford_mean(decrypt_stats)/1E+3,
         welford_stddev(decrypt_stats)/1E+3,
         (welford_mean(encrypt_stats) + welford_mean(decrypt_stats))/1E+3,
         (welford_mean(encrypt_stats)+
          welford_mean(decrypt_stats)+
          (CRYPTO_BYTES)*METRIC_ALPHA_FACTOR)/1E+3);
}

void collect_and_print_cycles_metrics(long unsigned int NumTests){
   char seed[40];
   initialize_pseudo_random_generator_seed(0, seed);

   unsigned char *pk= calloc(1, CRYPTO_PUBLICKEYBYTES);
   unsigned char *sk= calloc(1, CRYPTO_SECRETKEYBYTES);
   unsigned long long mlen=1024;
   unsigned char *plaintext = calloc(1, mlen);
   memset(plaintext,0xFF,mlen);
   unsigned char *ct= calloc(1, mlen+CRYPTO_BYTES);
   unsigned char *plaintextDec = calloc(1, CRYPTO_BYTES);
   unsigned long long clen=0;
   unsigned long long mlenDec=0;
   uint64_t beginning, end;

   welford_t keygen_time_stat;
   welford_init(&keygen_time_stat);

   for (unsigned trials = 0; trials < NumTests; trials++) {
      beginning = x86_64_rtdsc();
      crypto_encrypt_keypair(pk, sk);
      end = x86_64_rtdsc();
      welford_update(&keygen_time_stat, (end-beginning));
   }

   welford_t encrypt_time_stat;
   welford_t decrypt_time_stat;
   welford_init(&encrypt_time_stat);
   welford_init(&decrypt_time_stat);
   long int ok, decodeOk = 0, memcmpOk = 0;

   crypto_encrypt_keypair(pk, sk);
   for (unsigned trials = 0; trials < NumTests; trials++) {
      for (unsigned i = 0 ; i < mlen; i++){
            plaintext[i] = rand();
      }
      beginning = x86_64_rtdsc();
      crypto_encrypt( ct,&clen,plaintext,mlen,pk );
      end = x86_64_rtdsc();
      welford_update(&encrypt_time_stat, (end-beginning));

      beginning = x86_64_rtdsc();
      ok = crypto_encrypt_open(plaintextDec,&mlenDec,ct,clen,sk);
      end = x86_64_rtdsc();

      decodeOk += (ok == 0); /*NIST API has 0 = success, 1= failure */
      ok        = ( memcmp(plaintext,plaintextDec,mlenDec) == 0 );
      memcmpOk += ok;
      if (ok){
          welford_update(&decrypt_time_stat, (end-beginning));
      }
   } // end for trials
   print_cycles_metrics("CYCLES",
                        keygen_time_stat,
                        encrypt_time_stat,
                        decrypt_time_stat);
}
#endif


/*----------------------------------------------------------------------------*/

void test_PKC_mceliece_cca2_code(int isSeedFixed,
                                char* seed,
                                long unsigned int NumTests,
                                int machineReadable)
{

   if (NumTests < 2)
      fprintf(stderr, "\n\nA number of tests less than 2 is invalid !!!");
   initialize_pseudo_random_generator_seed(isSeedFixed, seed);

   unsigned char *pk= calloc(1, CRYPTO_PUBLICKEYBYTES);
   unsigned char *sk= calloc(1, CRYPTO_SECRETKEYBYTES);
   unsigned long long mlen=1024;
   unsigned char *plaintext = calloc(1, mlen);
   memset(plaintext,0xFF,mlen);
   unsigned char *ct= calloc(1, mlen+CRYPTO_BYTES);
   unsigned char *plaintextDec = calloc(1, CRYPTO_BYTES);
   unsigned long long clen=0;
   unsigned long long mlenDec=0;

   struct timespec beginning, end;
   long double time_elapsed_nanos;

   welford_t keygen_time_stat;
   welford_init(&keygen_time_stat);

   for (unsigned trials = 0; trials < NumTests; trials++) {
      clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &beginning);
      crypto_encrypt_keypair(pk, sk);
      clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end);
      time_elapsed_nanos = compute_time_interval(&beginning, &end)*1E+3;
      welford_update(&keygen_time_stat, time_elapsed_nanos);
   }

   welford_t encrypt_time_stat;
   welford_t decrypt_time_stat;
   welford_init(&encrypt_time_stat);
   welford_init(&decrypt_time_stat);

   long int ok, decodeOk = 0, memcmpOk = 0;

   for (unsigned trials = 0; trials < NumTests; trials++) {
      for (unsigned i = 0 ; i < mlen; i++){
            plaintext[i] = rand();
      }
      clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &beginning);
      crypto_encrypt( ct,&clen,plaintext,mlen,pk );
      clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end);
      time_elapsed_nanos = compute_time_interval(&beginning,&end)*1E+3;
      welford_update(&encrypt_time_stat, time_elapsed_nanos);

      clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &beginning);
      ok = !crypto_encrypt_open(plaintextDec,&mlenDec,ct,clen,sk);
      clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end);

      decodeOk += ok;
      ok        = ( memcmp(plaintext,plaintextDec,mlenDec) == 0 );
      memcmpOk += ok;

      if (ok){
          time_elapsed_nanos = compute_time_interval(&beginning,&end)*1E+3;
          welford_update(&decrypt_time_stat, time_elapsed_nanos);
      }

   } // end for trials

   if (machineReadable) {
        print_results_machine_readable("TIME",
                                       keygen_time_stat,
                                       encrypt_time_stat,
                                       decrypt_time_stat);
   } else {
        print_results_human_readable(keygen_time_stat,
                                     encrypt_time_stat,
                                     decrypt_time_stat,
                                     NumTests,
                                     decodeOk,
                                     memcmpOk);
   }
#if (defined HIGH_PERFORMANCE_X86_64)
if (machineReadable) {
   collect_and_print_cycles_metrics(NumTests);
}   
#endif
   publicKey_deletion_McEliece((publicKeyMcEliece_t *) pk);
   free(pk);
   privateKey_deletion_McEliece((privateKeyMcEliece_t *) sk);
   free(sk);
   free(plaintext);
   free(ct);
   free(plaintextDec);

} // end test_PKC_mceliece_cca2_code

/*----------------------------------------------------------------------------*/
