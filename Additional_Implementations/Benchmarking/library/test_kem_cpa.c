#include "test_kem.h"

#include "qc_ldpc_parameters.h"
#include "gf2x_limbs.h"
#include "gf2x_arith_mod_xPplusOne.h"

#include "niederreiter.h"
#include "niederreiter_keygen.h"
#include "niederreiter_encrypt.h"
#include "niederreiter_decrypt.h"

#include "rng.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>    // C99 sqrtl(...) for sample variance computation


#include <time.h> // struct timespec; clock_gettime(...); CLOCK_REALTIME
#include "api.h"
#include "timing_and_stat.h"

/*----------------------------------------------------------------------------*/

void print_KEM_parameters(int machine_readable)
{
if(machine_readable == 0) {
   fprintf(stderr,"\n --------------------------------------------------------");
   fprintf(stderr,
           "\n  %1sLEDAcrypt KEM-CPA (key encapsulation mechanism \n %1s with emphemeral public/private keypair) %1s ", " ", " ", " ");
   fprintf(stderr,"\n --------------------------------------------------------");
   fprintf(stderr,"\n  %54s "," ");
   fprintf(stderr,"\n                  CATEGORY:......%6d%17s ", CATEGORY," ");
   fprintf(stderr,"\n                        N0:......%6d%17s ", N0," ");
   fprintf(stderr,"\n                         P:......%6d(b)%14s ", P," ");
   fprintf(stderr,"\n                 Rate(K/N):......%8.1lf%15s ", K/((double)N),
           " ");
   fprintf(stderr,"\n  %54s "," ");
   fprintf(stderr,"\n   H circ. block weight V:......%6d%17s ", V, " ");

   fprintf(stderr,"\n  %54s "," ");
   fprintf(stderr,"\n        number of errors T:......%6d%17s ", NUM_ERRORS_T,
           " ");
   fprintf(stderr,"\n  %54s "," ");
   fprintf(stderr,
           "\n          private key size:......%6lu(B) = %2.1lf(KiB)%2s ",
                                               (long unsigned int) (CRYPTO_SECRETKEYBYTES),
                                               ((CRYPTO_SECRETKEYBYTES)/((double) 1024)),
                                               " " );
   fprintf(stderr,
           "\n          public  key size:......%6lu(B) = %2.1lf(KiB)%2s ",
           (long unsigned int) (CRYPTO_PUBLICKEYBYTES),
           ((CRYPTO_PUBLICKEYBYTES)/((double) 1024))," ");
   fprintf(stderr,"\n          encapsulated-key size:.%6lu(B) = %2.1lf(KiB)%2s ",
           (long unsigned int)((CRYPTO_CIPHERTEXTBYTES)),
           (((long unsigned int)((CRYPTO_CIPHERTEXTBYTES))/((double) 1024))), "  ");
   fprintf(stderr,"\n          shared secret:.......%8lu(B) = %2.1lf(KiB)%2s ",
           (long unsigned int)((CRYPTO_BYTES)),
           (((long unsigned int)((CRYPTO_BYTES))/((double) 1024))), "  ");
   fprintf(stderr,"\n  %54s "," ");
   fprintf(stderr,"\n --------------------------------------------------------");
} else {
   fprintf(stderr,"\n SIZE: %d & %d & %5lu & %5lu & %5lu & %5lu & %5lu \\\\\n",
          CATEGORY, N0,
          (long unsigned int) (CRYPTO_SECRETKEYBYTES),
          (long unsigned int) (CRYPTO_PUBLICKEYBYTES),
          (long unsigned int) ((CRYPTO_CIPHERTEXTBYTES)),
          (long unsigned int) ((CRYPTO_BYTES)),
          (long unsigned int) (CRYPTO_PUBLICKEYBYTES+CRYPTO_CIPHERTEXTBYTES));
}
} // end print_KEM_parameters

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
         "\n %s: %d & %d ",
         tag,
         CATEGORY,
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
         (welford_mean(keygen_stats)+ welford_mean(encrypt_stats) + welford_mean(decrypt_stats)));
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
         "\n %s: %d & %d ",
         tag,
         CATEGORY,
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
         (welford_mean(keygen_stats)+ welford_mean(encrypt_stats) + welford_mean(decrypt_stats))/1E+3,
         (welford_mean(keygen_stats)+
          welford_mean(encrypt_stats)+
          welford_mean(decrypt_stats)+
          (CRYPTO_PUBLICKEYBYTES+CRYPTO_CIPHERTEXTBYTES)*METRIC_ALPHA_FACTOR)/1E+3);
}

void collect_and_print_cycles_metrics(long unsigned int NumTests){
       unsigned char  *pk       = calloc(1, CRYPTO_PUBLICKEYBYTES),
                  *sk       = calloc(1, CRYPTO_SECRETKEYBYTES),
                  *ss       = calloc(1, CRYPTO_BYTES),
                  *decap_ss = calloc(1, CRYPTO_BYTES),
                  *ct       = calloc(1, CRYPTO_CIPHERTEXTBYTES);

   uint64_t beginning, end;

   welford_t keygen_time_stat;
   welford_init(&keygen_time_stat);

   for (unsigned trials = 0; trials < NumTests; trials++) {
      beginning = x86_64_rtdsc();
      crypto_kem_keypair(pk, sk);
      end = x86_64_rtdsc();
      welford_update(&keygen_time_stat, (end-beginning));
   }

   welford_t encrypt_time_stat;
   welford_t decrypt_time_stat;
   welford_init(&encrypt_time_stat);
   welford_init(&decrypt_time_stat);
   long int ok, decodeOk = 0, memcmpOk = 0;

   crypto_kem_keypair(pk, sk);
   for (unsigned trials = 0; trials < NumTests; trials++) {
      beginning = x86_64_rtdsc();
      crypto_kem_enc(ct, ss, pk);
      end = x86_64_rtdsc();
      welford_update(&encrypt_time_stat, (end-beginning));

      beginning = x86_64_rtdsc();
      ok =  crypto_kem_dec(decap_ss, ct, sk);
      end = x86_64_rtdsc();

      decodeOk += (ok == 0); /*NIST API has 0 = success, 1= failure */
      ok        = (memcmp(ss, decap_ss, CRYPTO_BYTES) == 0);
      memcmpOk += ok;
      if (ok){
          welford_update(&decrypt_time_stat, (end-beginning));
      }
   } // end for trials
   print_cycles_metrics("CYCLES",
                        keygen_time_stat,
                        encrypt_time_stat,
                        decrypt_time_stat);
   publicKey_deletion_niederreiter((publicKeyNiederreiter_t *) pk);
   free(pk);
   privateKey_deletion_niederreiter((privateKeyNiederreiter_t *)sk);
   free(sk);
   free(ss);
   free(decap_ss);
   free(ct);
}
#endif

/*----------------------------------------------------------------------------*/

void test_KEM_niederreiter_code(int isSeedFixed,
                                char* seed,
                                long unsigned int NumTests,
                                int machineReadable)
{
   initialize_pseudo_random_generator_seed(isSeedFixed, seed);
   unsigned char  *pk       = calloc(1, CRYPTO_PUBLICKEYBYTES),
                  *sk       = calloc(1, CRYPTO_SECRETKEYBYTES),
                  *ss       = calloc(1, CRYPTO_BYTES),
                  *decap_ss = calloc(1, CRYPTO_BYTES),
                  *ct       = calloc(1, CRYPTO_CIPHERTEXTBYTES);

   struct timespec beginning, end;
   long double time_elapsed_nanos;

   welford_t keygen_time_stat;
   welford_init(&keygen_time_stat);

   for (unsigned trials = 0; trials < NumTests; trials++) {
      clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &beginning);
      crypto_kem_keypair(pk, sk);
      clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end);
      time_elapsed_nanos = compute_time_interval(&beginning, &end)*1E+3;
      welford_update(&keygen_time_stat, time_elapsed_nanos);
   }

   welford_t encrypt_time_stat;
   welford_t decrypt_time_stat;
   welford_init(&encrypt_time_stat);
   welford_init(&decrypt_time_stat);
   long int ok, decodeOk = 0, memcmpOk = 0;

   crypto_kem_keypair(pk, sk);
   for (unsigned trials = 0; trials < NumTests; trials++) {

      clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &beginning);
      crypto_kem_enc(ct, ss, pk);
      clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end);
      time_elapsed_nanos = compute_time_interval(&beginning,&end)*1E+3;
      welford_update(&encrypt_time_stat, time_elapsed_nanos);

      clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &beginning);
      ok =  crypto_kem_dec(decap_ss, ct, sk);
      clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end);

      decodeOk += (ok == 0); /*NIST API has 0 = success, 1= failure */
      ok        = (memcmp(ss, decap_ss, CRYPTO_BYTES) == 0);
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
#if (defined HIGH_PERFORMANCE_X86_64)
        collect_and_print_cycles_metrics(NumTests);
#endif
   } else {
     print_results_human_readable(keygen_time_stat,
                                  encrypt_time_stat,
                                  decrypt_time_stat,
                                  NumTests,
                                  decodeOk,
                                  memcmpOk);
   }
   publicKey_deletion_niederreiter((publicKeyNiederreiter_t *) pk);
   free(pk);
   privateKey_deletion_niederreiter((privateKeyNiederreiter_t *)sk);
   free(sk);
   free(ss);
   free(decap_ss);
   free(ct);

} // end test_KEM_niederreiter_code

/*----------------------------------------------------------------------------*/
