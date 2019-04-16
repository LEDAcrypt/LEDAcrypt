#include "test_mceliece_cca2.h"

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

#define __USE_POSIX199309
#include <time.h> // struct timespec; clock_gettime(...); CLOCK_REALTIME

/*----------------------------------------------------------------------------*/

void print_PKC_mceliece_cca2_parameters(void)
{

   fprintf(stderr,
           "\n -------------------------------------------------------------");
   fprintf(stderr,
           "\n |        QC-LDPC PKC Kobara-Imai gamma (KI-g) Scheme        |");
   fprintf(stderr,
           "\n |                                                           |");
   fprintf(stderr,"\n | N0:............................................%5d      |",
           N0);
   fprintf(stderr,"\n | P:.............................................%5d(b)   |",
           P);
   fprintf(stderr,
           "\n | K/N:...............................%5d/%5d = %2.1lf      |", K, N,
           K/((double)N));
   fprintf(stderr,
           "\n |                                                           |");
   fprintf(stderr,"\n | H circ. block weight DV:.......................%5d      |",
           DV);
# if N0 == 2
   fprintf(stderr,"\n | Q  -- circ. block weight ... M0:%5d                     |",
           M0);
   fprintf(stderr,"\n | Q  -- circ. block weight ... M1:%5d                     |",
           M1);
#elif N0 == 3
   fprintf(stderr,"\n | Q  -- circ. block weight ... M0:%5d                     |",
           M0);
   fprintf(stderr,"\n | Q  -- circ. block weight ... M1:%5d                     |",
           M1);
   fprintf(stderr,"\n | Q  -- circ. block weight ... M2:%5d                     |",
           M2);
#elif N0 == 4
   fprintf(stderr,"\n | Q  -- circ. block weight ...............M0:%5d          |",
           M0);
   fprintf(stderr,"\n | Q  -- circ. block weight ...............M1:%5d          |",
           M1);
   fprintf(stderr,"\n | Q  -- circ. block weight ...............M2:%5d          |",
           M2);
   fprintf(stderr,"\n | Q  -- circ. block weight ...............M3:%5d          |",
           M3);
#else
#error "Unsupported number of Q blocks"
#endif
   fprintf(stderr,"\n | numb. of errors. T:............................%5d      |",
           NUM_ERRORS_T);
   fprintf(stderr,
           "\n |                                                           |");
   fprintf(stderr,"\n | True-RNG material byte length:.................%5d (B)  |",
           TRNG_BYTE_LENGTH);
   fprintf(stderr,"\n | Hash digest byte length:.......................%5d (B)  |",
           HASH_BYTE_LENGTH);
   fprintf(stderr,"\n | Min. Ciphertext Overhead.......................%5lu (B)  |",
           N0*NUM_DIGITS_GF2X_ELEMENT*DIGIT_SIZE_B-MAX_BYTES_IN_IWORD);
    fprintf(stderr,"\n | Max. Ciphertext Overhead.......................%5d (B)  |",
           N0*NUM_DIGITS_GF2X_ELEMENT*DIGIT_SIZE_B);   fprintf(stderr,
           "\n |                                                           |");
   fprintf(stderr,
           "\n | private key size:....................%5u(B) = %2.1lf(KiB)  |",
           (unsigned int) sizeof(privateKeyMcEliece_t),
           sizeof(privateKeyMcEliece_t) / ((double) 1024));
   fprintf(stderr,
           "\n | public  key size:....................%5u(B) = %2.1lf(KiB)  |",
           (unsigned int) sizeof(publicKeyMcEliece_t),
           sizeof(publicKeyMcEliece_t) / ((double) 1024));
   fprintf(stderr,
           "\n |                                                           |");
   fprintf(stderr,
           "\n -------------------------------------------------------------\n");

   fprintf(stderr,
           "\n KEYSIZETABLE: %d/%d & %5lu & %5lu & %5lu & %5lu & %5d \\\\ \n",
           CATEGORY, DFR_SL_LEVEL,
           sizeof(privateKeyMcEliece_t),
           (N0*(DV+M)+ DV*M)*sizeof(uint32_t), /* expanded private */
           sizeof(publicKeyMcEliece_t), /* pubkey */
           N0*NUM_DIGITS_GF2X_ELEMENT*DIGIT_SIZE_B-MAX_BYTES_IN_IWORD,/* min ovh */
           N0*NUM_DIGITS_GF2X_ELEMENT*DIGIT_SIZE_B); /*max ovh*/ 
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

static
void online_mean_and_variance(long double *accMean,
                              long double *accVar,
                              long int accCounter,
                              long double x)
{
   long double delta = x - (*accMean);
   (*accMean) += delta/(accCounter);
   long double delta2 = x - (*accMean);
   (*accVar) += delta*delta2;

} // end online_mean_and_variance

void test_PKC_mceliece_cca2_code(int ac, char *av[],
                                 long unsigned int NumTests)
{

   if (NumTests < 2)
      fprintf(stderr, "\n\nA number of tests less than 2 is invalid !!!");

   initialize_pseudo_random_generator_seed(ac, av);


   struct timespec timerBegin_kg, timerEnd_kg;
   long double time_elapsed_nsec_kg;
   long double sm_keygen = 0.0L, sm2_keygen = 0.0L;

   unsigned char *pk= calloc(1, CRYPTO_PUBLICKEYBYTES);
   unsigned char *sk= calloc(1, CRYPTO_SECRETKEYBYTES);
   for (unsigned trials = 0; trials < NumTests; trials++) {
      clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &timerBegin_kg);
      crypto_encrypt_keypair(pk,sk);
      clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &timerEnd_kg);
      time_elapsed_nsec_kg = compute_time_interval(&timerBegin_kg, &timerEnd_kg)*1E+9;
      online_mean_and_variance(&sm_keygen,
                               &sm2_keygen,
                               1+trials,
                               time_elapsed_nsec_kg/(1E+6));
   }


   struct timespec timerBegin, timerEnd;
   long double time_elapsed_nanos,
        sm_enc = 0.0L, sm2_enc = 0.0L,
        sm_decode = 0.0L, sm2_decode = 0.0L,
        sm_cmp = 0.0L, sm2_cmp = 0.0L;

   long int ok, decodeOk = 0, memcmpOk = 0;
   unsigned long long mlen=1024;
   unsigned char *plaintext = calloc(1, mlen);
   memset(plaintext,0xFF,mlen);
   unsigned char *ct= calloc(1, mlen+CRYPTO_BYTES);

   unsigned char *plaintextDec = calloc(1, CRYPTO_BYTES);
   unsigned long long clen=0;
   unsigned long long mlenDec=0;

   for (unsigned trials = 0; trials < NumTests; trials++) {
      for (unsigned i = 0 ; i < mlen; i++){
            plaintext[i] = rand();
      }
      clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &timerBegin);
      crypto_encrypt( ct,&clen,plaintext,mlen,pk );
      clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &timerEnd);
      time_elapsed_nanos = compute_time_interval(&timerBegin,&timerEnd)*1E+9;
      online_mean_and_variance(&sm_enc,
                               &sm2_enc,
                               1+trials,
                               time_elapsed_nanos/(1E+6));
      clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &timerBegin);
      ok = !crypto_encrypt_open(plaintextDec,&mlenDec,ct,clen,sk);
      clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &timerEnd);
      time_elapsed_nanos = compute_time_interval(&timerBegin,&timerEnd)*1E+9;
      decodeOk += ok;
      if (ok) online_mean_and_variance(&sm_decode,
                                          &sm2_decode,
                                          decodeOk,
                                          time_elapsed_nanos/(1E+6));
//       fprintf(stderr,"\nMLEN: %llu MLENDEC: %llu\n", mlen,mlenDec);

      ok        = ( memcmp(plaintext,plaintextDec,mlenDec) == 0 );
      memcmpOk += ok;
      if (ok) online_mean_and_variance(&sm_cmp,
                                          &sm2_cmp,
                                          memcmpOk,
                                          time_elapsed_nanos/(1E+6));

   } // end for trials

   fprintf(stderr,"\n\nPerformance tests.");
   fprintf(stderr, "\n\n%12sNumber of tests: %10lu", " ",
           (long unsigned int)NumTests);
   fprintf(stderr, "\nAverage key generation time: %14.3Lf (+,- %.3Lf) millisec",
           sm_keygen, sqrtl(sm2_keygen/(NumTests-1)));
   fprintf(stderr,"\nNumber of correct decodings: %10ld \
                  \n  ---  ---      decryptions: %10ld",
           decodeOk, memcmpOk);

   fprintf(stderr,"\n%4sAverage encryption time: %14.3Lf (+,- %.3Lf) millisec",
           " ",
           sm_enc,
           sqrtl(sm2_enc/(NumTests-1))
          );
   if (decodeOk > 1)
      fprintf(stderr,
              "\n%4sAverage decryption time: %14.3Lf (+,- %.3Lf) millisec ... with successful decoding",
              " ",
              sm_decode,
              sqrtl(sm2_decode/(decodeOk-1))
             );
   else
      fprintf(stderr,"\nNo successful decoding !!!\n");

   if (memcmpOk > 1) {
      fprintf(stderr,
              "\n%4sAverage decryption time: %14.3Lf (+,- %.3Lf) millisec ... with successful decoding, and ptx and decrypted ctx match\n",
              " ",
              sm_cmp,
              sqrtl(sm2_cmp/(memcmpOk-1))
             );
   } else
      fprintf(stderr,"\nNo successful enc/dec match !!!\n");

   fprintf(stderr,
           "\n TIMETABLE: & %d & %.2Lf  ($\\pm$ %.2Lf ) & %.2Lf   ($\\pm$ %.2Lf ) & %.2Lf   ($\\pm$ %.2Lf )  &   %.2Lf  \\\\ \n",
           DFR_SL_LEVEL,
           sm_keygen, sqrtl(sm2_keygen/(NumTests-1)),
           sm_enc,    sqrtl(sm2_enc/(NumTests-1)),
           sm_decode, sqrtl(sm2_decode/(decodeOk-1)),
           sm_keygen+sm_enc+sm_decode);

   publicKey_deletion_McEliece((publicKeyMcEliece_t *) pk);
   free(pk);
   pk = NULL;
   privateKey_deletion_McEliece((privateKeyMcEliece_t *) sk);
   free(sk);
   sk = NULL;
   free(plaintext);
   free(ct);
   free(plaintextDec);

} // end test_PKC_mceliece_cca2_code

/*----------------------------------------------------------------------------*/
