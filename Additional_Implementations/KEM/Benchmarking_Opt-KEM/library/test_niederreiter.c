#include "test_niederreiter.h"

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

#define __USE_POSIX199309

#include <time.h> // struct timespec; clock_gettime(...); CLOCK_REALTIME
#include "api.h"

/*----------------------------------------------------------------------------*/

void print_KEM_parameters(void)
{

   fprintf(stderr,"\n --------------------------------------------------------");
   fprintf(stderr,
           "\n |%1sLEDAkem (QC-LDPC based key encapsulation mechanism) %1s|", " ", " ");
   fprintf(stderr,"\n --------------------------------------------------------");
   fprintf(stderr,"\n |%54s|"," ");
   fprintf(stderr,"\n |                CATEGORY:......%6d%17s|", CATEGORY," ");
   fprintf(stderr,"\n |                      N0:......%6d%17s|", N0," ");
   fprintf(stderr,"\n |                       P:......%6d(b)%14s|", P," ");
   fprintf(stderr,"\n |               Rate(K/N):......%8.1lf%15s|", K/((double)N),
           " ");
   fprintf(stderr,"\n |%54s|"," ");
   fprintf(stderr,"\n | H circ. block weight DV:......%6d%17s|", DV, " ");
#if N0 == 2
   fprintf(stderr,"\n | Q circ. block weight M0:......%6d%17s|", M0, " ");
   fprintf(stderr,"\n |                      M1:......%6d%17s|", M1, " ");
#elif N0 == 3
   fprintf(stderr,"\n | Q circ. block weight M0:......%6d%17s|", M0, " ");
   fprintf(stderr,"\n |                      M1:......%6d%17s|", M1, " ");
   fprintf(stderr,"\n |                      M2:......%6d%17s|", M2, " ");
#elif N0 == 4
   fprintf(stderr,"\n | Q circ. block weight M0:......%6d%17s|", M0, " ");
   fprintf(stderr,"\n |                      M1:......%6d%17s|", M1, " ");
   fprintf(stderr,"\n |                      M2:......%6d%17s|", M2, " ");
   fprintf(stderr,"\n |                      M3:......%6d%17s|", M3, " ");
#else
#error "Unsupported number of Q blocks"
#endif
   fprintf(stderr,"\n |%54s|"," ");
   fprintf(stderr,"\n |      number of errors T:......%6d%17s|", NUM_ERRORS_T,
           " ");
   fprintf(stderr,"\n |%54s|"," ");
   fprintf(stderr,
           "\n |        private key size:......%6lu(B) = %2lu.%1lu(KiB)%2s|",
           sizeof(privateKeyNiederreiter_t),
           ((unsigned long)(sizeof(privateKeyNiederreiter_t)/((double) 1024))),
           (unsigned long)(10*((sizeof(privateKeyNiederreiter_t)/((double) 1024))-
                               (unsigned long)(sizeof(privateKeyNiederreiter_t)/((double) 1024)))), " ");
   fprintf(stderr,
           "\n |        public  key size:......%6lu(B) = %2lu.%1lu(KiB)%2s|",
           sizeof(publicKeyNiederreiter_t),
           (unsigned long )(sizeof(publicKeyNiederreiter_t)/((double) 1024)),
           (unsigned long)((sizeof(publicKeyNiederreiter_t)/((double) 1024))-
                           (unsigned long )(sizeof(publicKeyNiederreiter_t)/((double) 1024)))*10," ");
   fprintf(stderr,"\n |   encapsulated key size:......%6u(B) = %2lu.%1lu(KiB)%2s|",
           NUM_DIGITS_GF2X_ELEMENT*DIGIT_SIZE_B,
           (unsigned long)(NUM_DIGITS_GF2X_ELEMENT*DIGIT_SIZE_B/((double) 1024)),
           (unsigned long)(10*((NUM_DIGITS_GF2X_ELEMENT*DIGIT_SIZE_B/((double) 1024))-
                               (unsigned long)(NUM_DIGITS_GF2X_ELEMENT*DIGIT_SIZE_B/((double) 1024))))," ");
   fprintf(stderr,"\n |%54s|"," ");
   fprintf(stderr,"\n --------------------------------------------------------");
     fprintf(stderr,"\n KTABLE: %d/%d & %5lu & %5lu & %5lu & %5d & %5d \\ \n",
          CATEGORY, N0,
          sizeof(privateKeyNiederreiter_t), 
          (N0*(DV+M)+ DV*M)*sizeof(uint32_t),
          sizeof(publicKeyNiederreiter_t),
          NUM_DIGITS_GF2X_ELEMENT*DIGIT_SIZE_B,
          HASH_BYTE_LENGTH );
   
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

void test_KEM_niederreiter_code(int ac, char *av[], long unsigned int NumTests)
{

   if (NumTests < 2)
      fprintf(stderr, "\n\nA number of tests less than 2 is invalid !!!");

   initialize_pseudo_random_generator_seed(ac, av);

   unsigned char  *pk       = calloc(1, CRYPTO_PUBLICKEYBYTES),
                  *sk       = calloc(1, CRYPTO_SECRETKEYBYTES),
                  *ss       = calloc(1, CRYPTO_BYTES),
                  *decap_ss = calloc(1, CRYPTO_BYTES),
                  *ct       = calloc(1, CRYPTO_CIPHERTEXTBYTES);

   struct timespec timerBegin_kg, timerEnd_kg;
   long double time_elapsed_nsec_kg;
   long double sm_keygen = 0.0L, sm2_keygen = 0.0L;
   for (unsigned trials = 0; trials < NumTests; trials++) {

      clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &timerBegin_kg);
      crypto_kem_keypair(pk, sk);
      clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &timerEnd_kg);
      time_elapsed_nsec_kg = compute_time_interval(&timerBegin_kg, &timerEnd_kg)*1E+9;
      online_mean_and_variance(&sm_keygen,
                               &sm2_keygen,
                               1+trials,
                               time_elapsed_nsec_kg/(1E+6));
   }
    fprintf(stderr,
            "\n STABLE: %d/%d & %4.3Lf ($\\pm$ %.3Lf)" ,
            CATEGORY, 
            N0,
            sm_keygen, 
            sqrtl(sm2_keygen/(NumTests-1)));


   struct timespec timerBegin, timerEnd;
   long double time_elapsed_nanos,
        sm_enc = 0.0L, sm2_enc = 0.0L,
        sm_decode = 0.0L, sm2_decode = 0.0L,
        sm_cmp = 0.0L, sm2_cmp = 0.0L;
   long int ok, decodeOk = 0, memcmpOk = 0;

   crypto_kem_keypair(pk, sk);

   for (unsigned trials = 0; trials < NumTests; trials++) {


      clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &timerBegin);
      crypto_kem_enc(ct, ss, pk);
      clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &timerEnd);
      time_elapsed_nanos = compute_time_interval(&timerBegin,&timerEnd)*1E+9;
      online_mean_and_variance(&sm_enc,
                               &sm2_enc,
                               1+trials,
                               time_elapsed_nanos/(1E+6));

      clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &timerBegin);
      ok =  crypto_kem_dec(decap_ss, ct, sk);
      clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &timerEnd);
      time_elapsed_nanos = compute_time_interval(&timerBegin,&timerEnd)*1E+9;
      decodeOk += (ok == 0); /*NIST API has 0 = success, 1= failure*/

      if (ok == 0) online_mean_and_variance(&sm_decode,
                                          &sm2_decode,
                                          decodeOk,
                                          time_elapsed_nanos/(1E+6));
      ok        = (memcmp(ss, decap_ss, CRYPTO_BYTES) == 0);
      memcmpOk += ok;
      if (ok) online_mean_and_variance(&sm_cmp,
                                          &sm2_cmp,
                                          memcmpOk,
                                          time_elapsed_nanos/(1E+6));
   } // end for trials

    fprintf(stderr,
          " & %4.3Lf ($\\pm$ %.3Lf) & %4.3Lf ($\\pm$ %.3Lf) & %4.3Lf \\\\ \n",
          sm_enc, 
          sqrtl(sm2_enc/(NumTests-1)),
          sm_decode, 
          sqrtl(sm2_decode/(decodeOk-1)),
          sm_enc+sm_keygen+sm_decode );
   
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
      fprintf(stderr,"\nNo successful decoding !!!");

   if (memcmpOk > 1)
      fprintf(stderr,
              "\n%4sAverage decryption time: %14.3Lf (+,- %.3Lf) millisec ... with successful decoding, and ptx and decrypted ctx match",
              " ",
              sm_cmp,
              sqrtl(sm2_cmp/(memcmpOk-1))
             );
   else
      fprintf(stderr,"\nNo successful enc/dec match !!!");

   publicKey_deletion_niederreiter((publicKeyNiederreiter_t *) pk);
   free(pk);
   pk = NULL;
   privateKey_deletion_niederreiter((privateKeyNiederreiter_t *)sk);
   free(sk);
   sk = NULL;
   free(ss);
   free(decap_ss);
   free(ct);
   
} // end test_KEM_niederreiter_code

/*----------------------------------------------------------------------------*/
