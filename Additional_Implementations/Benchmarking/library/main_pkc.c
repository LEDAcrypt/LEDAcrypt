#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>    // C99 sqrtl(...) for sample variance computation
#include <assert.h>

#include "qc_ldpc_parameters.h"
#include "gf2x_limbs.h"
#include "gf2x_arith_mod_xPplusOne.h"

#include "test_pkc.h"
#include "api.h"


int main(int argc, char *argv[])
{
   if(argc == 1){
     printf("Benchmarking utility for LEDAcrypt\n");
     printf("Usage: %s <machine readable> [PRNG Seed]\n machine_readable: 0 for human readable results, 1 for machine readable\n PRNG Seed: optional seed for PRNG (provided as decimal number)\n",argv[0]);
     return -1;
   }
   print_PKC_mceliece_cca2_parameters(atoi(argv[1]));
   int fixed_seed = (argc == 3);
   test_PKC_mceliece_cca2_code(fixed_seed,argv[2], NUM_TESTS, atoi(argv[1]));

   /* NIST API compliance */
   unsigned char *pk= calloc(1, CRYPTO_PUBLICKEYBYTES);
   unsigned char *sk= calloc(1, CRYPTO_SECRETKEYBYTES);
   unsigned long long mlen=256;
   unsigned char *plaintext = calloc(1, mlen);
   memset(plaintext,0xFF,mlen);

   unsigned char *plaintextDec = calloc(1, CRYPTO_BYTES);
   unsigned char *ct= calloc(1, mlen+CRYPTO_BYTES);
   unsigned long long clen=0;
   unsigned long long mlenDec=0;

   fprintf(stderr,"Keygen returned: %d\n", crypto_encrypt_keypair( pk,sk ));
   fprintf(stderr,"Encrypt returned: %d\n", crypto_encrypt( ct,
           &clen,
           plaintext,
           mlen,
           pk ));
   fprintf(stderr,"Decrypt returned: %d\n", crypto_encrypt_open( plaintextDec,
           &mlenDec,
           ct,
           clen,
           sk ));
   fprintf(stderr,"Memcmp returned: %d\n", memcmp(plaintext,plaintextDec,mlen));

   free(pk);
   free(sk);
   free(plaintext);
   free(plaintextDec);
   free(ct);
   /* performance test */

   return EXIT_SUCCESS;
} // end main

/*----------------------------------------------------------------------------*/


