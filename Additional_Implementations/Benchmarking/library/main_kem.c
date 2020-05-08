#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "test_kem.h"
#include "api.h"

/*----------------------------------------------------------------------------*/

int main(int argc, char *argv[])
{
   if(argc == 1){
     printf("Benchmarking utility for LEDAcrypt\n");
     printf("Usage: %s <machine readable> [PRNG Seed]\n machine_readable: 0 for human readable results, 1 for machine readable\n PRNG Seed: optional seed for PRNG (provided as decimal number)\n",argv[0]);
     return -1;
   }
   print_KEM_parameters(atoi(argv[1]));
   int fixed_seed = (argc == 3);
   test_KEM_niederreiter_code(fixed_seed,argv[2], NUM_TESTS, atoi(argv[1]));

   unsigned char  *pk       = calloc(1, CRYPTO_PUBLICKEYBYTES),
                  *sk       = calloc(1, CRYPTO_SECRETKEYBYTES),
                  *ss       = calloc(1, CRYPTO_BYTES),
                  *decap_ss = calloc(1, CRYPTO_BYTES),
                  *ct       = calloc(1, CRYPTO_CIPHERTEXTBYTES);

   fprintf(stderr,"\n\nTesting NIST API.\n");
   fprintf(stderr,"\ncrypto_kem_keypair returned: %10d",
           crypto_kem_keypair(pk, sk));
   fprintf(stderr,"\n    crypto_kem_enc returned: %10d",
           crypto_kem_enc(ct, ss, pk));
   fprintf(stderr,"\n    crypto_kem_dec returned: %10d",
           crypto_kem_dec(decap_ss, ct, sk));
   int cmp = memcmp(ss, decap_ss, CRYPTO_BYTES);
   fprintf(stderr,"\n            memcmp returned: %10d", cmp);
   if (cmp == 0)
      fprintf(stderr,
              " ... shared_secret and decapsulated_shared_secret match!");
   else
      fprintf(stderr,
              " ... shared_secret and decapsulated_shared_secret do not match!");
   fprintf(stderr,"\n\n");

   free(pk);
   free(sk);
   free(ss);
   free(decap_ss);
   free(ct);

   return EXIT_SUCCESS;

} // end main

/*----------------------------------------------------------------------------*/


