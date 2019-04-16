/**
 *
 * @version 2.0 (March 2019)
 *
 * Reference ISO-C11 Implementation of the LEDAcrypt PKC cipher using GCC built-ins.
 *
 * In alphabetical order:
 *
 * @author Marco Baldi <m.baldi@univpm.it>
 * @author Alessandro Barenghi <alessandro.barenghi@polimi.it>
 * @author Franco Chiaraluce <f.chiaraluce@univpm.it>
 * @author Gerardo Pelosi <gerardo.pelosi@polimi.it>
 * @author Paolo Santini <p.santini@pm.univpm.it>
 *
 * This code is hereby placed in the public domain.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHORS ''AS IS'' AND ANY EXPRESS
 * OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHORS OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
 * BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 * WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
 * OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
 * EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 **/

#include "mceliece_keygen.h"
#include "H_Q_matrices_generation.h"
#include "gf2x_arith_mod_xPplusOne.h"
#include "rng.h"
#include "dfr_test.h"

/*----------------------------------------------------------------------------*/

void mceliece_keygen(publicKeyMcEliece_t   *const pk,
                     privateKeyMcEliece_t *const sk,
                     AES_XOF_struct *keys_expander)
{
   POSITION_T HtrPosOnes[N0][DV];
   POSITION_T HPosOnes[N0][DV];
   POSITION_T QPosOnes[N0][M];
   POSITION_T LPosOnes[N0][DV*M];

   int is_L_full;
   int isDFRok;
   sk->rejections = (int8_t) 0;
   do {
     generateHPosOnes_HtrPosOnes(HPosOnes,
                                 HtrPosOnes,
                                 keys_expander);

     generateQsparse(QPosOnes,
                     keys_expander);
     for (int i = 0; i < N0; i++) {
        for (int j = 0; j< DV*M; j++) {
           LPosOnes[i][j]=INVALID_POS_VALUE;
        }
     }

     POSITION_T auxPosOnes[DV*M];
     unsigned char processedQOnes[N0] = {0};
     for (int colQ = 0; colQ < N0; colQ++) {
        for (int i = 0; i < N0; i++) {
           gf2x_mod_mul_sparse(DV*M, auxPosOnes,
                               DV, HPosOnes[i],
                               qBlockWeights[i][colQ], QPosOnes[i]+processedQOnes[i]);
           gf2x_mod_add_sparse(DV*M, LPosOnes[colQ],
                               DV*M, LPosOnes[colQ],
                               DV*M, auxPosOnes);
           processedQOnes[i] += qBlockWeights[i][colQ];
        }
     }
     is_L_full = 1;
     for (int i = 0; i < N0; i++) {
        is_L_full = is_L_full && (LPosOnes[i][DV*M-1] != INVALID_POS_VALUE);
     }
     sk->rejections = sk->rejections + 1;
    
     if(is_L_full) {isDFRok = DFR_test(LPosOnes);}
   } while(!is_L_full || !isDFRok );
   sk->rejections = sk->rejections - 1;

   DIGIT Ln0dense[NUM_DIGITS_GF2X_ELEMENT] = {0x00};
   for(int j = 0; j < DV*M; j++) {
      if(LPosOnes[N0-1][j] != INVALID_POS_VALUE)
         gf2x_set_coeff(Ln0dense,LPosOnes[N0-1][j],1);
   }

   DIGIT Ln0Inv[NUM_DIGITS_GF2X_ELEMENT] = {0x00};
   gf2x_mod_inverse(Ln0Inv, Ln0dense);

   for (int i = 0; i < N0-1; i++) {
      gf2x_mod_mul_dense_to_sparse(pk->Mtr+i*NUM_DIGITS_GF2X_ELEMENT,
                                   Ln0Inv,
                                   LPosOnes[i],
                                   DV*M);
   }

   for (int i = 0; i < N0-1; i++) {
      gf2x_transpose_in_place(pk->Mtr+i*NUM_DIGITS_GF2X_ELEMENT);
   }
} // end mceliece_keygen


/*----------------------------------------------------------------------------*/
/* Implementation that should never be optimized out by the compiler */

static inline void zeroize( void *v, size_t n )
{
   volatile unsigned char *p = v;
   while( n-- ) *p++ = 0;
} // end zeroize

void publicKey_deletion_McEliece(publicKeyMcEliece_t   *const pk)
{
   zeroize(pk->Mtr, (N0-1)*NUM_DIGITS_GF2X_ELEMENT*DIGIT_SIZE_B);
} // publicKey_deletion_McEliece


/*----------------------------------------------------------------------------*/

void privateKey_deletion_McEliece(privateKeyMcEliece_t *const sk)
{
   zeroize(sk->prng_seed, TRNG_BYTE_LENGTH);
} // privateKey_deletion_McEliece

/*----------------------------------------------------------------------------*/
