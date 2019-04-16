/**
 *
 * <bf_decoding.c>
 *
 * @version 2.0 (March 2019)
 *
 * Reference ISO-C11 Implementation of the LEDAcrypt KEM-LT cipher using GCC built-ins.
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
#include "bf_decoding.h"
#include "gf2x_arith_mod_xPplusOne.h"
#include <string.h>
#include <assert.h>

#define ROTBYTE(a)   ( (a << 8) | (a >> (DIGIT_SIZE_b - 8)) )
#define ROTUPC(a)   ( (a >> 8) | (a << (DIGIT_SIZE_b - 8)) )
#define ROUND_UP(amount, round_amt) ( ((amount+round_amt-1)/round_amt)*round_amt )

int thresholds[2] = {B0 , (DV*M)/2+1};

#ifdef HIGH_PERFORMANCE_X86_64
int bf_decoding(DIGIT out[], // N0 polynomials
                const POSITION_T HtrPosOnes[N0][DV],
                const POSITION_T QtrPosOnes[N0][M],
                DIGIT privateSyndrome[]  //  1 polynomial
               )
{
#if P < 64
#error The circulant block size should exceed 64
#endif

   /* The unsatisfied parity checks are kept as N0 vectors , each one having
    * a length rounded up to the closest multiple of the coalescing factor 
    * DIGIT_SIZE_B, plus one entire digit, and with the last digit padded with zeroes.
    * This strategy allows to replicate an entire digit at the end of the computed UPCs 
    * effectively removing altogether the (key-dependent) checks which would happen 
    * in the correlation computation process */
#define SIZE_OF_UPC_VECTORIZED_READ 32
#define SIZE_OF_UPC_VECTORIZED_WRITE 64
   uint8_t unsatParityChecks[N0][ ROUND_UP(P,SIZE_OF_UPC_VECTORIZED_READ)+SIZE_OF_UPC_VECTORIZED_READ] = {{0}};
   POSITION_T currQBitPos[M];
   /* syndrome is endowed with cyclic padding in the leading word to avoid 
    * boundary checks. The pad should be at least as long as the bit-len of one 
    * vector of bits which is read during UPC computation */
   DIGIT currSyndrome[NUM_DIGITS_GF2X_ELEMENT+1];
   int check;
   int iteration = 0;

   do {
      gf2x_copy(currSyndrome+1, privateSyndrome);

/*position of the first set bit in the word, counting 63 ... 0*/
#if (MSb_POSITION_IN_MSB_DIGIT_OF_ELEMENT == (DIGIT_SIZE_b-1))
      currSyndrome[0] = currSyndrome[NUM_DIGITS_GF2X_ELEMENT];
#else
       currSyndrome[1] |= (currSyndrome[NUM_DIGITS_GF2X_ELEMENT] << (MSb_POSITION_IN_MSB_DIGIT_OF_ELEMENT+1) );
       currSyndrome[0] = currSyndrome[NUM_DIGITS_GF2X_ELEMENT] >> (DIGIT_SIZE_b- (MSb_POSITION_IN_MSB_DIGIT_OF_ELEMENT+1));
#endif

#define VECTYPE DIGIT
#define VECTYPE_SIZE_b DIGIT_SIZE_b
#define VECTYPE_SIZE_B DIGIT_SIZE_B
#define BIT_SIZE_UPC 8
#define VEC_COMB_MASK 0x0101010101010101ULL

    VECTYPE upcMat[VECTYPE_SIZE_b/BIT_SIZE_UPC];
    VECTYPE packedSynBits;
    for (int i = 0; i < N0; i++) {
       for (int valueIdx = 0; valueIdx < P; valueIdx = valueIdx + DIGIT_SIZE_b) { /* this index will be vectorized*/
          memset(upcMat, 0 , VECTYPE_SIZE_B*(VECTYPE_SIZE_b/BIT_SIZE_UPC));
          /* this fetches DIGIT_SIZE_b bits from each Htrpos, packed, and adds them to the 64 upc counters in upcmat */
          for(int HtrOneIdx = 0; HtrOneIdx < DV; HtrOneIdx++) {
                POSITION_T basePos = (HtrPosOnes[i][HtrOneIdx]+valueIdx);
                basePos = basePos %P ;
                /* lsb here is the one in base pos, others are subseq*/
                packedSynBits = gf2x_get_DIGIT_SIZE_coeff_vector_boundless(currSyndrome,basePos); 
                for(int upcMatRow = 0; upcMatRow < VECTYPE_SIZE_b/BIT_SIZE_UPC; upcMatRow++){
                    upcMat[upcMatRow] += packedSynBits & VEC_COMB_MASK;
                    packedSynBits = packedSynBits >> 1;
                }

              }/* end of for computing 64 upcs*/
         /* commit computed UPCs in the upc vector, in the proper order. 
          * UPCmat essentially needs transposition and linearization by row,
          * starting from the last row */
          for(int upcMatCol = 0; upcMatCol < VECTYPE_SIZE_b/BIT_SIZE_UPC; upcMatCol++){
              VECTYPE upcBuf = 0;
              for(int upcMatRow = 0; upcMatRow < VECTYPE_SIZE_b/BIT_SIZE_UPC; upcMatRow++){
               uint8_t matByte = upcMat[upcMatRow];
#if __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
               upcBuf |=  ((VECTYPE)(matByte)) << (8*upcMatRow);
#else
               upcBuf |=  (upcBuf << 8) + ((VECTYPE)(matByte));
#endif
                  upcMat[upcMatRow] = upcMat[upcMatRow] >> 8;
              }
              VECTYPE* vp = (VECTYPE *)(&unsatParityChecks[i][valueIdx+8*upcMatCol]);
              *(vp) = upcBuf;
          } 
         }
     }
     
      /* circular padding of unsatisfiedParityChecks so that the vector 
       * correlation computation does not need to wraparound loads*/
      for (int i = 0; i < N0; i++) {
          for(int j = 0; j < SIZE_OF_UPC_VECTORIZED_READ; j++)
           unsatParityChecks[i][P+j] = unsatParityChecks[i][j];
      }

      /* iteration based threshold determination*/
      int corrt_syndrome_based= thresholds[iteration];

      //Computation of correlation  with a full Q matrix

      for (int i = 0; i < N0; i++) {
            uint8_t qBlockIdxes[M];
            uint16_t qBlockOffsets[M];
            uint8_t qBlockIdxesCursor = 0, linearized_cursor = 0;
            int endQblockIdx = 0;
            for (int blockIdx = 0; blockIdx < N0; blockIdx++) {
              endQblockIdx += qBlockWeights[blockIdx][i];
              for (; linearized_cursor < endQblockIdx; linearized_cursor++) {
                qBlockIdxes[linearized_cursor] = qBlockIdxesCursor;
                qBlockOffsets[linearized_cursor] =  (uint16_t) qBlockIdxesCursor;
              }
              qBlockIdxesCursor++;
            }
         for (int j = 0; j <= P; j+=SIZE_OF_UPC_VECTORIZED_READ) {
            int currQoneIdx; // position in the column of QtrPosOnes[][...]
 
            /* vector of correlations computed in a single sweep of AVX2 vectorized
             * computation.  Stored with type punning as a 4 item uint64_t
             * vectorized corrs. */
            uint64_t correlation = 0; 
            uint64_t correlation_avx[4] = {0}; 
            __m256i vecCorr = _mm256_setzero_si256();
            __m256i vecUpcToAdd;
            /* a single correlation value is the sum of m "block circulant positioned" values of UPC */
            /* I need a 8x storage of UPCs which I can load as 8x bytes from upc vector starting from
             * the block circulant position. 
             * A trivial implementation will thus have m variables containing 8 consecutive UPCs.
             * 
             * A more reasonable one uses a single variable which gets filled 
             * with 8x UPCs per shot, bar circulant cases related to 8 consecutive correlation
             * computations, and a single circulant "one" from Q (ideal opts will unroll). 
             * Since the outer loop (j indexed, at line 183) is vectorized, a trailer of the last P%DIGIT_SIZE_b should be added 
             * (or not: a single vector computation discarding part of the results in correlation is equally fine) */

            /* This vectorization strategy will benefit from having a vector with the block circulant pos.s in Q
             * and another one with the block idxes. Both can be vector packed: pos.s are 16 bit wide, block idx.es are 2b wide
             * FTTB I will generate both of them at runtime, but the latter can be generated as a #define.*/

            for (currQoneIdx = 0; currQoneIdx < M; currQoneIdx++) {
                int tmp = QtrPosOnes[i][currQoneIdx]+j; /* this must turn into a vector load of 8x UPCs with consecutive j indices */
                tmp = tmp % P ;  /* this will need to be checked before and have a split vector load */
                currQBitPos[currQoneIdx] = tmp;  /* only a single base index must be kept per vector load*/
                vecUpcToAdd = _mm256_lddqu_si256((__m256i *) (&unsatParityChecks[qBlockOffsets[currQoneIdx]][tmp]));
                vecCorr = _mm256_add_epi8(vecUpcToAdd,vecCorr);
            }
            _mm256_storeu_si256( (__m256i *) correlation_avx,vecCorr);

            
            for(int vecCorrDeltaJ = 0 ; vecCorrDeltaJ < SIZE_OF_UPC_VECTORIZED_READ/DIGIT_SIZE_B; vecCorrDeltaJ++) { 
            /* Correlation based flipping */
            correlation = correlation_avx[vecCorrDeltaJ];
               for (int delta_j = 0 ; 
                    delta_j < DIGIT_SIZE_B && (j+vecCorrDeltaJ*DIGIT_SIZE_B+delta_j < P) ; 
                    delta_j++){
#if __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
                  uint64_t single_correl = ( correlation >> (8*delta_j) ) & 0xFF; 
#else          
                  uint64_t single_correl = ( correlation >> (8*(DIGIT_SIZE_B-1-delta_j)) ) & 0xFF; 
#endif         
                  if ((single_correl > corrt_syndrome_based)) { /* correlation will always be null for garbage bits */
                     gf2x_toggle_coeff(out+NUM_DIGITS_GF2X_ELEMENT*i, (j+vecCorrDeltaJ*DIGIT_SIZE_B+delta_j)%P); /*this can be vectorized if a vector of correlations are available*/
                     for (int v = 0; v < M; v++) { /* the same goes here */
                        unsigned syndromePosToFlip;
                        for (int HtrOneIdx = 0; HtrOneIdx < DV; HtrOneIdx++) {
                           syndromePosToFlip = (HtrPosOnes[qBlockIdxes[v]][HtrOneIdx] + 
                                                (currQBitPos[v]+vecCorrDeltaJ*DIGIT_SIZE_B+delta_j)%P );
                           syndromePosToFlip = syndromePosToFlip >= P ? syndromePosToFlip - P : syndromePosToFlip;
                           gf2x_toggle_coeff(privateSyndrome, syndromePosToFlip);
                        }
                     } // end for v
                  } // end if
               } // end for flipping correlations exceeding threshold
            } // end of for selecting part of correlation vector to analyze
         } // end for j
      } // end for i

      iteration = iteration + 1;
      check = 0;
      while (check < NUM_DIGITS_GF2X_ELEMENT && privateSyndrome[check] == 0) {
          check++;
      }

   } while (iteration < ITERATIONS_MAX && check < NUM_DIGITS_GF2X_ELEMENT);

   return (check == NUM_DIGITS_GF2X_ELEMENT);
}  // end QdecodeSyndromeThresh_bitFlip_sparse

/*C fallback if there are no AVX2*/
#else

int bf_decoding(DIGIT out[], // N0 polynomials
                const POSITION_T HtrPosOnes[N0][DV],
                const POSITION_T QtrPosOnes[N0][M],
                DIGIT privateSyndrome[]  //  1 polynomial
               )
{
#if P < 64
#error The circulant block size should exceed 64
#endif

   /* The unsatisfied parity checks are kept as N0 vectors , each one having
    * a length rounded up to the closest multiple of the coalescing factor 
    * DIGIT_SIZE_B, plus one entire digit, and with the last digit padded with zeroes.
    * This strategy allows to replicate an entire digit at the end of the computed UPCs 
    * effectively removing altogether the (key-dependent) checks which would happen 
    * in the correlation computation process */
#define SIZE_OF_UPC_VECTORIZED_READ 32
#define SIZE_OF_UPC_VECTORIZED_WRITE 64
   uint8_t unsatParityChecks[N0][ ROUND_UP(P,SIZE_OF_UPC_VECTORIZED_READ)+SIZE_OF_UPC_VECTORIZED_READ] = {{0}};
   POSITION_T currQBitPos[M];
   /* syndrome is endowed with cyclic padding in the leading word to avoid 
    * boundary checks. The pad should be at least as long as the bit-len of one 
    * vector of bits which is read during UPC computation */
   DIGIT currSyndrome[NUM_DIGITS_GF2X_ELEMENT+1];
   int check;
   int iteration = 0;

   do {
      gf2x_copy(currSyndrome+1, privateSyndrome);

/*position of the first set bit in the word, counting 63 ... 0*/
#if (MSb_POSITION_IN_MSB_DIGIT_OF_ELEMENT == (DIGIT_SIZE_b-1))
      currSyndrome[0] = currSyndrome[NUM_DIGITS_GF2X_ELEMENT];
#else
       currSyndrome[1] |= (currSyndrome[NUM_DIGITS_GF2X_ELEMENT] << (MSb_POSITION_IN_MSB_DIGIT_OF_ELEMENT+1) );
       currSyndrome[0] = currSyndrome[NUM_DIGITS_GF2X_ELEMENT] >> (DIGIT_SIZE_b- (MSb_POSITION_IN_MSB_DIGIT_OF_ELEMENT+1));
#endif

#define VECTYPE DIGIT
#define VECTYPE_SIZE_b DIGIT_SIZE_b
#define VECTYPE_SIZE_B DIGIT_SIZE_B
#define BIT_SIZE_UPC 8
#define VEC_COMB_MASK 0x0101010101010101ULL

    VECTYPE upcMat[VECTYPE_SIZE_b/BIT_SIZE_UPC];
    VECTYPE packedSynBits;
    for (int i = 0; i < N0; i++) {
       for (int valueIdx = 0; valueIdx < P; valueIdx = valueIdx + DIGIT_SIZE_b) { /* this index will be vectorized*/
          memset(upcMat, 0 , VECTYPE_SIZE_B*(VECTYPE_SIZE_b/BIT_SIZE_UPC));
          /* this fetches DIGIT_SIZE_b bits from each Htrpos, packed, and adds them to the 64 upc counters in upcmat */
          for(int HtrOneIdx = 0; HtrOneIdx < DV; HtrOneIdx++) {
                POSITION_T basePos = (HtrPosOnes[i][HtrOneIdx]+valueIdx);
                basePos = basePos %P ;
                /* lsb here is the one in base pos, others are subseq*/
                packedSynBits = gf2x_get_DIGIT_SIZE_coeff_vector_boundless(currSyndrome,basePos); 
                for(int upcMatRow = 0; upcMatRow < VECTYPE_SIZE_b/BIT_SIZE_UPC; upcMatRow++){
                    upcMat[upcMatRow] += packedSynBits & VEC_COMB_MASK;
                    packedSynBits = packedSynBits >> 1;
                }

              }/* end of for computing 64 upcs*/
         /* commit computed UPCs in the upc vector, in the proper order. 
          * UPCmat essentially needs transposition and linearization by row,
          * starting from the last row */
          for(int upcMatCol = 0; upcMatCol < VECTYPE_SIZE_b/BIT_SIZE_UPC; upcMatCol++){
              VECTYPE upcBuf = 0;
              for(int upcMatRow = 0; upcMatRow < VECTYPE_SIZE_b/BIT_SIZE_UPC; upcMatRow++){
               uint8_t matByte = upcMat[upcMatRow];
#if __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
               upcBuf |=  ((VECTYPE)(matByte)) << (8*upcMatRow);
#else
               upcBuf |=  (upcBuf << 8) + ((VECTYPE)(matByte));
#endif
                  upcMat[upcMatRow] = upcMat[upcMatRow] >> 8;
              }
              VECTYPE* vp = (VECTYPE *)(&unsatParityChecks[i][valueIdx+8*upcMatCol]);
              *(vp) = upcBuf;
          } 
         }
     }

      /* padding unsatisfiedParityChecks */
      for (int i = 0; i < N0; i++) {
          for(int j = 0; j < DIGIT_SIZE_B; j++)
           unsatParityChecks[i][P+j] = unsatParityChecks[i][j];
      }

      /* iteration based threshold determination*/
      int corrt_syndrome_based= thresholds[iteration];

      //Computation of correlation  with a full Q matrix
      for (int i = 0; i < N0; i++) {
            uint8_t qBlockIdxes[M];
            uint16_t qBlockOffsets[M];
            uint8_t qBlockIdxesCursor = 0, linearized_cursor = 0;
            int endQblockIdx = 0;
            for (int blockIdx = 0; blockIdx < N0; blockIdx++) {
              endQblockIdx += qBlockWeights[blockIdx][i];
              for (; linearized_cursor < endQblockIdx; linearized_cursor++) {
                qBlockIdxes[linearized_cursor] = qBlockIdxesCursor;
                qBlockOffsets[linearized_cursor] =  (uint16_t) qBlockIdxesCursor;
              }
              qBlockIdxesCursor++;
            }
         for (int j = 0; j <= P; j+=DIGIT_SIZE_B) {
            int currQoneIdx; // position in the column of QtrPosOnes[][...]
            /* a vector of (single byte wide) correlation values is computed
             * at each iteration of the outer loop. 8x vectorization is possible
             * with uint64_t. The most significant byte in the uint64_t vector 
             * of uint8_t is the one matching the outer loop index position*/
            uint64_t correlation = 0; 

            /* a single correlation value is the sum of m "block circulant positioned" values of UPC */
            /* I need a 8x storage of UPCs which I can load as 8x bytes from upc vector starting from
             * the block circulant position. 
             * A trivial implementation will thus have m variables containing 8 consecutive UPCs.
             * 
             * A more reasonable one uses a single variable which gets filled 
             * with 8x UPCs per shot, bar circulant cases related to 8 consecutive correlation
             * computations, and a single circulant "one" from Q (ideal opts will unroll). 
             * Since the outer loop (j indexed, at line 183) is vectorized, a trailer of the last P%DIGIT_SIZE_b should be added 
             * (or not: a single vector computation discarding part of the results in correlation is equally fine) */

            /* This vectorization strategy will benefit from having a vector with the block circulant pos.s in Q
             * and another one with the block idxes. Both can be vector packed: pos.s are 16 bit wide, block idx.es are 2b wide
             * FTTB I will generate both of them at runtime, but the latter can be generated as a #define.*/

            for (currQoneIdx = 0; currQoneIdx < M; currQoneIdx++) {
                int tmp = QtrPosOnes[i][currQoneIdx]+j; /* this must turn into a vector load of 8x UPCs with consecutive j indices */
                tmp = tmp % P ;  /* this will need to be checked before and have a split vector load */
                currQBitPos[currQoneIdx] = tmp;  /* only a single base index must be kept per vector load*/
                
                uint64_t * vp = (uint64_t*)(&unsatParityChecks[qBlockOffsets[currQoneIdx]][tmp]);
                correlation += *vp;
            }

            /* Correlation based flipping */
            for (int delta_j = 0 ; delta_j < DIGIT_SIZE_B && (j+delta_j < P) ; delta_j++){
#if __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
               uint64_t single_correl = ( correlation >> (8*delta_j) ) & 0xFF; 
#else
               uint64_t single_correl = ( correlation >> (8*(DIGIT_SIZE_B-1-delta_j)) ) & 0xFF; 
#endif
               if ((single_correl > corrt_syndrome_based)) { /* correlation will always be null for garbage bits */
                  gf2x_toggle_coeff(out+NUM_DIGITS_GF2X_ELEMENT*i, (j+delta_j)%P); /*this can be vectorized if a vector of correlations are available*/
                  for (int v = 0; v < M; v++) { /* the same goes here */
                     unsigned syndromePosToFlip;
                     for (int HtrOneIdx = 0; HtrOneIdx < DV; HtrOneIdx++) {
                        syndromePosToFlip = (HtrPosOnes[qBlockIdxes[v]][HtrOneIdx] + (currQBitPos[v]+delta_j)%P );
                        syndromePosToFlip = syndromePosToFlip >= P ? syndromePosToFlip - P : syndromePosToFlip;
                        gf2x_toggle_coeff(privateSyndrome, syndromePosToFlip);
                     }
                  } // end for v
               } // end if
            } // end for flipping correlations exceeding threshold
         } // end for j
      } // end for i

      iteration = iteration + 1;
      check = 0;
      while (check < NUM_DIGITS_GF2X_ELEMENT && privateSyndrome[check] == 0) {
          check++;
      }

   } while (iteration < ITERATIONS_MAX && check < NUM_DIGITS_GF2X_ELEMENT);

   return (check == NUM_DIGITS_GF2X_ELEMENT);
}  // end QdecodeSyndromeThresh_bitFlip_sparse
#endif
