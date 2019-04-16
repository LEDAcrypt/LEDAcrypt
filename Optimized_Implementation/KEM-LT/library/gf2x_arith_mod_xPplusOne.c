/**
 *
 * <gf2x_arith_mod_xPplusOne.c>
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


#include "gf2x_arith_mod_xPplusOne.h"
#include "rng.h"
#include <string.h>  // memcpy(...), memset(...)
#include <assert.h>
#include "architecture_detect.h"
#include <stdalign.h>
#include <stdio.h>

/*----------------------------------------------------------------------------*/

/* specialized for nin == 2 * NUM_DIGITS_GF2X_ELEMENT, as it is only used
 * by gf2x_mul */
static inline
void gf2x_mod(DIGIT out[],
              const int nin, const DIGIT in[])
{
  DIGIT aux[NUM_DIGITS_GF2X_ELEMENT+1];
  memcpy(aux, in, (NUM_DIGITS_GF2X_ELEMENT+1)*DIGIT_SIZE_B);
#if MSb_POSITION_IN_MSB_DIGIT_OF_MODULUS != 0
  right_bit_shift_n(NUM_DIGITS_GF2X_ELEMENT+1, aux,
                    MSb_POSITION_IN_MSB_DIGIT_OF_MODULUS);
#endif
  gf2x_add(NUM_DIGITS_GF2X_ELEMENT,out,
           NUM_DIGITS_GF2X_ELEMENT,aux+1,
           NUM_DIGITS_GF2X_ELEMENT,in+NUM_DIGITS_GF2X_ELEMENT);
#if MSb_POSITION_IN_MSB_DIGIT_OF_MODULUS != 0
  out[0] &=  ((DIGIT)1 << MSb_POSITION_IN_MSB_DIGIT_OF_MODULUS) - 1 ;
#endif

} // end gf2x_mod


/*----------------------------------------------------------------------------*/
static
inline
void right_bit_shift(const int length, DIGIT in[])
{
#if (defined HIGH_PERFORMANCE_X86_64)
   int j;
   if(length>=5){
   __m256i v,u,x,v_tmp;
   __m128i one;
   one = _mm_setzero_si128();
   one = _mm_insert_epi64(one, 1,0);

   v = _mm256_lddqu_si256( (__m256i *)&in[(length-1)-3]);
   u = _mm256_setzero_si256();
   u = _mm256_insert_epi64(u, in[(length-1)-4], 0);
   x = v;
   x = _mm256_permute4x64_epi64(x,0x93);
   x = _mm256_blend_epi32 (x,u, 0x03);

   _mm256_storeu_si256( ((__m256i *)&in[(length-1)-3]),
                        _mm256_or_si256( _mm256_srl_epi64(v, one),
                                         _mm256_slli_epi64(x,DIGIT_SIZE_b-1)
                                       )
                      );
   u=_mm256_permute4x64_epi64(u,0x39);
   j=(length-1)-8;
   for(; j >= 0;j = j-4){
      v = _mm256_lddqu_si256( (__m256i *)&in[j]);
      /* shuffle V so that V = [8 7 6 5] -> V = [7 6 5 8] to be blended with
         u[X X X 0]*/
      v_tmp = _mm256_permute4x64_epi64(v,0x39);
      x = _mm256_blend_epi32 (v_tmp,u, 0xC0);
      _mm256_storeu_si256( ((__m256i *)&in[j+1]),
                           _mm256_or_si256( _mm256_srl_epi64(x, one),
                                            _mm256_slli_epi64(v,DIGIT_SIZE_b-1)
                                          )
                      );
      u=v_tmp;
      /* useless, I'll just not blend the non-relevant part
      u= _mm256_and_si256(u,u_hi_mask); */
   }
   j+=4;
   /*here the highest word of u contains the unshifted MSW before head*/
   if(j == 0) {
     u=_mm256_srli_epi64(u, 1);
     in[j] = _mm256_extract_epi64 (u, 3);
   } 
   if (j > 0){
       /* stops GCC loop optimizer from complaining from an UB due to signed 
        * integer underflow */
       unsigned x; 
       x= j;
       for(; x>0;x--){
       in[x] >>= 1;
       in[x] |= (in[x-1] & (DIGIT)0x01) << (DIGIT_SIZE_b-1);
     }
     in[x] >>= 1;
   }
} else {
  for(j=length-1; j > 0; j--){
       in[j] >>= 1;
       in[j] |= (in[j-1] & (DIGIT)0x01) << (DIGIT_SIZE_b-1);
     }
     in[j] >>= 1;
}
#elif (defined HIGH_COMPATIBILITY_X86_64)
#define UNR 3
   int j;
   __m128i a,b,c,d,e,f;

   for (j = length-1; j > UNR*2 ; j=j-(UNR*2)) {

      a = _mm_lddqu_si128( (__m128i *)&in[j-1]);  //load in[j-1] and in[j]
      b = _mm_lddqu_si128( (__m128i *)&in[j-2]);  //load in[j-2] and in[j-1]
      c = _mm_lddqu_si128( (__m128i *)&in[j-3]);  //load in[j-3] and in[j-2]
      d = _mm_lddqu_si128( (__m128i *)&in[j-4]);  //load in[j-4] and in[j-3]
      e = _mm_lddqu_si128( (__m128i *)&in[j-5]);  //load in[j-5] and in[j-4]
      f = _mm_lddqu_si128( (__m128i *)&in[j-6]);  //load in[j-5] and in[j-6]

      a = _mm_srli_epi64(a, 1);
      b = _mm_slli_epi64(b, (DIGIT_SIZE_b-1));
      c = _mm_srli_epi64(c, 1);
      d = _mm_slli_epi64(d, (DIGIT_SIZE_b-1));
      e = _mm_srli_epi64(e, 1);
      f = _mm_slli_epi64(f, (DIGIT_SIZE_b-1));


      _mm_storeu_si128(((__m128i *)&in[j-1]), _mm_or_si128(a, b));
      _mm_storeu_si128(((__m128i *)&in[j-3]), _mm_or_si128(c, d));
      _mm_storeu_si128(((__m128i *)&in[j-5]), _mm_or_si128(e, f));

   }

   for(; j > 0; j--){
     in[j] >>= 1;
     in[j] |= (in[j-1] & (DIGIT)0x01) << (DIGIT_SIZE_b-1);
   }
   in[j] >>= 1;
#else
   int j;
   for (j = length-1; j > 0 ; j--) {
      in[j] >>= 1;
      in[j] |=  (in[j-1] & (DIGIT)0x01) << (DIGIT_SIZE_b-1);
   }
   in[j] >>=1;
#endif
} // end right_bit_shift

/*----------------------------------------------------------------------------*/
/* shifts by whole digits */
void left_DIGIT_shift_n(const int length, DIGIT in[], int amount)
{
   int j;
   for (j = 0; (j + amount) < length; j++) {
      in[j] = in[j+amount];
   }
   for (; j < length; j++) {
      in[j] = (DIGIT)0;
   }
} // end left_bit_shift_n

/*----------------------------------------------------------------------------*/
/* may shift by an arbitrary amount*/
static inline
void left_bit_shift_wide_n(const int length, DIGIT in[], int amount)
{
   left_DIGIT_shift_n(length, in, amount / DIGIT_SIZE_b);
   left_bit_shift_n(length, in, amount % DIGIT_SIZE_b);
} // end left_bit_shift_n

/*----------------------------------------------------------------------------*/

#if (defined(DIGIT_IS_UINT8) || defined(DIGIT_IS_UINT16))
static
uint8_t byte_reverse_with_less32bitDIGIT(uint8_t b)
{
   uint8_t r = b;
   int s = (sizeof(b) << 3) - 1;
   for (b >>= 1; b; b >>= 1) {
      r <<= 1;
      r |= b & 1;
      s--;
   }
   r <<= s;
   return r;
} // end byte_reverse_less32bitDIGIT
#endif

#if defined(DIGIT_IS_UINT32)
static
uint8_t byte_reverse_with_32bitDIGIT(uint8_t b)
{
   b = ( (b * 0x0802LU & 0x22110LU) | (b * 0x8020LU & 0x88440LU)
       ) * 0x10101LU >> 16;
   return b;
} // end byte_reverse_32bitDIGIT
#endif

#if defined(DIGIT_IS_UINT64)
static
uint8_t byte_reverse_with_64bitDIGIT(uint8_t b)
{
   b = (b * 0x0202020202ULL & 0x010884422010ULL) % 1023;
   return b;
} // end byte_reverse_64bitDIGIT
#endif

/*----------------------------------------------------------------------------*/

static
DIGIT reverse_digit(const DIGIT b)
{
   int i;
   union toReverse_t {
      uint8_t inByte[DIGIT_SIZE_B];
      DIGIT digitValue;
   } toReverse;

   toReverse.digitValue = b;
#if defined(DIGIT_IS_UINT64)
   for (i = 0; i < DIGIT_SIZE_B; i++)
      toReverse.inByte[i] = byte_reverse_with_64bitDIGIT(toReverse.inByte[i]);
   return __builtin_bswap64(toReverse.digitValue);
#elif defined(DIGIT_IS_UINT32)
   for (i = 0; i < DIGIT_SIZE_B; i++)
      toReverse.inByte[i] = byte_reverse_with_32bitDIGIT(toReverse.inByte[i]);
   return __builtin_bswap32(toReverse.digitValue);
#elif defined(DIGIT_IS_UINT16)
   for (i = 0; i < DIGIT_SIZE_B; i++)
      toReverse.inByte[i] = byte_reverse_with_less32bitDIGIT(toReverse.inByte[i]);
   reversed = __builtin_bswap16(toReverse.digitValue);
#elif defined(DIGIT_IS_UINT8)
   return byte_reverse_with_less32bitDIGIT(toReverse.inByte[0]);
#else
#error "Missing implementation for reverse_digit(...) \
with this CPU word bitsize !!! "
#endif
   return toReverse.digitValue;
} // end reverse_digit


/*----------------------------------------------------------------------------*/

void gf2x_transpose_in_place(DIGIT A[])
{
   /* it keeps the lsb in the same position and
    * inverts the sequence of the remaining bits
    */

   DIGIT mask = (DIGIT)0x1;
   DIGIT rev1, rev2, a00;
   int i, slack_bits_amount = NUM_DIGITS_GF2X_ELEMENT*DIGIT_SIZE_b - P;

   if (NUM_DIGITS_GF2X_ELEMENT == 1) {
      a00 = A[0] & mask;
      right_bit_shift(1, A);
      rev1 = reverse_digit(A[0]);
#if (NUM_DIGITS_GF2X_MOD_P_ELEMENT*DIGIT_SIZE_b - P)
      rev1 >>= (DIGIT_SIZE_b-(P%DIGIT_SIZE_b));
#endif
      A[0] = (rev1 & (~mask)) | a00;
      return;
   }

   a00 = A[NUM_DIGITS_GF2X_ELEMENT-1] & mask;
   right_bit_shift(NUM_DIGITS_GF2X_ELEMENT, A);

   for (i = NUM_DIGITS_GF2X_ELEMENT-1; i >= (NUM_DIGITS_GF2X_ELEMENT+1)/2; i--) {
      rev1 = reverse_digit(A[i]);
      rev2 = reverse_digit(A[NUM_DIGITS_GF2X_ELEMENT-1-i]);
      A[i] = rev2;
      A[NUM_DIGITS_GF2X_ELEMENT-1-i] = rev1;
   }
   if (NUM_DIGITS_GF2X_ELEMENT % 2 == 1)
      A[NUM_DIGITS_GF2X_ELEMENT/2] = reverse_digit(A[NUM_DIGITS_GF2X_ELEMENT/2]);

   if (slack_bits_amount)
      right_bit_shift_n(NUM_DIGITS_GF2X_ELEMENT, A,slack_bits_amount);
   A[NUM_DIGITS_GF2X_ELEMENT-1] = (A[NUM_DIGITS_GF2X_ELEMENT-1] & (~mask)) | a00;
} // end transpose_in_place

/*----------------------------------------------------------------------------*/

/* computes the degree of a polynomial, returns -1 in case the polynomial */
/* is the null polynomial. Currently used only by inverse. */
static inline
int gf2x_degree(const int length, const DIGIT p[]){
     int i;
     for (i = 0; p[i] == 0 && i < length; i++);
     if (i == length) return -1; // in this case p == 0
     // from now on, It is sure that there is at least one bit of the word p[i] that is asserted
     int degree = (length-i-1)*DIGIT_SIZE_b;
     int posMSb = 0;
     DIGIT v = p[i];
#if (defined HIGH_PERFORMANCE_X86_64)
#include <x86intrin.h>
     posMSb = DIGIT_SIZE_b - (int) __lzcnt64( (unsigned long long) v) - 1;
#else
     /* SW method 1 */
     const DIGIT b[] = {0x2, 0xC, 0xF0, 0xFF00, 0xFFFF0000, 0xFFFFFFFF00000000L};
     const unsigned int S[] = {1, 2, 4, 8, 16, 32};
     for (int i = 5; i >= 0; i--) {
        if (v & b[i]) {
                        v >>= S[i];
                        posMSb |= S[i];
        }
     }
     if (v == 0) posMSb = -1;
#endif

     degree = degree + posMSb;
     return degree;
} /* gf2x_degree */


/*----------------------------------------------------------------------------*/
/* computes poly times digit multiplication as a support for KTT inverse */
/* PRE : nr = na + 1 */

#ifdef HIGH_PERFORMANCE_X86_64
#define GF2X_DIGIT_TIMES_POLY_MUL gf2x_digit_times_poly_mul_avx
static
void gf2x_digit_times_poly_mul_avx(const int nr, 
                                     DIGIT Res[NUM_DIGITS_GF2X_ELEMENT+1],
                               const int na, const DIGIT A[],
                               const DIGIT B){

    __m128i prodRes0,prodRes1,
            accumRes,loopCarriedWord,lowToHighWord,
            wideB,wideA;

    int i;
    wideB=_mm_set_epi64x(0, B);
    loopCarriedWord = _mm_set_epi64x(0,0);

    for (i = na-1; i >= 1 ; i=i-2){
      /*wideA contains [ A[i] A[i-1] ] */
      wideA = _mm_lddqu_si128((__m128i *)&A[i-1]);

      prodRes0 = _mm_clmulepi64_si128(wideA, wideB, 1);
      prodRes1 = _mm_clmulepi64_si128(wideA, wideB, 0);

      accumRes = _mm_xor_si128(loopCarriedWord,prodRes0);
      lowToHighWord = _mm_slli_si128(prodRes1,8);
      accumRes = _mm_xor_si128(accumRes,lowToHighWord);

      accumRes = (__m128i) _mm_shuffle_pd( (__m128d) accumRes, 
                                          (__m128d) accumRes, 1);
      _mm_storeu_si128((__m128i *)(&Res[i]), accumRes);

      loopCarriedWord = _mm_srli_si128(prodRes1,8);
    }
    if (i == 0){ /*skipped last iteration i=0, compensate*/
      prodRes0 = _mm_clmulepi64_si128(_mm_set_epi64x(0, A[0]), wideB, 0);

      accumRes = loopCarriedWord;
      accumRes = _mm_xor_si128(accumRes,prodRes0);
      accumRes = (__m128i) _mm_shuffle_pd( (__m128d) accumRes,
                                           (__m128d) accumRes, 1);
      _mm_storeu_si128((__m128i *)(&Res[0]), accumRes);
    } else { /*i == 1*/
        /*regular exit condition, do nothing*/
    }

}

#else
#define GF2X_DIGIT_TIMES_POLY_MUL gf2x_digit_times_poly_mul

void gf2x_digit_times_poly_mul(const int nr, DIGIT Res[NUM_DIGITS_GF2X_ELEMENT+1],
                               const int na, const DIGIT A[],
                               const DIGIT B){

    DIGIT pres[2];
    Res[nr-1]=0;
    for (int i = (nr-1)-1; i >= 0 ; i--){
       GF2X_MUL(2, pres, 1, &A[i], 1, &B);
       Res[i+1] = Res[i+1] ^ pres[1];
       Res[i] =  pres[0];
    }
}
#endif

/*----------------------------------------------------------------------------
*
* Based on: K. Kobayashi, N. Takagi and K. Takagi, "Fast inversion algorithm in 
* GF(2m) suitable for implementation with a polynomial multiply instruction on 
* GF(2)," in IET Computers & Digital Techniques, vol. 6, no. 3, pp. 180-185, 
* May 2012. doi: 10.1049/iet-cdt.2010.0006
*/

int gf2x_mod_inverse_KTT(DIGIT out[], const DIGIT in[]){  /* in^{-1} mod x^P-1 */

#if NUM_DIGITS_GF2X_MODULUS == NUM_DIGITS_GF2X_ELEMENT
 DIGIT s[NUM_DIGITS_GF2X_ELEMENT+1] = {0},
       r[NUM_DIGITS_GF2X_ELEMENT+1];
 r[0]=0;
 memcpy(r+1,in, NUM_DIGITS_GF2X_ELEMENT*DIGIT_SIZE_B);

 /* S starts set to the modulus */
 s[NUM_DIGITS_GF2X_ELEMENT+1-1] = 1;
 s[0+1] |= ((DIGIT)1) << MSb_POSITION_IN_MSB_DIGIT_OF_MODULUS;

 DIGIT v[2*NUM_DIGITS_GF2X_ELEMENT] = {0}, 
       u[2*NUM_DIGITS_GF2X_ELEMENT] = {0};

 u[2*NUM_DIGITS_GF2X_ELEMENT-1] = (DIGIT) 2; /* x */

 int deg_r = NUM_DIGITS_GF2X_ELEMENT*DIGIT_SIZE_b -1;
 int deg_s = NUM_DIGITS_GF2X_ELEMENT*DIGIT_SIZE_b -1;

 DIGIT c,d;
 DIGIT h00,h01,h10,h11;

 DIGIT hibitmask = ( (DIGIT) 1) << (DIGIT_SIZE_b-1);

 DIGIT r_h00[NUM_DIGITS_GF2X_ELEMENT+2];
 DIGIT s_h01[NUM_DIGITS_GF2X_ELEMENT+2];
 DIGIT r_h10[NUM_DIGITS_GF2X_ELEMENT+2];
 DIGIT s_h11[NUM_DIGITS_GF2X_ELEMENT+2];
 DIGIT u_h00[2*NUM_DIGITS_GF2X_ELEMENT+1];
 DIGIT v_h01[2*NUM_DIGITS_GF2X_ELEMENT+1];
 DIGIT u_h10[2*NUM_DIGITS_GF2X_ELEMENT+1];
 DIGIT v_h11[2*NUM_DIGITS_GF2X_ELEMENT+1];

 while(deg_r > 0){
     c=r[1];
     d=s[1];
     if(c == 0){
        left_DIGIT_shift_n(NUM_DIGITS_GF2X_ELEMENT+1,r,1);
        left_DIGIT_shift_n(2*NUM_DIGITS_GF2X_ELEMENT,u,1);
         deg_r = deg_r - DIGIT_SIZE_b;
     } else {
        /* H = I */
        h00 = 1; h01 = 0;
        h10 = 0; h11 = 1;
        for(int j = 1 ; (j < DIGIT_SIZE_b) && (deg_r > 0) ;j++) {
           if ( (c & hibitmask) == 0){ /* */
               c = c << 1;

               h00 = h00 << 1; 
               h01 = h01 << 1;
               deg_r--;
           } else { /* hibit r[0] set */
               if (deg_r == deg_s){
                 deg_r--;
                 if ( (d & hibitmask) == hibitmask){ /* hibit r[0],s[0] set, deg_r == deg_s */
                    DIGIT temp = c;
                    c = (c^d) << 1; /* (c-d)*x */
                    d = temp;
                    /*mult H*/
                    DIGIT r00;
                    r00 = (h00 << 1) ^ (h10 << 1);
                    DIGIT r01;
                    r01 = (h01 << 1) ^ (h11 << 1);
                    h10 = h00;
                    h11 = h01;
                    h00 = r00;
                    h01 = r01;
                 } else { /* hibit r[0] set, s[0] unset, deg_r == deg_s */
                    DIGIT temp;
                    temp = c;
                    c = d << 1;
                    d = temp;
                    /*mult H*/
                    DIGIT r00;
                    r00 = (h10 << 1);
                    DIGIT r01;
                    r01 = (h11 << 1);
                    h10 = h00; 
                    h11 = h01;
                    h00 = r00;
                    h01 = r01;
                 }
               } else { /* if (deg_r != deg_s) */
                  deg_s--;
                  if ( (d & hibitmask) == hibitmask){ /* hibit r[0],s[0] set, deg_r != deg_s */
                     d = (c^d) << 1; /* (c-d) * x*/
                     /* mult H */
                     h10 = (h00 << 1) ^ (h10 << 1);
                     h11 = (h01 << 1) ^ (h11 << 1);
                  } else { /* hibit r[0] set, s[0] unset, deg_r != deg_s */
                     d = d << 1;
                     /*mul H*/

                     h10 = h10 << 1; 
                     h11 = h11 << 1;
                  }
               } /*(deg_r == deg_s)*/
           } /* if ( (c & ((DIGIT 1) << (DIGIT_SIZE_b-1))) == 0) */
        } /* while */
        /*update r , s */

        GF2X_DIGIT_TIMES_POLY_MUL(NUM_DIGITS_GF2X_ELEMENT+2, r_h00,
                                  NUM_DIGITS_GF2X_ELEMENT+1, r,
                                  h00);
        GF2X_DIGIT_TIMES_POLY_MUL(NUM_DIGITS_GF2X_ELEMENT+2, s_h01,
                                  NUM_DIGITS_GF2X_ELEMENT+1, s,
                                  h01);
        GF2X_DIGIT_TIMES_POLY_MUL(NUM_DIGITS_GF2X_ELEMENT+2, r_h10,
                                  NUM_DIGITS_GF2X_ELEMENT+1, r,
                                  h10);
        GF2X_DIGIT_TIMES_POLY_MUL(NUM_DIGITS_GF2X_ELEMENT+2, s_h11,
                                  NUM_DIGITS_GF2X_ELEMENT+1, s,
                                  h11);

        gf2x_add(NUM_DIGITS_GF2X_ELEMENT+1, r,
                 NUM_DIGITS_GF2X_ELEMENT+1, r_h00+1,
                 NUM_DIGITS_GF2X_ELEMENT+1, s_h01+1);

        gf2x_add(NUM_DIGITS_GF2X_ELEMENT+1, s,
                 NUM_DIGITS_GF2X_ELEMENT+1, r_h10+1,
                 NUM_DIGITS_GF2X_ELEMENT+1, s_h11+1);

        /* *********************** update u, v *************************/
        GF2X_DIGIT_TIMES_POLY_MUL(2*NUM_DIGITS_GF2X_ELEMENT+1, u_h00,
                                  2*NUM_DIGITS_GF2X_ELEMENT, u,
                                  h00);
        GF2X_DIGIT_TIMES_POLY_MUL(2*NUM_DIGITS_GF2X_ELEMENT+1, v_h01,
                                  2*NUM_DIGITS_GF2X_ELEMENT, v,
                                  h01);
        GF2X_DIGIT_TIMES_POLY_MUL(2*NUM_DIGITS_GF2X_ELEMENT+1, u_h10,
                                  2*NUM_DIGITS_GF2X_ELEMENT, u,
                                  h10);
        GF2X_DIGIT_TIMES_POLY_MUL(2*NUM_DIGITS_GF2X_ELEMENT+1, v_h11,
                                  2*NUM_DIGITS_GF2X_ELEMENT, v,
                                  h11);

        gf2x_add(2*NUM_DIGITS_GF2X_ELEMENT, u,
                 2*NUM_DIGITS_GF2X_ELEMENT, u_h00+1,
                 2*NUM_DIGITS_GF2X_ELEMENT, v_h01+1);
        gf2x_add(2*NUM_DIGITS_GF2X_ELEMENT, v,
                 2*NUM_DIGITS_GF2X_ELEMENT, u_h10+1,
                 2*NUM_DIGITS_GF2X_ELEMENT, v_h11+1);
     }
 }
 if (deg_r == 0) {
  memcpy(out,u,NUM_DIGITS_GF2X_ELEMENT*DIGIT_SIZE_B);
 }
 else {
  memcpy(out,v,NUM_DIGITS_GF2X_ELEMENT*DIGIT_SIZE_B);
 }
#else
#error IMPLEMENT MEMCPY INTO A LARGER OPERAND
#endif

 return 0;
}

/*----------------------------------------------------------------------------*/
// From:
// Darrel Hankerson, Julio López Hernandez, and Alfred Menezes.
// "Software Implementation of Elliptic Curve Cryptography over Binary Fields".
// CHES 2000, LNCS 1965. Springer-Verlag Berlin Heidelberg 2000
//
// Algorithm 10. Modified Almost Inverse Algorithm for inversion in GF(2^m)
// Acts as a fallback implementation if AVX2 ISA ext.s are not available

int gf2x_mod_inverse(DIGIT out[], const DIGIT in[])     /* in^{-1} mod x^P-1 */
{
   alignas(32) DIGIT b[NUM_DIGITS_GF2X_MODULUS] = {0},
                     c[NUM_DIGITS_GF2X_MODULUS] = {0},
                     u[NUM_DIGITS_GF2X_MODULUS] = {0},
                     v[NUM_DIGITS_GF2X_MODULUS] = {0};
               DIGIT *ptru = &u[0],
                     *ptrv = &v[0],
                     *ptrb = &b[0],
                     *ptrc = &c[0],
                     *tmp  = NULL;
   int deg_b = 0,
       deg_c = 0,
       deg_u = gf2x_degree(NUM_DIGITS_GF2X_ELEMENT, in),
       deg_v = P,
       t;

   *(ptrb+NUM_DIGITS_GF2X_MODULUS-1) = (DIGIT)0x1;

   memcpy(ptru+(NUM_DIGITS_GF2X_MODULUS-NUM_DIGITS_GF2X_ELEMENT),
          in,
          DIGIT_SIZE_B*NUM_DIGITS_GF2X_ELEMENT);

   DIGIT mask;
   *(ptrv+NUM_DIGITS_GF2X_MODULUS-1) = (DIGIT)0x1;
#if (MSb_POSITION_IN_MSB_DIGIT_OF_MODULUS == 0)
   mask = 0x1;
#else
   mask = (((DIGIT)0x1) << MSb_POSITION_IN_MSB_DIGIT_OF_MODULUS);
#endif
   *ptrv = mask;


   while (1) {
     while ( (*(ptru+NUM_DIGITS_GF2X_MODULUS-1) & (DIGIT)0x1) == (DIGIT)0x0 ) {

          {
            int len_u = (deg_u+DIGIT_SIZE_b)/DIGIT_SIZE_b;
            int posLeadWord = NUM_DIGITS_GF2X_MODULUS - len_u;
            right_bit_shift(len_u, ptru+posLeadWord);
          }
          deg_u = deg_u - 1;
          if ( (*(ptrb+NUM_DIGITS_GF2X_MODULUS-1) & (DIGIT)0x1) != (DIGIT)0x0 ){
            *(ptrb+0) ^= mask;
            *(ptrb+NUM_DIGITS_GF2X_MODULUS-1) ^= (DIGIT)0x1;
            if (deg_b != P) {
                               deg_b = P;
            }
            else {
                   int len_b = (deg_b+DIGIT_SIZE_b)/DIGIT_SIZE_b;
                   int posLeadWord = NUM_DIGITS_GF2X_MODULUS - len_b;
                   deg_b = gf2x_degree(len_b, ptrb+posLeadWord);
            } // end if-else
          } // end if

          {
            int len_b = (deg_b+DIGIT_SIZE_b)/DIGIT_SIZE_b;
            int posLeadWord = NUM_DIGITS_GF2X_MODULUS - len_b;
            right_bit_shift(len_b, ptrb+posLeadWord);
          }
          deg_b = deg_b - 1;
     } // end while
     if (deg_u == 0) {
       memcpy(out,
              ptrb+(NUM_DIGITS_GF2X_MODULUS-NUM_DIGITS_GF2X_ELEMENT),
              DIGIT_SIZE_B*NUM_DIGITS_GF2X_ELEMENT);

       return 1; // success!
     }
     if (deg_u < 0) return 0; // failure !
     int max_deg_uv = deg_u;
     if (deg_u < deg_v) {
       tmp = ptru; ptru = ptrv; ptrv = tmp; t = deg_u; deg_u = deg_v; deg_v = t;
       tmp = ptrb; ptrb = ptrc; ptrc = tmp; t = deg_b; deg_b = deg_c; deg_c = t;
       max_deg_uv = deg_v;
     }

     {
       int max_len_uv = (max_deg_uv+DIGIT_SIZE_b)/DIGIT_SIZE_b;
       int posLeadWord = NUM_DIGITS_GF2X_MODULUS - max_len_uv;
       gf2x_add(max_len_uv, ptru+posLeadWord,
                max_len_uv, ptru+posLeadWord,
                max_len_uv, ptrv+posLeadWord);
     }

     if (deg_u < deg_v) { deg_u = deg_v; }
     if (deg_u == deg_v) {
       int len_u = (deg_u+DIGIT_SIZE_b)/DIGIT_SIZE_b;
       int posLeadWord = NUM_DIGITS_GF2X_MODULUS - len_u;
       deg_u = gf2x_degree(len_u, ptru+posLeadWord);
     }

    int max_deg_bc = deg_b;
    if (deg_b < deg_c) { deg_b = deg_c; max_deg_bc = deg_c; }
     {
       int max_len_bc = (max_deg_bc+DIGIT_SIZE_b)/DIGIT_SIZE_b;
       int posLeadWord = NUM_DIGITS_GF2X_MODULUS - max_len_bc;
       gf2x_add(max_len_bc, ptrb+posLeadWord,
                max_len_bc, ptrb+posLeadWord,
                max_len_bc, ptrc+posLeadWord);
     }
/*
     gf2x_add(NUM_DIGITS_GF2X_MODULUS, ptrb,
              NUM_DIGITS_GF2X_MODULUS, ptrb,
              NUM_DIGITS_GF2X_MODULUS, ptrc);
    */
     if (deg_b == deg_c) {
         int len_b = (deg_b+DIGIT_SIZE_b)/DIGIT_SIZE_b;
         int posLeadWord = NUM_DIGITS_GF2X_MODULUS - len_b;
         deg_b = gf2x_degree(len_b, ptrb+posLeadWord);
     }
   } // end while 1
}  // end gf2x_mod_inverse

/*----------------------------------------------------------------------------*/

void gf2x_mod_mul(DIGIT Res[], const DIGIT A[], const DIGIT B[])
{

   DIGIT aux[2*NUM_DIGITS_GF2X_ELEMENT];
   GF2X_MUL(2*NUM_DIGITS_GF2X_ELEMENT, aux,
                 NUM_DIGITS_GF2X_ELEMENT, A,
                 NUM_DIGITS_GF2X_ELEMENT, B);
   gf2x_mod(Res, 2*NUM_DIGITS_GF2X_ELEMENT, aux);

} // end gf2x_mod_mul

/*----------------------------------------------------------------------------*/

/*PRE: the representation of the sparse coefficients is sorted in increasing
 order of the coefficients themselves */
void gf2x_mod_mul_dense_to_sparse(DIGIT Res[],
                                  const DIGIT dense[],
                                  POSITION_T sparse[],
                                  unsigned int nPos)
{
   DIGIT aux[2*NUM_DIGITS_GF2X_ELEMENT];
   DIGIT resDouble[2*NUM_DIGITS_GF2X_ELEMENT];
   
   memset(aux,0,NUM_DIGITS_GF2X_ELEMENT*DIGIT_SIZE_B);
   memcpy(aux+NUM_DIGITS_GF2X_ELEMENT,dense,
          NUM_DIGITS_GF2X_ELEMENT*DIGIT_SIZE_B);
   memset(resDouble,0,NUM_DIGITS_GF2X_ELEMENT*DIGIT_SIZE_B);
   memcpy(resDouble+NUM_DIGITS_GF2X_ELEMENT,dense,
          NUM_DIGITS_GF2X_ELEMENT*DIGIT_SIZE_B);

   if(sparse[0] != INVALID_POS_VALUE) {
      left_bit_shift_wide_n(2*NUM_DIGITS_GF2X_ELEMENT,resDouble,sparse[0]);
      left_bit_shift_wide_n(2*NUM_DIGITS_GF2X_ELEMENT,aux,sparse[0]);

      for (unsigned int i = 1; i < nPos; i++) {
         if (sparse[i] != INVALID_POS_VALUE) {
            left_bit_shift_wide_n(2*NUM_DIGITS_GF2X_ELEMENT,aux, (sparse[i]-sparse[i-1]) );
            gf2x_add(2*NUM_DIGITS_GF2X_ELEMENT,resDouble,
                     2*NUM_DIGITS_GF2X_ELEMENT,aux,
                     2*NUM_DIGITS_GF2X_ELEMENT,resDouble);
         }
      }
   }
   gf2x_mod(Res, 2*NUM_DIGITS_GF2X_ELEMENT, resDouble);
} // end gf2x_mod_mul

/*----------------------------------------------------------------------------*/


void gf2x_transpose_in_place_sparse(int sizeA, POSITION_T A[])
{

   POSITION_T t;
   int i = 0, j;

   if (A[i] == 0) {
      i = 1;
   }
   j = i;

   for (; i < sizeA && A[i] != INVALID_POS_VALUE; i++) {
      A[i] = P-A[i];
   }

   for (i -= 1; j < i; j++, i--) {
      t = A[j];
      A[j] = A[i];
      A[i] = t;
   }

} // end gf2x_transpose_in_place_sparse

/*----------------------------------------------------------------------------*/

void gf2x_mod_mul_sparse(int sizeR, /*number of ones in the result, 
                                     * max sizeA*sizeB */
                         POSITION_T Res[],
                         int sizeA, /*number of ones in A*/
                         const POSITION_T A[],
                         int sizeB, /*number of ones in B*/
                         const POSITION_T B[])
{
   /* compute all the coefficients, filling invalid positions with P*/
   unsigned lastFilledPos=0;
   for(int i = 0 ; i < sizeA ; i++){
     for(int j = 0 ; j < sizeB ; j++){
          uint32_t prod = ((uint32_t) A[i]) + ((uint32_t) B[j]);
          prod = ( (prod >= P) ? prod - P : prod);
       if ((A[i] != INVALID_POS_VALUE) &&
           (B[j] != INVALID_POS_VALUE)){
            Res[lastFilledPos] = prod;
        } else{
            Res[lastFilledPos] = INVALID_POS_VALUE;
        }
        lastFilledPos++;
     }
   }
   while(lastFilledPos < sizeR){
        Res[lastFilledPos] = INVALID_POS_VALUE;
        lastFilledPos++;
   }
   quicksort(Res, sizeR);
   /* eliminate duplicates */
   POSITION_T lastReadPos = Res[0];
   int duplicateCount;
   int write_idx = 0;
   int read_idx = 0;
   while(read_idx < sizeR  && Res[read_idx] != INVALID_POS_VALUE){
      lastReadPos = Res[read_idx];
      read_idx++;
      duplicateCount=1;
      while( (Res[read_idx] == lastReadPos) && (Res[read_idx] != INVALID_POS_VALUE)){
        read_idx++;
        duplicateCount++;
      }
      if (duplicateCount % 2) {
        Res[write_idx] = lastReadPos;
        write_idx++;
      }
   }
   /* fill remaining cells with INVALID_POS_VALUE */
   for(; write_idx < sizeR; write_idx++) {
      Res[write_idx] = INVALID_POS_VALUE;
   }
} // end gf2x_mod_mul_sparse


/*----------------------------------------------------------------------------*/
/* the implementation is safe even in case A or B alias with the result */
/* PRE: A and B should be sorted and have INVALID_POS_VALUE at the end */
void gf2x_mod_add_sparse(int sizeR,
                         POSITION_T Res[],
                         int sizeA,
                         POSITION_T A[],
                         int sizeB,
                         POSITION_T B[])
{
   POSITION_T tmpRes[sizeR];
   int idxA = 0, idxB = 0, idxR = 0;
   while ( idxA < sizeA  &&
           idxB < sizeB  &&
           A[idxA] != INVALID_POS_VALUE &&
           B[idxB] != INVALID_POS_VALUE ) {
      if (A[idxA] == B[idxB]) {
         idxA++;
         idxB++;
      } else {
         if (A[idxA] < B[idxB]) {
            tmpRes[idxR] = A[idxA];
            idxA++;
         } else {
            tmpRes[idxR] = B[idxB];
            idxB++;
         }
         idxR++;
      }
   }

   while (idxA < sizeA && A[idxA] != INVALID_POS_VALUE) {
      tmpRes[idxR] = A[idxA];
      idxA++;
      idxR++;
   }
   while (idxB < sizeB && B[idxB] != INVALID_POS_VALUE) {
      tmpRes[idxR] = B[idxB];
      idxB++;
      idxR++;
   }
   while (idxR < sizeR) {
      tmpRes[idxR] = INVALID_POS_VALUE;
      idxR++;
   }
   memcpy(Res,tmpRes,sizeof(POSITION_T)*sizeR);
} // end gf2x_mod_add_sparse

/*----------------------------------------------------------------------------*/

/* Return a uniform random value in the range 0..n-1 inclusive,
 * applying a rejection sampling strategy and exploiting as a random source
 * the NIST seedexpander seeded with the proper key.
 * Assumes that the maximum value for the range n is 2^32-1
 */
static
int rand_range(const int n, const int logn, AES_XOF_struct *seed_expander_ctx)
{

   unsigned long required_rnd_bytes = (logn+7)/8;
   unsigned char rnd_char_buffer[4];
   uint32_t rnd_value;
   uint32_t mask = ( (uint32_t)1 << logn) - 1;

   do {
      seedexpander(seed_expander_ctx, rnd_char_buffer, required_rnd_bytes);
      /* obtain an endianness independent representation of the generated random
       bytes into an unsigned integer */
      rnd_value =  ((uint32_t)rnd_char_buffer[3] << 24) +
                   ((uint32_t)rnd_char_buffer[2] << 16) +
                   ((uint32_t)rnd_char_buffer[1] <<  8) +
                   ((uint32_t)rnd_char_buffer[0] <<  0) ;
      rnd_value = mask & rnd_value;
   } while (rnd_value >= n);

   return rnd_value;
} // end rand_range

/*----------------------------------------------------------------------------*/
/* Obtains fresh randomness and seed-expands it until all the required positions
 * for the '1's in the circulant block are obtained */

void rand_circulant_sparse_block(POSITION_T *pos_ones,
                                 const int countOnes,
                                 AES_XOF_struct *seed_expander_ctx)
{

   int duplicated, placedOnes = 0;

   while (placedOnes < countOnes) {
      int p = rand_range(NUM_BITS_GF2X_ELEMENT,
                         BITS_TO_REPRESENT(P),
                         seed_expander_ctx);
      duplicated = 0;
      for (int j = 0; j < placedOnes; j++) if (pos_ones[j] == p) duplicated = 1;
      if (duplicated == 0) {
         pos_ones[placedOnes] = p;
         placedOnes++;
      }
   }
} // rand_circulant_sparse_block

/*----------------------------------------------------------------------------*/

void rand_circulant_blocks_sequence(DIGIT sequence[N0*NUM_DIGITS_GF2X_ELEMENT],
                                    const int countOnes,
                                    AES_XOF_struct *seed_expander_ctx)
{

   int rndPos[countOnes],  duplicated, counter = 0;
   memset(sequence, 0x00, N0*NUM_DIGITS_GF2X_ELEMENT*DIGIT_SIZE_B);


   while (counter < countOnes) {
      int p = rand_range(N0*NUM_BITS_GF2X_ELEMENT,BITS_TO_REPRESENT(P),
                         seed_expander_ctx);
      duplicated = 0;
      for (int j = 0; j < counter; j++) if (rndPos[j] == p) duplicated = 1;
      if (duplicated == 0) {
         rndPos[counter] = p;
         counter++;
      }
   }
   for (int j = 0; j < counter; j++) {
      int polyIndex = rndPos[j] / P;
      int exponent = rndPos[j] % P;
      gf2x_set_coeff( sequence + NUM_DIGITS_GF2X_ELEMENT*polyIndex, exponent,
                      ( (DIGIT) 1));
   }

} // end rand_circulant_blocks_sequence

/*----------------------------------------------------------------------------*/
