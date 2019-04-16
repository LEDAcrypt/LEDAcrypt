/**
 *
 * <gf2x_arith.h>
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

#pragma once

#include "gf2x_limbs.h"
#include "architecture_detect.h"

/*----------------------------------------------------------------------------*/
/*
 * Elements of GF(2)[x] are stored in compact dense binary form.
 *
 * Each bit in a byte is assumed to be the coefficient of a binary
 * polynomial f(x), in Big-Endian format (i.e., reading everything from
 * left to right, the most significant element is met first):
 *
 * byte:(0000 0000) == 0x00 ... f(x) == 0
 * byte:(0000 0001) == 0x01 ... f(x) == 1
 * byte:(0000 0010) == 0x02 ... f(x) == x
 * byte:(0000 0011) == 0x03 ... f(x) == x+1
 * ...                      ... ...
 * byte:(0000 1111) == 0x0F ... f(x) == x^{3}+x^{2}+x+1
 * ...                      ... ...
 * byte:(1111 1111) == 0xFF ... f(x) == x^{7}+x^{6}+x^{5}+x^{4}+x^{3}+x^{2}+x+1
 *
 *
 * A "machine word" (A_i) is considered as a DIGIT.
 * Bytes in a DIGIT are assumed in Big-Endian format:
 * E.g., if sizeof(DIGIT) == 4:
 * A_i: A_{i,3} A_{i,2} A_{i,1} A_{i,0}.
 * A_{i,3} denotes the most significant byte, A_{i,0} the least significant one.
 * f(x) ==   x^{31} + ... + x^{24} +
 *         + x^{23} + ... + x^{16} +
 *         + x^{15} + ... + x^{8}  +
 *         + x^{7}  + ... + x^{0}
 *
 *
 * Multi-precision elements (i.e., with multiple DIGITs) are stored in
 * Big-endian format:
 *           A = A_{n-1} A_{n-2} ... A_1 A_0
 *
 *           position[A_{n-1}] == 0
 *           position[A_{n-2}] == 1
 *           ...
 *           position[A_{1}]  ==  n-2
 *           position[A_{0}]  ==  n-1
 */
/*----------------------------------------------------------------------------*/

#define TC3

#define MIN_KAR_DIGITS 9
#define MIN_TOOM_DIGITS 50

#if defined(TC3) 
#define GF2X_MUL gf2x_mul_TC3
#else
#if defined(HIGH_PERFORMANCE_X86_64) || defined(HIGH_COMPATIBILITY_X86_64)
#define GF2X_MUL gf2x_mul_avx
#else /* if no ISA-specific optimizations are available */
#define GF2X_MUL gf2x_mul_comb
#endif

#endif


/*----------------------------------------------------------------------------*/

static inline void gf2x_add(const int nr, DIGIT Res[],
                            const int na, const DIGIT A[],
                            const int nb, const DIGIT B[]) {
#if (defined HIGH_PERFORMANCE_X86_64)
 __m256i a, b;
 unsigned i;
 for(i = 0; i < nr/4; i++){
     a = _mm256_lddqu_si256( (__m256i *)A + i );
     b = _mm256_lddqu_si256( (__m256i *)B + i );
     _mm256_storeu_si256(((__m256i *)Res + i), _mm256_xor_si256(a, b));
 }
 i = i*2;
 if(nr %4 >= 2){
 __m128i c, d;
     c = _mm_lddqu_si128( (__m128i *)A + i );
     d = _mm_lddqu_si128( (__m128i *)B + i );
     _mm_storeu_si128(((__m128i *)Res + i), _mm_xor_si128(c, d));
     i++;
 }

 if( (nr & 1) == 1){
      Res[nr-1] = A[nr-1] ^ B[nr-1];
 }

#elif (defined HIGH_COMPATIBILITY_X86_64)
 __m128i a, b;
 for (unsigned i = 0; i < nr/2; i++){
     a = _mm_lddqu_si128( (__m128i *)A + i );
     b = _mm_lddqu_si128( (__m128i *)B + i );
     _mm_storeu_si128(((__m128i *)Res + i), _mm_xor_si128(a, b));
 }
 if( (nr & 1) != 0){
      Res[nr-1] = A[nr-1] ^ B[nr-1];
 }
#else
   for (unsigned i = 0; i < nr; i++)
      Res[i] = A[i] ^ B[i];
#endif
} // end gf2x_add

/*----------------------------------------------------------------------------*/

void GF2X_MUL(const int nr, DIGIT Res[],
              const int na, const DIGIT A[],
              const int nb, const DIGIT B[]
             );

/* PRE: MAX ALLOWED ROTATION AMOUNT : DIGIT_SIZE_b */
void right_bit_shift_n(const int length, DIGIT in[], const int amount);

/* PRE: MAX ALLOWED ROTATION AMOUNT : DIGIT_SIZE_b */
void left_bit_shift_n(const int length, DIGIT in[], const int amount);
/*----------------------------------------------------------------------------*/
