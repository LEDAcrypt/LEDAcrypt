#pragma once

#include "qc_ldpc_parameters.h"
#include "niederreiter.h"


#define NUM_TESTS 30
void print_KEM_parameters(void);
void test_KEM_niederreiter_code(int ac, char *av[], long unsigned int NumTests);

#define THRESH_MAX_NUM_DECRYPT_ERR 50
void test_KEM_niederreiter_FER(int ac, char *av[], int NumTests);
void print_FER(int numTests, publicKeyNiederreiter_t   *pk,
               privateKeyNiederreiter_t *sk);

