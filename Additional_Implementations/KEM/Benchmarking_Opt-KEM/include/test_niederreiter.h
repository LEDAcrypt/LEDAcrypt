#pragma once

#include "qc_ldpc_parameters.h"
#include "niederreiter.h"


#define NUM_TESTS 30
void print_KEM_parameters(void);
void test_KEM_niederreiter_code(int ac, char *av[], long unsigned int NumTests);

