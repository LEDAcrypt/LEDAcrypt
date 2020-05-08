#pragma once

#define NUM_TESTS 30

void print_PKC_mceliece_cca2_parameters(int machine_readable);
void test_PKC_mceliece_cca2_code(int isSeedFixed,
                                char* seed,
                                long unsigned int NumTests,
                                int machineReadable);
