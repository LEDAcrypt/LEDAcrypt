#pragma once

#define NUM_TESTS 30

void print_KEM_parameters(int machine_readable);
void test_KEM_niederreiter_code(int isSeedFixed,
                                char* seed,
                                long unsigned int NumTests,
                                int machineReadable);
