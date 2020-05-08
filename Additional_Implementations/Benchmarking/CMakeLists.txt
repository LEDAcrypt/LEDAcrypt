cmake_minimum_required(VERSION 3.7)

project(LEDAcrypt C)

set(CC gcc)
# set(CC clang)

set(CMAKE_C_FLAGS  "${CMAKE_C_FLAGS} -Wall -pedantic -Wuninitialized -march=native -O3 -g3")
set(SANITIZE "")

message("FLAGS" ${CMAKE_C_FLAGS})

set(CMAKE_C_STANDARD 11)

find_library(KECCAK_LIB keccak)
if(NOT KECCAK_LIB)
 set(STANDALONE_KECCAK 1)
endif()

set(REFERENCE 0)
# default compilation picks optimized codebase
if(REFERENCE)
set(BASE_DIR ../../Reference_Implementation)
set(COMMON_HEADERS
        ${BASE_DIR}/Common/include/H_Q_matrices_generation.h
        ${BASE_DIR}/Common/include/djbsort.h
        ${BASE_DIR}/Common/include/gf2x_limbs.h
        ${BASE_DIR}/Common/include/aes256.h
        ${BASE_DIR}/Common/include/gf2x_arith.h
        ${BASE_DIR}/Common/include/rng.h
        ${BASE_DIR}/Common/include/gf2x_arith_mod_xPplusOne.h
        ${BASE_DIR}/Common/include/sha3.h
)
else()
set(BASE_DIR ../../Optimized_Implementation)
set(COMMON_HEADERS
        ${BASE_DIR}/Common/include/architecture_detect.h
        ${BASE_DIR}/Common/include/H_Q_matrices_generation.h
        ${BASE_DIR}/Common/include/djbsort.h
        ${BASE_DIR}/Common/include/gf2x_limbs.h
        ${BASE_DIR}/Common/include/aes256.h
        ${BASE_DIR}/Common/include/gf2x_arith.h
        ${BASE_DIR}/Common/include/rng.h
        ${BASE_DIR}/Common/include/gf2x_arith_mod_xPplusOne.h
        ${BASE_DIR}/Common/include/sha3.h

)
endif()

set(COMMON_SOURCES
        ${BASE_DIR}/Common/library/H_Q_matrices_generation.c
        ${BASE_DIR}/Common/library/djbsort.c
        ${BASE_DIR}/Common/library/gf2x_arith_mod_xPplusOne.c
        ${BASE_DIR}/Common/library/aes256.c
        ${BASE_DIR}/Common/library/gf2x_arith.c
        ${BASE_DIR}/Common/library/rng.c
        ${BASE_DIR}/Common/library/bf_decoding.c
)

if(REFERENCE EQUAL 0)
   list(APPEND COMMON_SOURCES ${BASE_DIR}/Common/library/inverse_exp.c)
   list(APPEND COMMON_HEADERS ${BASE_DIR}/Common/include/inverse_powers_mod_p.h
                              ${BASE_DIR}/Common/include/inverse_perm_tables.h)
endif()


if(STANDALONE_KECCAK)
  set(KECCAK_EXTERNAL_LIB "")
  set(KECCAK_EXTERNAL_ENABLE "")
  list(APPEND COMMON_SOURCES ${BASE_DIR}/Common/library/sha3_fallback.c)
else()
  set(KECCAK_EXTERNAL_LIB keccak)
  set(KECCAK_EXTERNAL_ENABLE "-DSHA_3_LIBKECCAK")
endif()


foreach(cipher KEM-CPA KEM PKC)

   if(${cipher} STREQUAL "KEM-CPA")
     set(SOURCES
          ${COMMON_SOURCES}
          ${BASE_DIR}/${cipher}/library/kem.c
          ${BASE_DIR}/${cipher}/library/niederreiter_encrypt.c
          ${BASE_DIR}/${cipher}/library/niederreiter_decrypt.c
          ${BASE_DIR}/${cipher}/library/niederreiter_keygen.c
     )
     set(HEADERS
          ${COMMON_HEADERS}
          ${BASE_DIR}/${cipher}/include/api.h
          ${BASE_DIR}/${cipher}/include/niederreiter.h
          ${BASE_DIR}/${cipher}/include/niederreiter_keygen.h
          ${BASE_DIR}/${cipher}/include/bf_decoding.h
          ${BASE_DIR}/${cipher}/include/niederreiter_decrypt.h
          ${BASE_DIR}/${cipher}/include/qc_ldpc_parameters.h
          ${BASE_DIR}/${cipher}/include/niederreiter_encrypt.h
     )
     set(NZEROMAX 4)
   elseif(${cipher} STREQUAL "KEM")
     set(SOURCES
          ${COMMON_SOURCES}
          ${BASE_DIR}/${cipher}/library/kem.c
          ${BASE_DIR}/${cipher}/library/niederreiter_encrypt.c
          ${BASE_DIR}/${cipher}/library/dfr_test.c
          ${BASE_DIR}/${cipher}/library/niederreiter_decrypt.c
          ${BASE_DIR}/${cipher}/library/niederreiter_keygen.c
     )
     set(HEADERS
          ${COMMON_HEADERS}
          ${BASE_DIR}/${cipher}/include/api.h
          ${BASE_DIR}/${cipher}/include/niederreiter.h
          ${BASE_DIR}/${cipher}/include/niederreiter_keygen.h
          ${BASE_DIR}/${cipher}/include/bf_decoding.h
          ${BASE_DIR}/${cipher}/include/niederreiter_decrypt.h
          ${BASE_DIR}/${cipher}/include/qc_ldpc_parameters.h
          ${BASE_DIR}/${cipher}/include/dfr_test.h
          ${BASE_DIR}/${cipher}/include/niederreiter_encrypt.h
     )
     set(NZEROMAX 4)
   elseif(${cipher} STREQUAL "PKC")
     set(SOURCES
          ${COMMON_SOURCES}
          ${BASE_DIR}/${cipher}/library/constant_weight_codec.c
          ${BASE_DIR}/${cipher}/library/dfr_test.c
          ${BASE_DIR}/${cipher}/library/encrypt.c
          ${BASE_DIR}/${cipher}/library/mceliece_cca2_encrypt.c
          ${BASE_DIR}/${cipher}/library/mceliece_cca2_decrypt.c
          ${BASE_DIR}/${cipher}/library/mceliece_keygen.c
     )
     set(HEADERS
          ${COMMON_HEADERS}
          ${BASE_DIR}/${cipher}/include/api.h
          ${BASE_DIR}/${cipher}/include/bf_decoding.h
          ${BASE_DIR}/${cipher}/include/constant_weight_codec.h
          ${BASE_DIR}/${cipher}/include/dfr_test.h
          ${BASE_DIR}/${cipher}/include/mceliece.h
          ${BASE_DIR}/${cipher}/include/mceliece_keygen.h
          ${BASE_DIR}/${cipher}/include/mceliece_cca2_encrypt.h
          ${BASE_DIR}/${cipher}/include/mceliece_cca2_decrypt.h
          ${BASE_DIR}/${cipher}/include/qc_ldpc_parameters.h
     )
     set(NZEROMAX 4)
   endif()

   foreach(nzero RANGE 2 ${NZEROMAX})
     foreach(sl RANGE 1 5 2)
       if(${cipher} STREQUAL "KEM-CPA")
           set(TARGET_BINARY_NAME LEDAcrypt_${cipher}_SL_${sl}_N0_${nzero})
           add_executable(${TARGET_BINARY_NAME}
                              ${HEADERS}
                              ${SOURCES}
                              ./library/main_kem_cpa.c
                              ./library/test_kem_cpa.c
                              ./include/test_kem_cpa.h )
           target_include_directories(${TARGET_BINARY_NAME} PRIVATE
                                          ${BASE_DIR}/${cipher}/include
                                          ${BASE_DIR}/Common/include
                                          ./include)
           target_link_libraries(${TARGET_BINARY_NAME} m ${SANITIZE} ${KECCAK_EXTERNAL_LIB})
           set_target_properties(${TARGET_BINARY_NAME} PROPERTIES
                                    RUNTIME_OUTPUT_DIRECTORY ./bin/${cipher})
           set_property(TARGET ${TARGET_BINARY_NAME} APPEND PROPERTY
            COMPILE_FLAGS "-DCATEGORY=${sl} -DN0=${nzero} ${KECCAK_EXTERNAL_ENABLE} ")
       else()
          foreach(dfrlev 64 SL)
             set(TARGET_BINARY_NAME LEDAcrypt_${cipher}_SL_${sl}_N0_${nzero}_DFR_${dfrlev})
             if(${cipher} STREQUAL "PKC")
               add_executable(${TARGET_BINARY_NAME}
                                 ${HEADERS}
                                 ${SOURCES}
                                 ./library/main_pkc.c
                                 ./library/test_pkc.c
                                 ./include/test_pkc.h )
             else()
               add_executable(${TARGET_BINARY_NAME}
                                 ${HEADERS}
                                 ${SOURCES}
                                 ./library/main_kem.c
                                 ./library/test_kem.c
                                 ./include/test_kem.h )
             endif()
             target_include_directories(${TARGET_BINARY_NAME} PRIVATE
                                          ${BASE_DIR}/${cipher}/include
                                          ${BASE_DIR}/Common/include
                                          ./include)
             target_link_libraries(${TARGET_BINARY_NAME} m ${SANITIZE} ${KECCAK_EXTERNAL_LIB})
             set_target_properties(${TARGET_BINARY_NAME} PROPERTIES
                                    RUNTIME_OUTPUT_DIRECTORY ./bin/${cipher})
            if(dfrlev STREQUAL "SL")
              set_property(TARGET ${TARGET_BINARY_NAME} APPEND PROPERTY
              COMPILE_FLAGS "-DCONSTANT_TIME=1 -DN0=${nzero} -DCATEGORY=${sl} -DDFR_SL_LEVEL=1 ${KECCAK_EXTERNAL_ENABLE}")
            else()
              set_property(TARGET ${TARGET_BINARY_NAME} APPEND PROPERTY
              COMPILE_FLAGS "-DCONSTANT_TIME=1 -DN0=${nzero} -DCATEGORY=${sl} -DDFR_SL_LEVEL=0 ${KECCAK_EXTERNAL_ENABLE}")
            endif()
          endforeach(dfrlev)
       endif()
     endforeach(sl)
   endforeach(nzero)
endforeach(cipher)
