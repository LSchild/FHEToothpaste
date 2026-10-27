cmake_minimum_required(VERSION 3.26)
project(toothpaste)

set(CMAKE_CXX_STANDARD 20)
include(FetchContent)

# SETUP OpenFHE

find_package(OpenFHE CONFIG REQUIRED PATHS "${CMAKE_CURRENT_SOURCE_DIR}/dependencies/openfhe" NO_DEFAULT_PATH)
# find_package(OpenFHE CONFIG REQUIRED)

if (OpenFHE_FOUND)
    message(STATUS "FOUND PACKAGE OpenFHE")
    message(STATUS "OpenFHE Version: ${BASE_OPENFHE_VERSION}")
    message(STATUS "OpenFHE installed as shared libraries: ${OpenFHE_SHARED}")
    message(STATUS "OpenFHE include files location: ${OpenFHE_INCLUDE}")
    message(STATUS "OpenFHE lib files location: ${OpenFHE_LIBDIR}")
    message(STATUS "OpenFHE Native Backend size: ${OpenFHE_NATIVE_SIZE}")
else()
    message(FATAL_ERROR "PACKAGE OpenFHE NOT FOUND")
endif ()

#set( CMAKE_CXX_FLAGS "${OpenFHE_CXX_FLAGS} -flto" )
set(CMAKE_CXX_FLAGS "-march=native -O3 -flto")
#set(CMAKE_CXX_FLAGS "-mavx2 -g2")

include_directories( ${OPENMP_INCLUDES} )
include_directories( ${OpenFHE_INCLUDE} )
include_directories( ${OpenFHE_INCLUDE}/third-party/include )
include_directories( ${OpenFHE_INCLUDE}/core )
include_directories( ${OpenFHE_INCLUDE}/pke )
include_directories( ${OpenFHE_INCLUDE}/binfhe )

### add directories for other OpenFHE modules as needed for your project

link_directories( ${OpenFHE_LIBDIR} )
link_directories( ${OPENMP_LIBRARIES} )
if(BUILD_STATIC)
    set( CMAKE_EXE_LINKER_FLAGS "${OpenFHE_EXE_LINKER_FLAGS} -static")
    link_libraries( ${OpenFHE_STATIC_LIBRARIES} )
else()
    set( CMAKE_EXE_LINKER_FLAGS ${OpenFHE_EXE_LINKER_FLAGS} )
    link_libraries( ${OpenFHE_SHARED_LIBRARIES} )
endif()

# SETUP HEXL
set(HEXL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/dependencies/hexl")
include_directories("${HEXL_DIR}/include")

# SETUP GTEST
FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

# SETUP OMP
#find_package(OpenMP)
#if (OPENMP_FOUND)
#    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
#endif()

add_library(toothpaste SHARED
        lib/include/operators/automorphism_evaluation.h
        lib/src/operators/automorphism_evaluation.cpp
        lib/include/backend/aligned_vector.h
        lib/include/backend/solinas_worker.h
        lib/include/backend/solinas_internal.h
        lib/include/operators/cggi_blind_rotation.h
        lib/include/operators/bmmp_blind_rotation.h
        lib/include/operators/mux_operator.h
        lib/src/operators/mux_operator.cpp
        lib/src/operators/cggi_blind_rotation.cpp
        lib/include/base_crypto.h
        lib/src/base_crypto.cpp
        lib/src/operators/bmmp_blind_rotation.cpp
        lib/include/utils/generic_utils.h
        lib/include/operators/endo_glwe_conversion.h
        lib/src/operators/lwe_conversion.cpp
        lib/include/static/gadget_decomp.h
        lib/src/operators/rlwe_conversion.cpp
        lib/include/utils/math_utils.h
        lib/src/utils/math_utils.cpp
        lib/include/operators/lwe_to_rgsw_conversion.h
        lib/include/static/container_types.h
        lib/include/interfaces/operator_context.h
        lib/include/static/enum_ids.h
        lib/include/interfaces/blindrotation_operator.h
        lib/include/interfaces/conversion_operator.h
        lib/include/static/keys.h
        lib/src/static/key.cpp
        lib/include/interfaces/decomposition_operator.h
        lib/include/utils/math_utils.h
        lib/include/operators/trace_evaluation.h
        lib/src/operators/trace_evaluation.cpp
        lib/src/operators/lwe_to_rgsw_conversion.cpp
        lib/include/operators/lwe_to_rlwe_packing.h
        lib/src/operators/lwe_to_rlwe_packing.cpp
        lib/include/operators/sap_decomposition.h
        lib/src/operators/sap_decomposer.cpp
        lib/include/backend/backend.h
        lib/include/static/modulus_switching.h
        lib/src/static/modulus_switching.cpp
        lib/include/static/sample_extraction.h
        lib/include/backend/hexl_worker.h
        lib/src/backend/hexl_worker.cpp
        lib/src/backend/backend.cpp
        lib/include/interfaces/functional_bootstrap.h
        lib/include/backend/definitions.h
        lib/include/operators/negacyclic_function_evaluator.h
        lib/src/operators/negacyclic_function_evaluator.cpp
        lib/include/operators/ks_function_evaluator.h
        lib/src/operators/ks_function_evaluator.cpp
        lib/src/backend/solinas.cpp
        lib/include/backend/trivial_worker.h
        lib/src/backend/trivial_worker.cpp
)

target_include_directories(toothpaste PUBLIC lib/include)
target_link_libraries(toothpaste PUBLIC ${HEXL_DIR}/lib64/libhexl.a)

# Testing
enable_testing()

add_executable(test_library
        tests/scheme_switch_tests.cpp
        tests/mux_operator_tests.cpp
        tests/base_crypto_tests.cpp
        tests/cggi_blind_rotation_test.cpp
        tests/bmmp_blind_rotation_tests.cpp
        tests/glwe_conversion_test.cpp
        tests/noise_tests/blind_rotation_noise.cpp
        tests/automorphism_test.cpp
        tests/trace_tests.cpp
        tests/packing_tests.cpp
        tests/sap_decomposer_tests.cpp
        tests/noise_tests/lwe_to_lwe_noise.cpp
        tests/noise_tests/rlwe_to_rlwe_noise.cpp
        tests/backend/solinas_test.cpp)


target_link_libraries(test_library
        toothpaste
    GTest::gtest_main
)
target_include_directories(test_library PUBLIC lib/include)

include(GoogleTest)
gtest_discover_tests(test_library)


#add_executable(pirouette_encrypt iot/pirouette_encrypt.cpp)
#add_executable(respire_encrypt iot/respire_transciphering_encrypt.cpp)