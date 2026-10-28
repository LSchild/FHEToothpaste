//
// Created by leonard on 4/24/26.
//

#ifndef LARGE_FUNCTIONS_GADGET_DECOMP_H
#define LARGE_FUNCTIONS_GADGET_DECOMP_H

#include <iostream>
#include <cstdint>
#include <algorithm>

#include "hexl/hexl.hpp"
#include "backend/backend.h"
#include "utils/math_utils.h"

inline std::pair<uint64_t, uint64_t> GetCorrectorForSignedToUnsigned(uint64_t modulus, uint64_t modulus_bits, uint64_t basis_bits, uint64_t digits) {
    uint64_t corrector = 0;//digits * basis_bits == modulus_bits ? 0 : 1ull << (modulus_bits - digits * basis_bits - 1);

    int64_t shift = modulus_bits - basis_bits;
    while (shift > 0) {
        corrector += 1ull << (shift + basis_bits - 1);
        shift -= basis_bits;
    }
    uint64_t least_digit_correction = 1ull <<(shift + basis_bits - 1);
    corrector += least_digit_correction;

    return std::make_pair(corrector % modulus, least_digit_correction);
}

/**
 * Performs an UNSIGNED gadget decomposition and copies/repeats the result $k$ times.
 * Specifically, given a vector $source$ as input, together with a basis L = (1 << $basis_bits$)
 * writes the $digits$ most significant digits of $source$ in basis $L$ to result $k$ times
 * @tparam k Repetition parameter
 * @param result Output buffer, should have size at least n * k * digits
 * @param source Input buffer, expected to have size n
 * @param n Input dimension
 * @param digits Number of most significant digits to keep
 * @param basis_bits Log(2) of the basis/radix to use
 * @param modulus_bits bits of the ring Q
 */
template<uint32_t k>
__attribute_deprecated__ void UnsignedGadgetDecomposeRep(uint64_t *const result, const uint64_t *const source, uint64_t n, uint64_t digits, uint64_t basis_bits, uint64_t modulus_bits) {
    // we assume that basis_bits | modulus_bits OR modulus_bits/basis_bits > digits
    uint64_t mask = (1ull << basis_bits) - 1;
    uint64_t shift = modulus_bits - basis_bits;

    for(uint64_t d_i = 0; d_i < digits; d_i++) {
        uint64_t* result_row = result + d_i * (k + 1) * n;
        for (uint64_t i = 0; i < n; i++) {
            result_row[i] = (source[i] >> shift) & mask;
        }
        shift -= basis_bits;

        // todo: revisit
        if (shift >= 64) {
            // if we enter this branch, an overflow must have happened
            // implying that digits * base_bits > modulus_bits
            auto old_shift = (shift + basis_bits) % 64;
            mask = (1ull << old_shift) - 1;
        }

        for(uint32_t kk = 0; kk < k; kk++) {
            std::copy(result_row, result_row + n, result_row + (kk + 1) * n);
        }
    }
}
#define UnsignedDigitDecompose UnsignedGadgetDecomposeRep<0>
#define UnsignedDigitDecomposeRep2 UnsignedGadgetDecomposeRep<1>


/**
 * Performs a SIGNED gadget decomposition and copies/repeats the result $k$ times.
 * Specifically, given a vector $source$ as input, together with a basis L = (1 << $basis_bits$)
 * writes the $digits$ most significant digits of $source$ in basis $L$ to result $k$ times
 * @tparam k Repetition parameter
 * @param result Output buffer, should have size at least n * k * digits
 * @param source Input buffer, expected to have size n
 * @param n Input dimension
 * @param digits Number of most significant digits to keep
 * @param basis_bits Log(2) of the basis/radix to use
 * @param modulus_bits bits of the ring Q
 */
template<uint32_t k> __attribute_deprecated__ void SignedGadgetDecomposeRep(uint64_t *const result, uint64_t *const source, uint64_t n, uint64_t digits, uint64_t basis_bits, uint64_t modulus, uint64_t modulus_bits) {


    // We aim to obtain a decomposition w.r.t a basis B of a vector of x such that the digits lie in [-B/2,B/2)
    // Then, note that x in [-B/2, B/2) implies x + B/2 in [0, B) i.e. a normal unsigned digit decomp
    // Next, let corr = sum_i B^i B//2 and for x' = x + corr mod Q let its digits be [x0',x1',...]
    // Then, we obtain the signed digits by setting xi = xi' - B//2 which is correct as
    // \sum B^i xi = \sum B^i (xi' - B//2) = x + corr - \sum B^i B//2 = x

    // we assume that basis_bits | modulus_bits OR modulus_bits/basis_bits > digits
    auto [corrector, least_digit_cor] = GetCorrectorForSignedToUnsigned(modulus, modulus_bits, basis_bits, digits);
    intel::hexl::EltwiseAddMod(source, source, corrector, n, modulus);
    UnsignedGadgetDecomposeRep<k>(result, source, n, digits, basis_bits, modulus_bits);
    if (least_digit_cor == (1ull << (basis_bits - 1))) {
        intel::hexl::EltwiseSubMod(result, result, 1ull << (basis_bits - 1), n * k * digits, modulus);
    }
    else {
        intel::hexl::EltwiseSubMod(result, result, 1ull << (basis_bits - 1), n * k * (digits - 1), modulus);
        intel::hexl::EltwiseSubMod(result + n * k * (digits - 1), result + n * k * (digits - 1), least_digit_cor, n * k, modulus);
    }
}

#define SignedDigitDecompose SignedGadgetDecomposeRep<0>
#define SignedDigitDecomposeRep2 SignedGadgetDecomposeRep<1>

/**
 * Performs an UNSIGNED gadget decomposition with the output in NTT format, and copies/repeats the result $k$ times.
 * Specifically, given a vector $source$ as input, together with a basis L = (1 << $basis_bits$)
 * writes the $digits$ most significant digits of $source$ in basis $L$ to result $k$ times
 * @tparam k Repetition parameter
 * @param result Result buffer of size at least n * k * digits
 * @param source Source buffer assumed of size n
 * @param ntt Pointer to ntt engine (subject to change)
 * @param digits Number of most significant digits to keep
 * @param basis_bits Log(2) of the basis/radix to use
 * @param modulus_bits bits of the ring Q
 */
template<uint32_t k>
void UnsignedGadgetDecomposeRepNTT(uint64_t *const result, const uint64_t *const source, std::shared_ptr<MathWorker> ntt, uint64_t digits, uint64_t basis_bits, uint64_t modulus_bits) {
    // we assume that basis_bits | modulus_bits OR modulus_bits/basis_bits > digits
    const uint64_t mask = (1ull << basis_bits) - 1;
    const uint64_t n = ntt->GetDimension();
    uint64_t shift = modulus_bits - basis_bits;

    for(uint64_t d_i = 0; d_i < digits; d_i++) {
        uint64_t* result_row = result + d_i * (k + 1) * n;
        for (uint64_t i = 0; i < n; i++) {
            result_row[i] = (source[i] >> shift) & mask;
        }
        shift -= basis_bits;
        ntt->ForwardNTT(result_row, result_row);
        for(uint32_t kk = 0; kk < k; kk++) {
            std::copy(result_row, result_row + n, result_row + (kk + 1) * n);
        }
    }
}

#define UnsignedDigitDecomposeNTT UnsignedGadgetDecomposeRepNTT<0>
#define UnsignedDigitDecomposeRep2NTT UnsignedGadgetDecomposeRepNTT<1>

/**
 * Performs an SIGNED gadget decomposition with the output in NTT format, and copies/repeats the result $k$ times.
 * Specifically, given a vector $source$ as input, together with a basis L = (1 << $basis_bits$)
 * writes the $digits$ most significant digits of $source$ in basis $L$ to result $k$ times
 * @tparam k Repetition parameter
 * @param result Result buffer of size at least n * k * digits
 * @param source Source buffer assumed of size n
 * @param ntt Pointer to ntt engine (subject to change)
 * @param digits Number of most significant digits to keep
 * @param basis_bits Log(2) of the basis/radix to use
 * @param modulus_bits bits of the ring Q
 */
template<uint32_t k>
void SignedGadgetDecomposeRepNTT(uint64_t *const result, const uint64_t *const source, std::shared_ptr<MathWorker> ntt,  uint64_t digits, uint64_t basis_bits,uint64_t modulus_bits) {
    // We aim to obtain a decomposition w.r.t a basis B of a vector of x such that the digits lie in [-B/2,B/2)
    // Then, note that x in [-B/2, B/2) implies x + B/2 in [0, B) i.e. a normal unsigned digit decomp
    // Next, let corr = sum_i B^i B//2 and for x' = x + corr mod Q let its digits be [x0',x1',...]
    // Then, we obtain the signed digits by setting xi = xi' - B//2 which is correct as
    // \sum B^i xi = \sum B^i (xi' - B//2) = x + corr - \sum B^i B//2 = x

    // we assume that basis_bits | modulus_bits OR modulus_bits/basis_bits > digits
    const uint64_t modulus = ntt->GetModulus();
    const uint64_t n = ntt->GetDimension();
    const auto last_result_row = result + (digits - 1) * (k + 1) * n;

    auto [corrector, least_digit_correction] = GetCorrectorForSignedToUnsigned(modulus, modulus_bits, basis_bits, digits);
    // Should be fine though
    for(uint64_t i = 0; i < n; i++) {
        last_result_row[i] = ntt->AddMod(source[i], corrector);
    }

    uint64_t mask = (1ull << basis_bits) - 1;
    int64_t shift = modulus_bits - basis_bits;
    uint64_t max_digits = (modulus_bits / basis_bits) + (modulus_bits % basis_bits == 0 ? 0 : 1);

    for(uint64_t d_i = 0; d_i < digits; d_i++) {
        uint64_t* result_row = result + d_i * (k + 1) * n;

        auto delta = 1ull << (basis_bits - 1);

        if (d_i == max_digits - 1)  {
            delta = least_digit_correction;
            shift = 0;
            mask = (delta << 1) - 1;
        }

        for (uint64_t i = 0; i < n; i++) {
            result_row[i] = (last_result_row[i] >> shift) & mask;
        }
        shift -= basis_bits;

        //intel::hexl::EltwiseSubMod(result_row, result_row, delta, n, modulus);
        for(uint64_t i = 0; i < n; i++) {
            result_row[i] = ntt->SubMod(result_row[i], delta);
        }

        ntt->ForwardNTT(result_row, result_row);
        for(uint32_t kk = 0; kk < k; kk++) {
            std::copy(result_row, result_row + n, result_row + (kk + 1) * n);
        }
    }
}

#define SignedDigitDecomposeNTT SignedGadgetDecomposeRepNTT<0>
#define SignedDigitDecomposeRep2NTT SignedGadgetDecomposeRepNTT<1>

#endif //LARGE_FUNCTIONS_GADGET_DECOMP_H
