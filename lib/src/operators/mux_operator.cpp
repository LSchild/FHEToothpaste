//
// Created by leonard on 3/23/26.
//

#include "operators/mux_operator.h"
#include "static/gadget_decomp.h"

#include "utils/math_utils.h"
#include "utils/generic_utils.h"


#include <iostream>

MuxOperator::MuxOperator(std::shared_ptr<MathWorker> ntt, uint64_t basis_log2, uint64_t digits) :
    m_ntt(ntt),
    m_basis(1ull << basis_log2),
    m_basis_log2(basis_log2),
    m_digits(digits)
{
    auto deg = ntt->GetDimension();
    // space for 2
    m_scratch_space.resize(deg * (8 * digits + 2));
    m_mask = (1 << m_basis_log2) - 1;
    m_modulus_bits = IntLog2(m_ntt->GetModulus());
    double ceiled_digits = std::ceil(double(m_modulus_bits) / m_basis_log2);
    m_first_shift = m_basis_log2 * uint64_t(ceiled_digits - 1);

}

void MuxOperator::RLWEPrimeProduct(uint64_t *acc, const uint64_t *rlwe_prime, const uint64_t *polynomial, bool add_to_acc) {
    // TODO: Is it better to just zero out the acc if \add_to_acc = false ?
    auto N = m_ntt->GetDimension();
    auto Q = m_ntt->GetModulus();

    uint64_t* scratch_space = m_scratch_space.data();

    auto current_shift = m_first_shift;
    auto loop_start = 0;

    if (add_to_acc) {
        loop_start = 0;
    } else {

        loop_start = 1;

        for (uint64_t x_i = 0; x_i < N; x_i++) {
            scratch_space[x_i] = (polynomial[x_i] >> current_shift) & m_mask;
        }

        m_ntt->ForwardNTT(scratch_space, scratch_space);
        m_ntt->MulModElw(acc, rlwe_prime, scratch_space, N);
        m_ntt->MulModElw(acc + N, rlwe_prime + N, scratch_space, N);

        current_shift -= m_basis_log2;
        rlwe_prime += 2 * N;

    }


    for(uint32_t d_i = loop_start; d_i < m_digits; d_i++) {

        for (uint64_t x_i = 0; x_i < N; x_i++) {
            scratch_space[x_i] = (polynomial[x_i] >> current_shift) & m_mask;
        }



        // Apply NTT to current digit
        m_ntt->ForwardNTT(scratch_space, scratch_space);
        // multiply by A component of RLWE' sample
        m_ntt->MulModElw(scratch_space + N, rlwe_prime, scratch_space, N);
        // Add A * d_i to A component of accumulator
        m_ntt->AddModElw(acc, acc, scratch_space + N, N);
        // multiply digit with B component
        m_ntt->MulModElw(scratch_space, rlwe_prime + N, scratch_space, N);
        // add to B component
        m_ntt->AddModElw(acc + N, scratch_space, acc + N, N);

        current_shift -= m_basis_log2;
        rlwe_prime += 2 * N;


    }
}

std::shared_ptr<MathWorker> MuxOperator::GetNTT() const {
    return m_ntt;
}

uint64_t MuxOperator::GetModulusBits() const {
    return m_modulus_bits;
}

uint64_t MuxOperator::GetRGSWBasis() const {
    return m_basis;
}

uint64_t MuxOperator::GetRGSWBasisLog2() const {
    return m_basis_log2;
}

uint64_t MuxOperator::GetRGSWDigits() const {
    return m_digits;
}

uint64_t MuxOperator::GetMask() const {
    return m_mask;
}

uint64_t MuxOperator::GetFirstShift() const {
    return m_first_shift;
}

void MuxOperator::ExternalProduct(uint64_t * __restrict result, const uint64_t* __restrict rgsw, const uint64_t* __restrict rhs, bool add_to_result) {

    auto N = m_ntt->GetDimension();

    RLWEPrimeProduct(result, rgsw, rhs, add_to_result);
    RLWEPrimeProduct(result, rgsw + 2 * N * m_digits, rhs + N, true);

}

void MuxOperator::BinaryMux(uint64_t *result, const uint64_t *const rgsw_control, uint64_t *const case_0, uint64_t *const case_1) {
    auto N = m_ntt->GetDimension();
    auto Q = m_ntt->GetModulus();

    m_ntt->SubModElw(m_scratch_space.data() + 2 * N, case_1, case_0, 2 * N);

    m_ntt->BackwardNTT(m_scratch_space.data() + 2 * N, m_scratch_space.data() + 2 * N);
    m_ntt->BackwardNTT(m_scratch_space.data() + 3 * N, m_scratch_space.data() + 3 * N);

    ExternalProduct(result, rgsw_control, m_scratch_space.data() + 2 * N);
    m_ntt->AddModElw(result, result, case_0, 2 * N);
}

void MuxOperator::BinaryCMux(uint64_t *acc, const uint64_t *const rgsw_control, const uint64_t *const monomial) {
    // TODO optimize
    auto N = m_ntt->GetDimension();
    auto Q = m_ntt->GetModulus();

    AlignedBuffer tmp(2 * N);


    m_ntt->MulModElw(m_scratch_space.data() + 2 * N, acc, monomial, N);
    m_ntt->MulModElw(m_scratch_space.data() + 3 * N, acc + N, monomial, N);

    m_ntt->BackwardNTT(m_scratch_space.data() + 2 * N, m_scratch_space.data() + 2 * N);
    m_ntt->BackwardNTT(m_scratch_space.data() + 3 * N, m_scratch_space.data() + 3 * N);

    ExternalProduct(tmp.data(), rgsw_control, m_scratch_space.data() + 2 * N);
    m_ntt->AddModElw(acc, acc, tmp.data(), 2 * N);
}



void MuxOperator::TernaryCMux(uint64_t *acc, const uint64_t *const rgsw_controls, const uint64_t *const monomial) {


    // set up digit buffer as pair s.t.
    // db = [ NTT(d_a0) NTT(d_a0) NTT(d_a1) ... NTT(d_b0) NTT(d_b0) ...]
    // X = C0 \times db // 1 call to modmul
    // reduce rows // log2(m_digits) calls (technically)
    // reduce columnsc
    // mul-poly & repeat for C1
    // Next, collapse C1 into it [C0 | C1] * [db | db]
    auto N = m_ntt->GetDimension();
    auto Q = m_ntt->GetModulus();
    auto intt_buffer = m_scratch_space.data();
    auto digit_buffer = m_scratch_space.data() + 2 * N;
    auto mon_pos = monomial ;
    auto mon_neg = monomial + N ;

    // intt_buffer = [a]
    m_ntt->BackwardNTT(intt_buffer, acc);
    SignedDigitDecomposeRep2NTT(digit_buffer, intt_buffer, m_ntt, m_digits, m_basis_log2, m_modulus_bits);
    m_ntt->BackwardNTT(intt_buffer, acc + N);
    SignedDigitDecomposeRep2NTT(digit_buffer + m_digits * 2 * N, intt_buffer, m_ntt,  m_digits, m_basis_log2, m_modulus_bits);



    std::copy(digit_buffer, digit_buffer + 4 * N * m_digits, digit_buffer + 4 * N * m_digits);
    m_ntt->MulModElw(digit_buffer, rgsw_controls, digit_buffer, 8 * N * m_digits);

    // reduce left side, by carrying it forward
    for(uint64_t i = 0; i < 2 * m_digits - 1; i++) {
        m_ntt->AddModElw(digit_buffer + (i + 1) * 2 * N,digit_buffer + i * 2 * N, digit_buffer + (i + 1) * 2 * N, 2 * N);
    }
    m_ntt->MulModElw(intt_buffer,     digit_buffer + (2 * m_digits - 1) * 2 * N,     mon_pos, N);
    m_ntt->MulModElw(intt_buffer + N, digit_buffer + (2 * m_digits - 1) * 2 * N + N, mon_pos, N);
    m_ntt->AddModElw(acc, intt_buffer, acc, 2 * N);


    // reduce right side, by carrying it forward
    digit_buffer += 4 * N * m_digits;
    for(uint64_t i = 0; i < 2 * m_digits - 1; i++) {
        m_ntt->AddModElw(digit_buffer + (i + 1) * 2 * N,digit_buffer + i * 2 * N, digit_buffer + (i + 1) * 2 * N, 2 * N);
    }
    m_ntt->MulModElw(intt_buffer, digit_buffer + (2 * m_digits - 1) * 2 * N, mon_neg, N);
    m_ntt->MulModElw(intt_buffer + N, digit_buffer + (2 * m_digits - 1) * 2 * N + N, mon_neg, N);
    m_ntt->AddModElw(acc, intt_buffer, acc, 2 * N);

}

void MuxOperator::MultiMux(uint64_t *rlwe_sample, uint64_t k, const uint64_t *const rgsw_control_bits, const uint64_t *const weights) {
    const auto N = m_ntt->GetDimension();

    if (k == 2) {
        TernaryCMux(rlwe_sample, rgsw_control_bits, weights);
    }

    const auto Q = m_ntt->GetModulus();
    const auto rlwe_prime_size = 2 * N * m_digits;
    const auto rgsw_size = 2 * rlwe_prime_size;
    const auto number_of_additions_before_reduction = 1ull << (64 - m_modulus_bits);

    const auto scratch_space = m_scratch_space.data();
    const auto digit_buffer = m_scratch_space.data() + 2 * N;

    m_ntt->BackwardNTT(scratch_space, rlwe_sample);
    m_ntt->BackwardNTT(scratch_space + N, rlwe_sample + N);

    // digit decomp
    SignedDigitDecomposeRep2NTT(digit_buffer, scratch_space, m_ntt, m_digits, m_basis_log2, m_modulus_bits);
    SignedDigitDecomposeRep2NTT(digit_buffer + rlwe_prime_size, scratch_space + N, m_ntt, m_digits, m_basis_log2, m_modulus_bits);

    /*
    for(uint32_t jj = 0; jj < 2; jj++) {
        auto current_shift = m_first_shift;
        auto input_block = scratch_space + jj * N;
        auto current_output_block = digit_buffer + jj * rlwe_prime_size;
        for (uint64_t d_i = 0; d_i < m_digits; d_i++) {
            for (uint64_t i = 0; i < N; i++) {
                current_output_block[i] = (input_block[i] >> current_shift) & m_mask;
            }
            m_ntt->ComputeForward(current_output_block, current_output_block, 1, 1);
            std::copy(current_output_block, current_output_block + N, current_output_block + N);
            current_shift -= m_basis_log2;
            current_output_block += 2 * N;
        }
    } */

    const auto acc_target = digit_buffer + rgsw_size;
    auto red_ctr = 0;
    // main loop

    if (m_modulus_bits < 62) {

        for (uint64_t control_idx = 0; control_idx < k; control_idx++) {
            // do the main product
            m_ntt->MulModElw(acc_target, rgsw_control_bits + control_idx * rgsw_size, digit_buffer, rgsw_size);
            // ext-prod output is now given as sum of rows of prior rlwe_sample
            // note, we accumulate into the first row
            for (uint64_t d_i = 1; d_i < 2 * m_digits; d_i++) {
                auto addition_rhs = acc_target + d_i * 2 * N;
                for (uint64_t coef_i = 0; coef_i < 2 * N; coef_i++) {
                    acc_target[coef_i] += addition_rhs[coef_i];
                }
                if (++red_ctr >= number_of_additions_before_reduction) {
                    red_ctr = 0;
                    m_ntt->ReduceModElw(acc_target, acc_target, 2 * N);
                }
            }

            // multiply with monomial/weight
            auto weight_idx = weights + control_idx * N;
            m_ntt->MulModElw(acc_target, acc_target, weight_idx, N);
            m_ntt->MulModElw(acc_target + N, acc_target + N, weight_idx, N);

            // add to rlwe_sample
            m_ntt->AddModElw(rlwe_sample, acc_target, rlwe_sample, 2 * N);
        }

    } else {
        for (uint64_t control_idx = 0; control_idx < k; control_idx++) {
            // do the main product
            m_ntt->MulModElw(acc_target, rgsw_control_bits + control_idx * rgsw_size, digit_buffer, rgsw_size);
            // ext-prod output is now given as sum of rows of prior rlwe_sample
            // note, we accumulate into the first row
            for (uint64_t d_i = 1; d_i < 2 * m_digits; d_i++) {
                auto addition_rhs = acc_target + d_i * 2 * N;
                m_ntt->AddEqModElw(acc_target, addition_rhs, 2 * N);
            }

            // multiply with monomial/weight
            auto weight_idx = weights + control_idx * N;
            m_ntt->MulModElw(acc_target, acc_target, weight_idx, N);
            m_ntt->MulModElw(acc_target + N, acc_target + N, weight_idx, N);

            // add to rlwe_sample
            m_ntt->AddModElw(rlwe_sample, acc_target, rlwe_sample, 2 * N);
        }
    }

}
