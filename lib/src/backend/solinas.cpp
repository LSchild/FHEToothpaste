//
// Created by leonard on 8/8/26.
//

#include "backend/solinas_2_32_1_m1_1.h"
#include "utils/generic_utils.h"
#include "utils/math_utils.h"


Solinas_2_32_1_M1_1::Solinas_2_32_1_M1_1(uint64_t dimension) : m_ntt_dim(dimension) {
    Precompute();
}

uint64_t Solinas_2_32_1_M1_1::AddMod(uint64_t left, uint64_t right) {
    auto a = left + K;
    auto t = a + right;

    if (t < a) {
        return t;
    } else {
        return t - K;
    }
}

uint64_t Solinas_2_32_1_M1_1::SubMod(uint64_t left, uint64_t right) {
    auto z = left - right;
    if (left < right) {
        return z - K;
    } else {
        return z;
    }
}

uint64_t Solinas_2_32_1_M1_1::MulMod(uint64_t left, uint64_t right) {
    __uint128_t prod = __uint128_t(left) * right;
    return ReduceMod(prod);
}

uint64_t Solinas_2_32_1_M1_1::FMAMod(uint64_t left, uint64_t right, uint64_t summand) {
    __uint128_t prod = __uint128_t(left) * right;
    prod += summand;

    return ReduceMod(prod);
}

uint64_t Solinas_2_32_1_M1_1::PowMod(uint64_t base, uint64_t exponent) {
    if (exponent >= Q - 1)
        exponent -= (Q - 1);

    if (exponent == 0)
        return 1;
    if (exponent == 1)
        return base;

    uint64_t res = 1;
    uint64_t acc = base;

    while (exponent >= 1) {
        if (exponent & 1)
            res = MulMod(res, acc);
        acc = MulMod(acc, acc);
        exponent >>= 1;
    }

    return res;
}

uint64_t Solinas_2_32_1_M1_1::InvMod(uint64_t left) {
    return PowMod(left, Q - 2);
}


void Solinas_2_32_1_M1_1::AddModElw(uint64_t *result, const uint64_t *left, const uint64_t *right, uint64_t len) {

    for(uint64_t i = 0; i < len; i++) {
        result[i] = AddMod(left[i], right[i]);
    }

}

void Solinas_2_32_1_M1_1::AddEqModElw(uint64_t *result, const uint64_t *right, uint64_t len) {
    for(uint64_t i = 0; i < len; i++) {
        result[i] = AddMod(result[i], right[i]);
    }
}

void Solinas_2_32_1_M1_1::SubModElw(uint64_t *result, const uint64_t* __restrict  left, const uint64_t* __restrict right, uint64_t len) {
    for(uint64_t i = 0; i < len; i++) {
        result[i] = SubMod(left[i], right[i]);
    }
}

void Solinas_2_32_1_M1_1::SubEqModElw(uint64_t *left_out, const uint64_t* __restrict right, uint64_t len) {
    for(uint64_t i = 0; i < len; i++) {
        left_out[i] = SubMod(left_out[i], right[i]);
    }
}

void Solinas_2_32_1_M1_1::MulModElw(uint64_t *result, const uint64_t* __restrict left, const uint64_t* __restrict right, uint64_t len) {
    for(uint64_t i = 0; i < len; i++) {
        result[i] = MulMod(left[i], right[i]);
    }
}

void Solinas_2_32_1_M1_1::MulEqModElw(uint64_t *left_out, const uint64_t* __restrict right, uint64_t len) {
    for(uint64_t i = 0; i < len; i++) {
        left_out[i] = MulMod(left_out[i], right[i]);
    }
}

uint64_t Solinas_2_32_1_M1_1::GetModulus() const {
    return Q;
}

uint64_t Solinas_2_32_1_M1_1::GetDimension() const {
    return m_ntt_dim;
}

void
Solinas_2_32_1_M1_1::FMAModElw(uint64_t *result, const uint64_t* __restrict left, const uint64_t* __restrict right, const uint64_t* __restrict summand, uint64_t len) {
    for(uint64_t i = 0; i < len; i++) {
        result[i] = FMAMod(left[i], right[i], summand[i]);
    }
}

void Solinas_2_32_1_M1_1::FMAEqModElw(uint64_t *left_out, const uint64_t* __restrict right, const uint64_t* __restrict summand, uint64_t len) {
    for(uint64_t i = 0; i < len; i++) {
        left_out[i] = FMAMod(left_out[i], right[i], summand[i]);
    }
}

void Solinas_2_32_1_M1_1::MulSModElw(uint64_t *result, const uint64_t* __restrict left, uint64_t right, uint64_t len) {

    if (right == 0) {
        ZERO_UINT64_ARR(result, len);
        return;
    }

    if (right == 1) {
        std::memcpy(result, left, sizeof(uint64_t) * len);
        return;
    }

    if (right == 2) {
        AddModElw(result, left, left, len);
        return;
    }

    for(uint64_t i = 0; i < len; i++) {
        result[i] = MulMod(left[i], right);
    }
}

void Solinas_2_32_1_M1_1::MulSEqModElw(uint64_t *left_out, uint64_t right, uint64_t len) {
    for(uint64_t i = 0; i < len; i++) {
        left_out[i] = MulMod(left_out[i], right);
    }
}

void Solinas_2_32_1_M1_1::Precompute() {

    m_scratch = AlignedBuffer(m_ntt_dim);

    m_twiddle = AlignedBuffer(m_ntt_dim * 2);
    m_inverse_twiddle = AlignedBuffer(m_ntt_dim * 2);
    m_inverse_final_factor = AlignedBuffer(m_ntt_dim);

    auto log2_dim = IntLog2(m_ntt_dim);
    // we have a 2^{20}-th root of unity so z^{2^{20}} = 1
    // hence (z^{2^20 / 2N})^{2N} = 1
    auto root = PowMod(Solinas_2_32_1_M1_1::root_2_pow_20, 1ull << (19 - log2_dim));
    auto root_inverse = InvMod(root);

    m_twiddle[0] = 1;
    m_twiddle[1] = root;

    m_inverse_twiddle[0] = 1;
    m_inverse_twiddle[1] = root_inverse;

    for(uint64_t i = 2; i < m_ntt_dim; i++) {
        m_twiddle[i] = MulMod(m_twiddle[i - 1], root);
        m_inverse_twiddle[i] = MulMod(m_inverse_twiddle[i - 1], root_inverse);
    }

    auto dim_inverse = InvMod(m_ntt_dim);
    MulSModElw(m_inverse_final_factor.data(), m_inverse_twiddle.data(), dim_inverse, m_ntt_dim);
}

void Solinas_2_32_1_M1_1::ForwardNTT(std::vector<uint64_t> &result, std::vector<uint64_t> &input) {
    ForwardNTT(result.data(), input.data());
}

void Solinas_2_32_1_M1_1::ForwardNTT(uint64_t *output, const uint64_t *input) {

    const auto dim = m_ntt_dim;
   // bool copy_flag = output == input;

   //ZERO_UINT64_ARR(output, dim);

    auto log2dim = IntLog2(dim);

    for(uint64_t i =0; i < dim; i++) {
        auto idx_rev = ReverseBitsGeneric(i) >> (64 - log2dim);
        if (idx_rev >= i) {
            auto val_a = input[idx_rev];
            auto val_b = input[i];
            output[i] = MulMod(val_a, m_twiddle[idx_rev]);
            output[idx_rev] = MulMod(val_b, m_twiddle[i]);
        }
    }

    for(uint64_t l_m = 1; l_m <= log2dim; l_m++) {

        auto m = 1u << l_m;
        auto m_half = 1u << (l_m - 1);

        for(uint64_t k = 0; k < dim; k += m) {
            auto u = output[k];
            auto t = output[k + m_half];
            output[k] = AddMod(u, t);
            output[k + m_half] = SubMod(u, t);
        }

        for(uint64_t j = 1; j < m_half; j++) {
            auto omega = m_twiddle[(j * dim) >> (l_m - 1)];
            for (uint64_t k = 0; k < dim; k += m) {
                auto u = output[k + j];
                auto t = output[k + j + m_half];
                auto t_omega = MulMod(t, omega);
                output[k + j] = AddMod(u, t_omega);
                output[k + j + m_half] = SubMod(u, t_omega);
            }
        }
    }
}

void Solinas_2_32_1_M1_1::BackwardNTT(std::vector<uint64_t> &result, std::vector<uint64_t> &input) {
    BackwardNTT(result.data(), input.data());
}

void Solinas_2_32_1_M1_1::BackwardNTT(uint64_t *output, const uint64_t *input) {

    const auto dim = m_ntt_dim;

    //ZERO_UINT64_ARR(output, m_ntt_dim);
    //std::memcpy(output, input, sizeof(uint64_t) * m_ntt_dim);

    auto log2dim = IntLog2(dim);

    for(uint64_t i =0; i < dim; i++) {
        auto idx_rev = ReverseBitsGeneric(i) >> (64 - log2dim);
        if (idx_rev >= i) {
            auto val_a = input[idx_rev];
            auto val_b = input[i];
            output[i] = val_a;
            output[idx_rev] = val_b;
        }
    }

    for(uint64_t l_m = 1; l_m <= log2dim; l_m++) {

        auto m = 1u << l_m;
        auto m_half = 1u << (l_m - 1);

        for(uint64_t k = 0; k < dim; k += m) {
            auto u = output[k];
            auto t = output[k + m_half];
            output[k] = AddMod(u, t);
            output[k + m_half] = SubMod(u, t);
        }
        for(uint64_t j = 1; j < m_half; j++) {
            auto omega = m_inverse_twiddle[(j * dim) >> (l_m - 1)];

            for (uint64_t k = 0; k < dim; k += m) {
                auto u = output[k + j];
                auto t = output[k + j + m_half];
                auto t_omega = MulMod(t, omega);

                output[k + j] = AddMod(u, t_omega);
                output[k + j + m_half] = SubMod(u, t_omega);
            }
        }
    }

    MulEqModElw(output, m_inverse_final_factor.data(), dim);
}

void Solinas_2_32_1_M1_1::ReduceModElw(uint64_t *out, const uint64_t *in, uint64_t len) {
    for(uint64_t i =0; i < len; i++) {
        out[i] = ReduceMod(in[i]);
    }
}