//
// Created by leonard on 6/18/26.
//

#ifndef LARGE_FUNCTIONS_SOLINAS_2_32_1_M1_1_H
#define LARGE_FUNCTIONS_SOLINAS_2_32_1_M1_1_H

#include "backend/backend.h"
#include "backend/aligned_vector.h"

#include <cstdint>



struct Solinas_2_32_1_M1_1 : public MathWorker {

    static const uint64_t Q = 18446744069414584321ull;
    static const uint64_t K = (1ull << 32) - 1;
    static const uint64_t root_2_pow_20 = 3511170319078647661;

    Solinas_2_32_1_M1_1(uint64_t dimension);

    inline uint64_t ReduceMod(__uint128_t val) {
        uint64_t z_lo = (val << 64) >> 64;
        if (z_lo >= Q)
            z_lo -= Q;
        uint64_t z_hi = val >> 64;
        uint64_t z_hi_lo = (z_hi << 32) >> 32;
        uint64_t z_hi_hi = z_hi >> 32;
        auto res = SubMod(z_lo, z_hi_lo+z_hi_hi);
        return AddMod(res, z_hi_lo << 32);
    }

    void ForwardNTT(uint64_t *output, const uint64_t * input) override;

    void ForwardNTT(std::vector<uint64_t> &result, std::vector<uint64_t> &input) override;

    void BackwardNTT(uint64_t *result, const uint64_t *input) override;

    void BackwardNTT(std::vector<uint64_t> &result, std::vector<uint64_t> &input) override;

    uint64_t AddMod(uint64_t left, uint64_t right) override;

    uint64_t SubMod(uint64_t left, uint64_t right) override;

    uint64_t MulMod(uint64_t left, uint64_t right) override;

    uint64_t FMAMod(uint64_t left, uint64_t right, uint64_t summand) override;

    uint64_t InvMod(uint64_t left) override;

    uint64_t PowMod(uint64_t base, uint64_t exponent) override;

    void AddModElw(uint64_t *result, const uint64_t* __restrict left, const uint64_t* __restrict right, uint64_t len) override;

    void AddEqModElw(uint64_t *left_out, const uint64_t* __restrict right, uint64_t len) override;

    void SubModElw(uint64_t *result, const uint64_t* __restrict left, const uint64_t* __restrict right, uint64_t len) override;

    void SubEqModElw(uint64_t *left_out, const uint64_t* __restrict right, uint64_t len) override;

    void MulModElw(uint64_t *result, const uint64_t* __restrict left, const uint64_t* __restrict right, uint64_t len) override;

    void MulEqModElw(uint64_t *left_out, const uint64_t* __restrict right, uint64_t len) override;

    void FMAModElw(uint64_t *result, const uint64_t* __restrict left, const uint64_t* __restrict right, const uint64_t* __restrict summand, uint64_t len) override;

    void FMAEqModElw(uint64_t *left_out, const uint64_t* __restrict right, const uint64_t* __restrict summand, uint64_t len) override;

    void MulSModElw(uint64_t *result, const uint64_t* __restrict left, uint64_t right, uint64_t len) override;

    void MulSEqModElw(uint64_t *left_out, uint64_t right, uint64_t len) override;

    void ReduceModElw(uint64_t *out, const uint64_t *in, uint64_t len) override;

    [[nodiscard]] uint64_t GetModulus() const override;

    [[nodiscard]] uint64_t GetDimension() const override;

private:

    void Precompute();

    AlignedBuffer m_scratch;
    AlignedBuffer m_twiddle;
    AlignedBuffer m_inverse_twiddle;
    AlignedBuffer m_inverse_final_factor;

    uint64_t m_ntt_dim;

};

#endif //LARGE_FUNCTIONS_SOLINAS_2_32_1_M1_1_H
