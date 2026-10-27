//
// Created by leonard on 10/27/26.
//

#ifndef TOOTHPASTE_TRIVIAL_WORKER_H
#define TOOTHPASTE_TRIVIAL_WORKER_H

#include "backend/backend.h"

struct TrivialWorker :  MathWorker {

    TrivialWorker(uint64_t modulus, uint64_t dimension);

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

    void AddModElw(uint64_t *result, const uint64_t *left, const uint64_t *right, uint64_t len) override;

    void AddEqModElw(uint64_t *left_out, const uint64_t *right, uint64_t len) override;

    void SubModElw(uint64_t *result, const uint64_t* __restrict left, const uint64_t* __restrict right, uint64_t len) override;

    void SubEqModElw(uint64_t *left_out, const uint64_t* __restrictright, uint64_t len) override;

    void MulModElw(uint64_t *result,const uint64_t* __restrict left, const uint64_t* __restrict right, uint64_t len) override;

    void MulEqModElw(uint64_t *left_out, const uint64_t* __restrict right, uint64_t len) override;

    void FMAModElw(uint64_t *result, const uint64_t* __restrict left, const uint64_t* __restrict right, const uint64_t* __restrict summand, uint64_t len) override;

    void FMAEqModElw(uint64_t *left_out, const uint64_t* __restrict right, const uint64_t* __restrict summand, uint64_t len) override;

    void MulSModElw(uint64_t *result, const uint64_t* __restrict left, uint64_t right, uint64_t len) override;

    void MulSEqModElw(uint64_t *left_out, uint64_t right, uint64_t len) override;

    void ReduceModElw(uint64_t *out, const uint64_t *in, uint64_t len) override;

    uint64_t GetModulus() const override;

    uint64_t GetDimension() const override;

    ~TrivialWorker() override = default;

private:

    uint64_t m_modulus;
    uint64_t m_dimension;
    uint64_t m_barrett_factor;

};

#endif //TOOTHPASTE_TRIVIAL_WORKER_H
