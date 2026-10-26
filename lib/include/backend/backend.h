//
// Created by leonard on 6/18/26.
//

#ifndef LARGE_FUNCTIONS_BACKEND_H
#define LARGE_FUNCTIONS_BACKEND_H

#include <cstdint>
#include <concepts>
#include <vector>
#include <memory>

struct MathWorker {

    virtual void ForwardNTT(uint64_t* result, const uint64_t* input) = 0;

    virtual void ForwardNTT(std::vector<uint64_t>& result, std::vector<uint64_t>& input) = 0;

    virtual void BackwardNTT(uint64_t* result, const uint64_t* input) = 0;

    virtual void BackwardNTT(std::vector<uint64_t>& result, std::vector<uint64_t>& input) = 0;

    virtual uint64_t AddMod(uint64_t left, uint64_t right) = 0;

    virtual uint64_t SubMod(uint64_t left, uint64_t right) = 0;

    virtual uint64_t MulMod(uint64_t left, uint64_t right) = 0;

    virtual uint64_t FMAMod(uint64_t left, uint64_t right, uint64_t summand) = 0;

    virtual uint64_t InvMod(uint64_t left) = 0;

    virtual uint64_t PowMod(uint64_t base, uint64_t exponent) = 0;

    virtual void AddModElw(uint64_t* __restrict result, const uint64_t* __restrict left, const uint64_t* __restrict right, uint64_t len) = 0;

    virtual void AddEqModElw(uint64_t* __restrict left_out, const uint64_t* __restrict right, uint64_t len) = 0;

    virtual void SubModElw(uint64_t* __restrict result, const uint64_t* __restrict left, const uint64_t* __restrict right, uint64_t len) = 0;

    virtual void SubEqModElw(uint64_t* __restrict left_out, const uint64_t* __restrict right, uint64_t len) = 0;

    virtual void MulModElw(uint64_t* __restrict result, const uint64_t* __restrict left, const uint64_t* __restrict right, uint64_t len) = 0;

    virtual void MulEqModElw(uint64_t* __restrict left_out, const uint64_t* __restrict, uint64_t len) = 0;

    virtual void FMAModElw(uint64_t* __restrict result, const uint64_t* __restrict left, const uint64_t* __restrict, const uint64_t* __restrict summand, uint64_t len) = 0;

    virtual void FMAEqModElw(uint64_t* __restrict left_out, const uint64_t* __restrict right, const uint64_t* __restrict summand, uint64_t len) = 0;

    virtual void MulSModElw(uint64_t* __restrict result, const uint64_t* __restrict left, uint64_t right, uint64_t len) = 0;

    virtual void MulSEqModElw(uint64_t* __restrict left_out, uint64_t right, uint64_t len) = 0;

    virtual void ReduceModElw(uint64_t *out, const uint64_t *in, uint64_t len) = 0;

    [[nodiscard]] virtual uint64_t GetModulus() const = 0;

    [[nodiscard]] virtual uint64_t GetDimension() const = 0;

    virtual ~MathWorker() = default;

};

std::shared_ptr<MathWorker> SelectWorker(uint64_t modulus, uint64_t dimension);

#endif //LARGE_FUNCTIONS_BACKEND_H
