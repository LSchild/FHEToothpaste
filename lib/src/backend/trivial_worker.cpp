//
// Created by leonard on 10/27/26.
//

#include "backend/trivial_worker.h"
#include "hexl/hexl.hpp"

TrivialWorker::TrivialWorker(uint64_t modulus, uint64_t dimension) : m_modulus(modulus), m_dimension(dimension) {

    auto factor = intel::hexl::MultiplyFactor(1, 64, modulus);
    m_barrett_factor = factor.BarrettFactor();
}

void TrivialWorker::ForwardNTT(uint64_t * output, const uint64_t * input) {

    throw std::invalid_argument("Trivial Worker should not be used for NTT !!!");

}

void TrivialWorker::ForwardNTT(std::vector<uint64_t> &result, std::vector<uint64_t> &input) {

    throw std::invalid_argument("Trivial Worker should not be used for NTT !!!");

}

void TrivialWorker::BackwardNTT(uint64_t *result, const uint64_t *input) {

    throw std::invalid_argument("Trivial Worker should not be used for NTT !!!");
}

void TrivialWorker::BackwardNTT(std::vector<uint64_t> &result, std::vector<uint64_t> &input) {

    throw std::invalid_argument("Trivial Worker should not be used for NTT !!!");

}

uint64_t TrivialWorker::AddMod(uint64_t left, uint64_t right) {
    return intel::hexl::AddUIntMod(left, right, m_modulus);
}

uint64_t TrivialWorker::SubMod(uint64_t left, uint64_t right) {
    return intel::hexl::SubUIntMod(left, right, m_modulus);
}

uint64_t TrivialWorker::MulMod(uint64_t left, uint64_t right) {

    return intel::hexl::MultiplyMod(left, right, m_barrett_factor, m_modulus);

}

uint64_t TrivialWorker::FMAMod(uint64_t left, uint64_t right, uint64_t summand) {
    auto tmp = intel::hexl::MultiplyMod(left, right, m_barrett_factor,m_modulus);
    return intel::hexl::AddUIntMod(tmp, summand, m_modulus);
}

uint64_t TrivialWorker::InvMod(uint64_t left) {
    return intel::hexl::InverseMod(left, m_modulus);
}

uint64_t TrivialWorker::PowMod(uint64_t base, uint64_t exponent) {
    return intel::hexl::PowMod(base, exponent, m_modulus);
}

void TrivialWorker::AddModElw(uint64_t *result, const uint64_t *left, const uint64_t *right, uint64_t len) {
    intel::hexl::EltwiseAddMod(result, left, right, len, m_modulus);
}

void TrivialWorker::AddEqModElw(uint64_t *left_out, const uint64_t *right, uint64_t len) {
    intel::hexl::EltwiseAddMod(left_out, left_out, right, len, m_modulus);
}

void TrivialWorker::SubModElw(uint64_t *result, const uint64_t* __restrict left, const uint64_t* __restrict right, uint64_t len) {
    intel::hexl::EltwiseSubMod(result, left, right, len, m_modulus);
}

void TrivialWorker::SubEqModElw(uint64_t *left_out, const uint64_t* __restrict right, uint64_t len) {
    intel::hexl::EltwiseSubMod(left_out, left_out, right, len, m_modulus);
}

void TrivialWorker::MulModElw(uint64_t *result, const uint64_t* __restrict left, const uint64_t* __restrict right, uint64_t len) {
    intel::hexl::EltwiseMultMod(result, left, right, len, m_modulus, 1);
}

void TrivialWorker::MulEqModElw(uint64_t *left_out, const uint64_t* __restrict right, uint64_t len) {
    intel::hexl::EltwiseMultMod(left_out, left_out, right, len, m_modulus, 1);
}

void TrivialWorker::FMAModElw(uint64_t *result, const uint64_t* __restrict left, const uint64_t* __restrict right, const uint64_t* __restrict summand, uint64_t len) {
    for(uint64_t i = 0; i < len; i++) {
        result[i] = FMAMod(left[i], right[i], summand[i]);
    }
}

void TrivialWorker::FMAEqModElw(uint64_t *left_out, const uint64_t* __restrict right, const uint64_t* __restrict summand, uint64_t len) {
    FMAModElw(left_out, left_out, right, summand, len);
}

void TrivialWorker::MulSModElw(uint64_t *result, const uint64_t* __restrict left, uint64_t right, uint64_t len) {
    intel::hexl::EltwiseFMAMod(result, left, right, nullptr, len, m_modulus, 1);
}

void TrivialWorker::MulSEqModElw(uint64_t *left_out, uint64_t right, uint64_t len) {
    intel::hexl::EltwiseFMAMod(left_out, left_out, right, nullptr, len, m_modulus, 1);
}

void TrivialWorker::ReduceModElw(uint64_t *out, const uint64_t *in, uint64_t len) {
    intel::hexl::EltwiseReduceMod(out, in, len, m_modulus, m_modulus, 1);
}

uint64_t TrivialWorker::GetModulus() const {
    return m_modulus;
}

uint64_t TrivialWorker::GetDimension() const {
    return m_dimension;
}
