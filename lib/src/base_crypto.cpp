//
// Created by leonard on 3/24/26.
//
#include "base_crypto.h"

#include "backend/backend.h"
#include "backend/aligned_vector.h"

#include "openfhe.h"

RLWEEncryptor::RLWEEncryptor(uint64_t modulus, uint32_t ring_dimension, double std) {
    m_std = std;
    m_ntt = SelectWorker(modulus, ring_dimension);
}

RLWEEncryptor::RLWEEncryptor(uint64_t modulus, uint32_t ring_dimension, uint64_t root_of_unity, double std) {
    m_std = std;
    m_ntt = SelectWorker(modulus, ring_dimension); // TODO std::make_shared<intel::hexl::NTT>(ring_dimension, Q, root_of_unity);
}

RLWEEncryptor::RLWEEncryptor(std::shared_ptr<MathWorker> ntt, double std) {
    m_ntt = ntt;
    m_std = std;
}


void RLWEEncryptor::MakeRLWE(uint64_t *result, uint64_t *msg, uint64_t *secret_ntt, bool msg_is_ntt) {
    // TODO: for now keep using the OPENFHE samplers so we don't mess up
    auto N = m_ntt->GetDimension();
    auto Q = m_ntt->GetModulus();

    // TODO: double check
    AlignedBuffer tmp = AlignedBuffer(N);
    lbcrypto::DiscreteUniformGeneratorImpl<NativeVector> sampler_unif(Q);
    lbcrypto::DiscreteGaussianGeneratorImpl<NativeVector> sampler_gauss(m_std);

    // uniformly random = is already in ntt format ;)
    auto A = sampler_unif.GenerateVector(N);
    auto E = sampler_gauss.GenerateVector(N, Q);

    for (uint32_t i = 0; i < m_ntt->GetDimension(); i++) {
        result[i] = A[i].ConvertToInt<uint64_t>();
        tmp[i] = E[i].ConvertToInt<uint64_t>();
    }


    m_ntt->MulModElw(result + N, result, secret_ntt, N);

    if (msg_is_ntt) {
        m_ntt->ForwardNTT(tmp.data(),tmp.data());
        m_ntt->AddModElw(tmp.data(), tmp.data(), msg, N);
    } else {
        m_ntt->AddModElw(tmp.data(), tmp.data(), msg, N);
        m_ntt->ForwardNTT(tmp.data(),tmp.data());
    }

    m_ntt->AddModElw(result + N, result + N, tmp.data(), N);
}

void RLWEEncryptor::MakeRGSW(uint64_t *result, uint64_t *msg, uint64_t *secret_ntt, uint64_t rgsw_basis, uint64_t rgsw_digits, bool msg_is_ntt) {

    auto N = m_ntt->GetDimension();
    auto Q = m_ntt->GetModulus();
    AlignedBuffer tmp = AlignedBuffer(3 * N);
    AlignedBuffer basis = AlignedBuffer(2 * N);
    std::fill(tmp.begin(), tmp.end(), 0);

    uint64_t max_digits = std::ceil(std::log2((long double)Q)/std::log2(rgsw_basis));

    // make tmp = [msg, -s * msg]
    if (msg_is_ntt) {
        m_ntt->MulModElw(tmp.data(), msg, secret_ntt, N);
        m_ntt->SubModElw(tmp.data() + N, tmp.data() + N, tmp.data(), N);
        std::copy(msg, msg + N, tmp.data());
    } else {
        // tmp = [ NTT(msg) | 0 | 0]
        m_ntt->ForwardNTT(tmp.data(), msg);
        // tmp = [ NTT(msg) | 0 | NTT(secret * msg)]
        m_ntt->MulModElw(tmp.data() + 2 * N, tmp.data(), secret_ntt, N);
        // tmp = [ NTT(msg) | -NTT(secret * msg) | NTT(secret * msg)]
        m_ntt->SubModElw(tmp.data() + N, tmp.data() + N, tmp.data() + 2 * N, N);
    }


    for (uint64_t d_i = 0; d_i < max_digits - rgsw_digits; d_i++) {
        m_ntt->MulSEqModElw(tmp.data(), rgsw_basis, 2 * N);
        //intel::hexl::EltwiseFMAMod(tmp.data(), tmp.data(), rgsw_basis, nullptr, 2 * N, Q, 1);
    }

    for (uint64_t d_i = 0; d_i < rgsw_digits; d_i++) {
        uint64_t true_offset = (rgsw_digits - d_i - 1);
        MakeRLWE(result + true_offset * 2 * N, tmp.data() + N, secret_ntt, true);
        MakeRLWE(result + true_offset * 2 * N + 2 * N * rgsw_digits, tmp.data(), secret_ntt, true);
        m_ntt->MulSEqModElw(tmp.data(), rgsw_basis, 2 * N);
    }

}

void RLWEEncryptor::PhaseRLWE(uint64_t *result, uint64_t *rlwe, uint64_t *secret_ntt) {
    auto N = m_ntt->GetDimension();
    auto Q = m_ntt->GetModulus();
    AlignedBuffer tmp = AlignedBuffer(N);

    m_ntt->MulModElw(tmp.data(), rlwe, secret_ntt, N);
    std::copy(rlwe + N, rlwe + 2 * N, result);
    m_ntt->SubModElw(result, result, tmp.data(), N);
    m_ntt->BackwardNTT(result, result);
}

void RLWEEncryptor::PhaseRGSW(uint64_t *result, uint64_t *rgsw, uint64_t *secret_ntt, uint64_t rgsw_digits) {
    auto N = m_ntt->GetDimension();

    for (uint64_t i = 0; i < 2 * rgsw_digits; i++) {
        PhaseRLWE(result + i * 2 * N, rgsw + i * 2 * N, secret_ntt);
    }
}

uint64_t RLWEEncryptor::GetDimension() {
    return m_ntt->GetDimension();
}

uint64_t RLWEEncryptor::GetModulus() {
    return m_ntt->GetModulus();
}

double RLWEEncryptor::GetStd() {
    return m_std;
}

std::shared_ptr<MathWorker> RLWEEncryptor::GetNTT() {
    return m_ntt;
}


LWEEncryptor::LWEEncryptor(uint64_t modulus, uint64_t lwe_dimension, double std) : m_modulus(modulus), m_dimension(lwe_dimension), m_std(std) {}

void LWEEncryptor::MakeLWE(uint64_t *result, uint64_t msg, const uint64_t *const secret) const {

    lbcrypto::DiscreteUniformGeneratorImpl<NativeVector> sampler_unif(m_modulus);
    lbcrypto::DiscreteGaussianGeneratorImpl<NativeVector> sampler_gauss(m_std);
    auto vec = sampler_unif.GenerateVector(m_dimension);

    auto b = sampler_gauss.GenerateInteger(m_modulus).ConvertToInt<uint64_t>();

    auto precon = intel::hexl::MultiplyFactor(1, 64, m_modulus).BarrettFactor();

    for (uint64_t i = 0; i < m_dimension; i++) {
        auto a_i = vec[i].ConvertToInt<uint64_t>();
        auto prod = intel::hexl::MultiplyMod(a_i, secret[i], precon, m_modulus);
        result[i] = a_i;
        b = intel::hexl::AddUIntMod(b, prod, m_modulus);
    }
    result[m_dimension] = intel::hexl::AddUIntMod(b, msg, m_modulus);
}

void LWEEncryptor::PhaseLWE(uint64_t *result, uint64_t *lwe_vec, const uint64_t *const secret) const {
    uint64_t  acc = 0;
    auto precon = intel::hexl::MultiplyFactor(1, 64, m_modulus).BarrettFactor();

    for (uint64_t i = 0; i < m_dimension; i++) {
        auto prod = intel::hexl::MultiplyMod(lwe_vec[i], secret[i], precon, m_modulus);
        acc = intel::hexl::AddUIntMod(acc, prod, m_modulus);
    }

    *result = intel::hexl::SubUIntMod(lwe_vec[m_dimension], acc, m_modulus);
}

uint64_t LWEEncryptor::PhaseLWE(uint64_t *lwe_vec,const uint64_t *const secret) const {
    uint64_t res;
    PhaseLWE(&res, lwe_vec, secret);
    return res;
}

double LWEEncryptor::GetStd() {
    return m_std;
}

uint64_t LWEEncryptor::GetModulus() {
    return m_modulus;
}

uint64_t LWEEncryptor::GetDimension() {
    return m_dimension;
}
