//
// Created by leonard on 3/24/26.
//
#include "operators/cggi_blind_rotation.h"
#include "base_crypto.h"
#include "utils/generic_utils.h"
#include "utils/math_utils.h"
#include "operators/mux_operator.h"

#include <cassert>
#include <iostream>
#include <utility>

CGGIBlindRotationContext::CGGIBlindRotationContext(KeyDistribution distribution, uint64_t modulus, uint64_t ring_dim, uint64_t lwe_dim, uint64_t basis, uint64_t digits, double std) {
    SetKeyDistribution(distribution);
    SetModulus(modulus);
    SetRingDimension(ring_dim);
    SetLWEDimension(lwe_dim);
    SetStd(std);
    SetBlindRotationBasis(basis);
    SetBlindRotationRGSWDigits(digits);

    auto ntt = SelectWorker(modulus, ring_dim);
    SetNTT(ntt);
}

CGGIBlindRotationContext::CGGIBlindRotationContext(KeyDistribution distr, std::shared_ptr<MathWorker> ntt,
                                                   uint64_t lwe_dim, uint64_t basis, uint64_t digits, double std) {
    SetKeyDistribution(distr);
    SetModulus(ntt->GetModulus());
    SetRingDimension(ntt->GetDimension());
    SetLWEDimension(lwe_dim);
    SetStd(std);
    SetBlindRotationBasis(basis);
    SetBlindRotationRGSWDigits(digits);
    SetNTT(ntt);

}

CGGIBlindRotationContext::CGGIBlindRotationContext(const CGGIBlindRotationContext &other)  : enable_shared_from_this(other) {
    SetKeyDistribution(other.m_distribution);
    SetModulus(other.m_modulus);
    SetRingDimension(other.m_N);
    SetStd(other.m_std);
    SetLWEDimension(other.m_n);
    SetBlindRotationBasis(other.m_basis);
    SetBlindRotationRGSWDigits(other.m_digits);
    SetNTT(other.GetNTT());
}

void CGGIBlindRotationContext::SetKeyDistribution(KeyDistribution dist) {
    m_distribution = dist;
}


void CGGIBlindRotationContext::SetModulus(uint64_t modulus) {
    m_modulus = modulus;
}

void CGGIBlindRotationContext::SetStd(double std) {
    m_std = std;
}

void CGGIBlindRotationContext::SetRingDimension(uint64_t ring_dim) {
    m_N = ring_dim;
}

void CGGIBlindRotationContext::SetNTT(std::shared_ptr<MathWorker> ntt) {
    m_ntt = std::move(ntt);
}

void CGGIBlindRotationContext::SetLWEDimension(uint64_t lwe_dim) {
    m_n = lwe_dim;
}

void CGGIBlindRotationContext::SetBlindRotationBasis(uint64_t L) {
    m_basis = L;
    m_basis_log2 = IntLog2(L);
}

void CGGIBlindRotationContext::SetBlindRotationRGSWDigits(uint64_t digits) {
    m_digits = digits;
}


long double CGGIBlindRotationContext::ComputeOutputVariance(long double input_variance) const {

    long double delta = std::pow(2.0, IntLog2(m_modulus)) / std::pow(m_basis, m_digits);
    if (delta <= 1) {
        delta = 0.0;
    }
    delta *= delta / 12.0;

    // assume worst case
    long double message_norm_square = 0.0;
    switch (m_distribution) {
        case BINARY:
            message_norm_square = m_N;
            break;
        case TERNARY:
            message_norm_square = m_N;
            break;
        case GAUSSIAN:
            std::cerr << "Gaussian keys are not supported yet. Setting key element norm to 5" << std::endl;
            message_norm_square = (5 * 5) * m_N;
            break;
    }

    long double fresh_var = m_std * m_std;
    long double rlwe_p_prod_var = m_digits * m_N;
    rlwe_p_prod_var *= m_basis * m_basis;
    rlwe_p_prod_var *= fresh_var;
    rlwe_p_prod_var /= 12.0;

    rlwe_p_prod_var += delta * message_norm_square;

    // update performs RGSW(s_i) * (X^a - 1) * acc + acc
    // -> var((X^a - 1) * acc) = 2 * var(acc)
    // -> var(GSW(s_i) * (X^a - 1) * acc + acc) = 2 * var(acc) + var_rgsw + acc

    long double var = m_n * 8 * rlwe_p_prod_var;

    return var;
}

std::shared_ptr<MathWorker> CGGIBlindRotationContext::GetNTT() const {
    return m_ntt;
}

uint64_t CGGIBlindRotationContext::GetBlindRotationBasis() const {
    return m_basis;
}

uint64_t CGGIBlindRotationContext::GetBlindRotationBasisLog2() const {
    return m_basis_log2;
}

uint64_t CGGIBlindRotationContext::GetBlindRotationRGSWDigits() const {
    return m_digits;
}

uint64_t CGGIBlindRotationContext::GetModulus() const {
    return m_modulus;
}

double CGGIBlindRotationContext::GetStd() const {
    return m_std;
}

KeyDistribution CGGIBlindRotationContext::GetKeyDistribution() const {
    return m_distribution;
}

uint64_t CGGIBlindRotationContext::GetRingDimension() const {
    return m_N;
}

uint64_t CGGIBlindRotationContext::GetLWEDimension() const {
    return m_n;
}

OperatorID CGGIBlindRotationContext::GetOperatorID() const {
    return OperatorID::BR_CGGI;
}

Container CGGIBlindRotationContext::GetInputContainer() const {
    Container lwecont = std::make_shared<LWEContainerImpl>(m_n, 2 * m_N, 0.0);
    auto o_var = ComputeOutputVariance(0.0);
    Container rlwecont = std::make_shared<RLWEContainerImpl>(m_N, m_modulus, o_var);
    std::vector<Container> args = {lwecont, rlwecont};
    return std::make_shared<TupleContainerImpl>(args);
}

Container CGGIBlindRotationContext::GetOutputContainer(Container input) const {

    auto tp_container = std::dynamic_pointer_cast<TupleContainerImpl>(input);
    auto rlwe_container_gen = tp_container->GetElem(1);
    auto rlwe_container = std::dynamic_pointer_cast<RLWEContainerImpl>(rlwe_container_gen);

    auto var = ComputeOutputVariance();
    return std::make_shared<RLWEContainerImpl>(m_N, m_modulus, var + rlwe_container->GetVariance());
}

std::unique_ptr<CGGIBlindRotator> CGGIBlindRotationContext::ConstructOperator(const std::vector<GenericKey>& bundle) const {

    auto op = std::make_unique<CGGIBlindRotator>(shared_from_this());
    op->KeyGen(bundle[0].GetKeyPtr(), bundle[1].GetKeyPtr());

    return std::move(op);
}


CGGIBlindRotator::CGGIBlindRotator(const std::shared_ptr<const CGGIBlindRotationContext>& params) : m_params(params), m_params_set(true) {
    m_engine = params->GetNTT();
    m_mux = std::make_unique<MuxOperator>(m_engine, m_params->GetBlindRotationBasisLog2(), m_params->GetBlindRotationRGSWDigits());
    m_encryptor = std::make_shared<RLWEEncryptor>(m_engine, m_params->GetStd());
    m_accumulator.resize(2 * m_params->GetRingDimension());
    SetupMonomials();
}


void CGGIBlindRotator::KeyGen(const std::vector<uint64_t> &lwe_key, const std::vector<uint64_t> &rlwe_key) {
    assert(m_params_set);
    KeyGen(lwe_key.data(), rlwe_key.data());
}


void CGGIBlindRotator::BlindRotate(std::vector<uint64_t> &result, const std::vector<uint64_t> &lwe_vec,
                                   std::vector<uint64_t> &rlwe_acc_vec,
                                   bool output_as_coefficients) {
    assert(m_keys_generated);
    BlindRotate(result.data(), lwe_vec.data(), rlwe_acc_vec.data(), output_as_coefficients);
}

void CGGIBlindRotator::KeyGen(const uint64_t *lwe_key, const uint64_t *rlwe_key) {
    assert(m_params_set);


    switch (m_params->GetKeyDistribution()) {
        case BINARY: {
            KeyGenBinary(lwe_key, rlwe_key);
            break;
        }
        case TERNARY: {
            KeyGenTernary(lwe_key, rlwe_key);
            break;
        }
        default: {
            std::cerr << "Unsupported key distribution" << std::endl;
            assert(false);
        }
    }
    m_keys_generated = true;
}

void CGGIBlindRotator::BlindRotate(uint64_t *result, const uint64_t *lwe_vec, uint64_t *rlwe_acc_vec,
                                   bool output_as_coefficients) {

    assert(m_keys_generated);
    auto N = m_params->GetRingDimension();

    std::copy(rlwe_acc_vec, rlwe_acc_vec + 2 * N, result);

    switch (m_params->GetKeyDistribution()) {
        case BINARY: {
            BlindRotateBinary(lwe_vec, result);
            break;
        }
        case TERNARY: {
            BlindRotateTernary(lwe_vec, result);
            break;
        }
        default: {
            std::cerr << "Unsupported key distribution" << std::endl;
            assert(false);
        }
    }
    if (output_as_coefficients) {
        m_params->GetNTT()->BackwardNTT(result, result);
        m_params->GetNTT()->BackwardNTT(result + N, result + N);
    }
}

void CGGIBlindRotator::KeyGenBinary(const uint64_t *__restrict lwe_key, const uint64_t *__restrict rlwe_key) {
    auto n = m_params->GetLWEDimension();
    auto N = m_params->GetRingDimension();
    auto L = m_params->GetBlindRotationBasis();
    auto l = m_params->GetBlindRotationRGSWDigits();

    auto sk_ntt = AlignedBuffer(N);
    m_engine->ForwardNTT(sk_ntt.data(), rlwe_key);
    auto data = AlignedBuffer(N, 0);

    m_brk.resize(n * (4 * N * l));
    for(uint64_t i = 0; i < n; i++) {
        data[0] = lwe_key[i];
        auto rgsw_location = m_brk.data() + i * 4 * N * l;
        m_encryptor->MakeRGSW(rgsw_location, data.data(), sk_ntt.data(), L, l, false);
    }
}

void CGGIBlindRotator::KeyGenTernary(const uint64_t *__restrict lwe_key, const uint64_t *__restrict rlwe_key) {
    auto n = m_params->GetLWEDimension();
    auto N = m_params->GetRingDimension();
    auto L = m_params->GetBlindRotationBasis();
    auto l = m_params->GetBlindRotationRGSWDigits();

    auto sk_ntt = AlignedBuffer(N);
    m_engine->ForwardNTT(sk_ntt.data(), rlwe_key);
    auto data = AlignedBuffer(2 * N, 0);

    m_brk.resize(n * (2 * 4 * N * l));
    for(uint64_t i = 0; i < n; i++) {
        if (lwe_key[i] == 0) {
            data[0] = 0;
            data[N] = 0;
        } else {
            if (lwe_key[i] == 1) {
                data[0] = 1;
                data[N] = 0;
            } else {
                data[0] = 0;
                data[N] = 1;
            }
        }
        auto rgsw_location_lo = m_brk.data() + i * 8 * N * l;
        auto rgsw_location_hi = rgsw_location_lo + 4 * N * l;
        m_encryptor->MakeRGSW(rgsw_location_lo, data.data(), sk_ntt.data(), L, l, false);
        m_encryptor->MakeRGSW(rgsw_location_hi, data.data() + N, sk_ntt.data(), L, l, false);
    }
}

void CGGIBlindRotator::BlindRotateBinary(const uint64_t *const __restrict lwe_vec,  uint64_t *__restrict rlwe_vec) {
    // LWE dimension
    auto n = m_params->GetLWEDimension();
    // RLWE dimension
    auto N = m_params->GetRingDimension();
    // RLWE Q
    auto Q = m_params->GetModulus();
    // RGSW digits
    auto l = m_params->GetBlindRotationRGSWDigits();

    // pointer to scratch space
    auto m_acc_p = m_accumulator.data();
    // pointer to array containing [NTT(1), NTT(X), NTT(X^2), NTT(X^3), ...
    auto m_mon_p = m_monomials.data();
    // pointer to NTT(1)
    auto one_const_p = m_mon_p;
    // pointer to NTT(0)
    auto m_zero_p = m_mon_p + N * N;

    auto ntt = m_params->GetNTT();
    ntt->ForwardNTT(rlwe_vec, rlwe_vec);
    ntt->ForwardNTT(rlwe_vec + N, rlwe_vec + N);

    // Multiply by X^b
    uint64_t b = lwe_vec[n];
    if (b >= N) {
        // b >= N => X^b = -X^{b - N}
        b -= N;
        // m_acc = [-rlwe_A, -rlwe_B]
        ntt->SubModElw(m_acc_p, m_zero_p, rlwe_vec, 2 * N);
        // monomial ptr to X^{b - N}
        auto monomial = m_monomials.data() + b * N;
        // m_acc = [-rlwe_A * X^{b - N}, -rlwe_B * X^{b - N}]
        //       = [rlwe_A * X^b, rlwe_B * X^b]
        ntt->MulModElw(m_acc_p, m_acc_p, monomial, N);
        ntt->MulModElw(m_acc_p + N, m_acc_p + N, monomial, N);
    } else {
        // monomial ptr to X^b
        auto monomial = m_monomials.data() + b * N;
        // m_acc = [rlwe_A * X^b, rlwe_B * X^b]
        ntt->MulModElw(m_acc_p, rlwe_vec, monomial, N);
        ntt->MulModElw(m_acc_p + N, rlwe_vec + N, monomial, N);
    }

    // copy back to recycle scratch space
    std::copy(m_accumulator.begin(), m_accumulator.end(), rlwe_vec);

    // std::fill(m_accumulator.begin(), m_accumulator.end(), 0);
    ZERO_UINT64_ARR(m_accumulator.data(), 2 * N);

    // pointer to blind-rotation key
    auto m_brk_p = m_brk.data();

    for(uint64_t i = 0; i < n; i++) {
        auto a_i = lwe_vec[i];
        auto m_brk_i_p = m_brk_p + i * 4 * N * l;
        if (a_i == 0)
            continue;
        // since we do b - <a, s>, we negate
        auto a_i_neg = 2 * N - a_i;

        // set up the monomial
        // TODO is precomputing much faster ?
        if (a_i_neg >= N) {
            a_i_neg -= N;
            ntt->SubModElw(m_acc_p, m_zero_p, m_mon_p + N * a_i_neg, N);
        } else {
            std::copy(m_mon_p + a_i_neg * N, m_mon_p + a_i_neg * N + N, m_acc_p);
        }
        ntt->SubModElw(m_acc_p,m_acc_p, one_const_p, N);

        m_mux->BinaryCMux(rlwe_vec, m_brk_i_p, m_acc_p);
    }
}

void CGGIBlindRotator::BlindRotateTernary(const uint64_t *const __restrict lwe_vec, uint64_t *__restrict rlwe_vec) {
    // LWE dimension
    auto n = m_params->GetLWEDimension();
    // RLWE dimension
    auto N = m_params->GetRingDimension();
    // RLWE Q
    auto Q = m_params->GetModulus();
    // RGSW digits
    auto l = m_params->GetBlindRotationRGSWDigits();

    // pointer to scratch space
    auto m_acc_p = m_accumulator.data();
    // pointer to array containing [NTT(1), NTT(X), NTT(X^2), NTT(X^3), ...
    auto m_mon_p = m_monomials.data();
    // pointer to NTT(1)
    auto one_const_p = m_mon_p;
    // pointer to NTT(0)
    auto m_zero_p = m_mon_p + N * N;

    auto ntt = m_params->GetNTT();
    ntt->ForwardNTT(rlwe_vec, rlwe_vec);
    ntt->ForwardNTT(rlwe_vec + N, rlwe_vec + N);

    // Multiply by X^b
    uint64_t b = lwe_vec[n];
    if (b >= N) {
        // b >= N => X^b = -X^{b - N}
        b -= N;
        // m_acc = [-rlwe_A, -rlwe_B]
        ntt->SubModElw(m_acc_p, m_zero_p, rlwe_vec, 2 * N);
        // monomial ptr to X^{b - N}
        auto monomial = m_monomials.data() + b * N;
        // m_acc = [-rlwe_A * X^{b - N}, -rlwe_B * X^{b - N}]
        //       = [rlwe_A * X^b, rlwe_B * X^b]
        ntt->MulModElw(m_acc_p, m_acc_p, monomial, N);
        ntt->MulModElw(m_acc_p + N, m_acc_p + N, monomial, N);
    } else {
        // monomial ptr to X^b
        auto monomial = m_monomials.data() + b * N;
        // m_acc = [rlwe_A * X^b, rlwe_B * X^b]
        ntt->MulModElw(m_acc_p, rlwe_vec, monomial, N);
        ntt->MulModElw(m_acc_p + N, rlwe_vec + N, monomial, N);
    }

    // copy back to recycle scratch space
    std::copy(m_accumulator.begin(), m_accumulator.end(), rlwe_vec);
    //std::fill(m_accumulator.begin(), m_accumulator.end(), 0);
    ZERO_UINT64_ARR(m_accumulator.data(), 2 * N);

    // pointer to blind-rotation key
    auto m_brk_p = m_brk.data();

    for(uint64_t i = 0; i < n; i++) {
        auto a_i = lwe_vec[i];
        auto m_brk_i_p = m_brk_p + i * 8 * N * l;
        if (a_i == 0)
            continue;
        // since we do b - <a, s>, we negate
        auto a_i_neg = 2 * N - a_i;

        // set up the monomials, we aim for [X^{-a_i} - 1, X^{a_i} - 1]
        // TODO is precomputing much faster ?
        if (a_i_neg >= N) {
            a_i_neg -= N;
            // m_acc_p ptr to [-X^{a_i_neg - N}, 0 ]
            ntt->SubModElw(m_acc_p, m_zero_p, m_mon_p + N * a_i_neg, N);
            // m_acc_p ptr to [-X^{a_i_neg - N}, X^{a_i} ]
            //              = [X^{a_i_neg}, X^{a_i}]
            if (a_i == N) [[unlikely]] {
                std::copy(m_acc_p, m_acc_p + N, m_acc_p + N);
            } else {
                std::copy(m_mon_p + N * a_i, m_mon_p + N * a_i + N, m_acc_p + N);
            }
        } else {
            a_i -= N;
            std::copy(m_mon_p + a_i_neg * N, m_mon_p + a_i_neg * N + N, m_acc_p);
            ntt->SubModElw(m_acc_p + N, m_zero_p + N, m_mon_p + N * a_i, N);
        }
        // Subtract 1
        ntt->SubModElw(m_acc_p,m_acc_p, one_const_p, N);
        ntt->SubModElw(m_acc_p + N,m_acc_p + N, one_const_p, N);

        m_mux->TernaryCMux(rlwe_vec, m_brk_i_p, m_acc_p);
    }
}

void CGGIBlindRotator::SetupMonomials() {
    auto N = m_params->GetRingDimension();

    m_monomials.resize(N * (N + 2));
    std::fill(m_monomials.begin(), m_monomials.end(), 0);
    for(uint64_t i = 0; i < N; i++) {
        auto idx = i * N;
        m_monomials[idx + i] = 1;
        m_engine->ForwardNTT(m_monomials.data() + idx, m_monomials.data() + idx);
    }

}

std::shared_ptr<RLWEEncryptor> CGGIBlindRotator::GetEncryptor() const {
    return m_encryptor;
}

const std::shared_ptr<const OperatorContext<BlindRotator>> CGGIBlindRotator::GetContext() const {
    return std::reinterpret_pointer_cast<const OperatorContext<BlindRotator>>(m_params);
}