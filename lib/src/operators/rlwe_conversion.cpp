//
// Created by leonard on 4/26/26.
//

#include <cassert>
#include <iostream>
#include "base_crypto.h"
#include "operators/endo_glwe_conversion.h"
#include "utils/generic_utils.h"
#include "utils/math_utils.h"
#include "static/gadget_decomp.h"

RLWEConversionContext::RLWEConversionContext(KeyDistribution source_key_distribution, uint64_t modulus, uint64_t N,
                                             uint64_t basis, uint64_t digits, double std) :
                                                    m_source_distribution(source_key_distribution),
                                                    m_modulus(modulus),
                                                    m_N(N),
                                                    m_basis(basis),
                                                    m_digits(digits),
                                                    m_std(std) {

    m_basis_log2 = IntLog2(m_basis);
    m_ntt = SelectWorker(modulus, N);

}

RLWEConversionContext::RLWEConversionContext(KeyDistribution source_key_distribution,
                                             std::shared_ptr<MathWorker> ntt, uint64_t basis, uint64_t digits,
                                             double std) :
        m_source_distribution(source_key_distribution),
        m_modulus(ntt->GetModulus()),
        m_N(ntt->GetDimension()),
        m_basis(basis),
        m_digits(digits),
        m_std(std),
        m_ntt(ntt) {

    m_basis_log2 = IntLog2(m_basis);

}

RLWEConversionContext::RLWEConversionContext(const RLWEConversionContext &other) :
        enable_shared_from_this(other),
m_source_distribution(other.GetSourceKeyDistribution()),
m_modulus(other.GetModulus()),
m_N(other.GetDimension()),
m_basis(other.GetGadgetBasis()),
m_basis_log2(other.GetGadgetBasisLog2()),
m_digits(other.GetGadgetDigits()),
m_std(other.GetStd()),
m_ntt(other.GetNTT()){
}


long double RLWEConversionContext::ComputeOutputVariance(long double input_variance) const {

    long double delta = std::pow(2.0, IntLog2(m_modulus)) / std::pow(m_basis, m_digits);
    if (delta <= 1) {
        delta = 0.0;
    }
    delta *= delta / 12.0;

    // assume worst case
    long double message_norm_square = 0.0;
    switch (m_source_distribution) {

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

    long double var = rlwe_p_prod_var + delta * message_norm_square;
    return input_variance + var;
}

KeyDistribution RLWEConversionContext::GetSourceKeyDistribution() const {
    return m_source_distribution;
}

uint64_t RLWEConversionContext::GetModulus() const {
    return m_modulus;
}
uint64_t RLWEConversionContext::GetDimension() const {
    return m_N;
}

uint64_t RLWEConversionContext::GetGadgetBasis() const {
    return m_basis;
}

uint64_t RLWEConversionContext::GetGadgetBasisLog2() const {
    return m_basis_log2;
}

uint64_t RLWEConversionContext::GetGadgetDigits() const {
    return m_digits;
}

double RLWEConversionContext::GetStd() const {
    return m_std;
}

std::shared_ptr<MathWorker> RLWEConversionContext::GetNTT() const {
    return m_ntt;
}

void RLWEConversionContext::SetSourceKeyDistribution(KeyDistribution distribution) {
    m_source_distribution = distribution;
}

void RLWEConversionContext::SetStd(double std) {
    m_std = std;
}

void RLWEConversionContext::SetDimension(uint64_t input_dimension) {
    m_N = input_dimension;
}

void RLWEConversionContext::SetNTT(std::shared_ptr<MathWorker> ntt) {
    m_ntt = std::move(ntt);
}

void RLWEConversionContext::SetModulus(uint64_t modulus) {
    m_modulus = modulus;
}

void RLWEConversionContext::SetGadgetDigits(uint64_t digits) {
    m_digits = digits;
}

void RLWEConversionContext::SetGadgetBasis(uint64_t basis) {
    m_basis = basis;
    m_basis_log2 = IntLog2(m_basis);
}

Container RLWEConversionContext::GetInputContainer() const {
    return std::make_shared<RLWEContainerImpl>(m_modulus, m_N, 0.0);
}

Container RLWEConversionContext::GetOutputContainer(Container input) const {
    auto in = std::dynamic_pointer_cast<RLWEContainerImpl>(input);
    auto out_var = ComputeOutputVariance(in->GetVariance());
    return std::make_shared<RLWEContainerImpl>(m_N, m_modulus, out_var);
}

OperatorID RLWEConversionContext::GetOperatorID() const {
    return CONV_RLWE_RLWE;
}

std::unique_ptr<RLWEtoRLWEConverter> RLWEConversionContext::ConstructOperator(const std::vector<GenericKey> &keys) const {
    auto op = std::unique_ptr<RLWEtoRLWEConverter>(new RLWEtoRLWEConverter(this->shared_from_this()));
    op->KeyGen(keys[0].GetKeyPtr(), keys[1].GetKeyPtr());

    return std::move(op);
}


RLWEtoRLWEConverter::RLWEtoRLWEConverter(std::shared_ptr<const RLWEConversionContext> params) : m_params(params) {
    auto N = params->GetDimension();
    auto digits = params->GetGadgetDigits();

    m_acc.resize(2 * N * digits);
    m_ksk.resize(2 * N * digits);
    m_params_set = true;
}


void RLWEtoRLWEConverter::KeyGen(const std::vector<uint64_t> &source_key, const std::vector<uint64_t> &target_key) {
    assert(m_params_set);
    m_keys_generated = true;
    KeyGen(source_key.data(), target_key.data());
}

void RLWEtoRLWEConverter::KeyGen(const uint64_t *const source_key, const uint64_t *const target_key) {
    assert(m_params_set);
    m_keys_generated = true;
    // assumption: both keys in non-ntt form

    auto N = m_params->GetDimension();
    auto modulus = m_params->GetModulus();
    auto modulus_bits = IntLog2(m_params->GetModulus());
    auto basis_bits = m_params->GetGadgetBasisLog2();
    auto digits = m_params->GetGadgetDigits();
    auto ntt = m_params->GetNTT();
    RLWEEncryptor encryptor(ntt, m_params->GetStd());

    AlignedBuffer target_key_ntt(N, 0);
    AlignedBuffer gadget_entry(N, 0);
    ntt->ForwardNTT(target_key_ntt.data(), target_key);

    auto ksk = m_ksk.data();
    auto g_ij = 1ull << (modulus_bits - basis_bits);

    for(uint64_t i = 0; i < digits; i++) {
        auto chunk = ksk + i * 2 * N;
        ntt->MulSModElw(gadget_entry.data(), source_key, g_ij, N);
        //intel::hexl::EltwiseFMAMod(gadget_entry.data(), source_key, g_ij, nullptr, N, modulus, 1);
        encryptor.MakeRLWE(chunk, gadget_entry.data(), target_key_ntt.data());
        g_ij >>= basis_bits;
        if (g_ij == 0) {
            g_ij = 1;
        }
    }

}

void RLWEtoRLWEConverter::Convert(std::vector<uint64_t> &output, const std::vector<uint64_t> &input) {
    Convert(output.data(), input.data());
}

void RLWEtoRLWEConverter::Convert(uint64_t *output, const uint64_t *const input) {
    auto modulus = m_params->GetModulus();
    auto N = m_params->GetDimension();
    auto ntt = m_params->GetNTT();
    auto digits = m_params->GetGadgetDigits();
    auto basis_bits = m_params->GetGadgetBasisLog2();
    auto modulus_bits = IntLog2(modulus);
    auto acc = m_acc.data();
    ZERO_UINT64_ARR(acc, m_acc.size());

    SignedDigitDecomposeRep2NTT(acc, input, ntt, digits, basis_bits, modulus_bits);

    ntt->MulEqModElw(acc, m_ksk.data(), 2 * N * digits);
    //intel::hexl::EltwiseMultMod(acc, acc, m_ksk.data(), 2 * N * digits, modulus, 1);

    for(uint64_t i = 1; i < digits; i++) {
        ntt->AddEqModElw(acc, acc + i * 2 * N, 2 * N);
        //intel::hexl::EltwiseAddMod(acc, acc, acc + i * 2 * N, 2 * N, modulus);
    }
    ntt->SubEqModElw(output, acc, 2 * N);
    //intel::hexl::EltwiseSubMod(output, output, acc, 2 * N, modulus);

    // DANGER: Correct provided b @input + N is in NTT format
    ntt->AddEqModElw(output + N, input + N, N);
    //intel::hexl::EltwiseAddMod(output + N, output + N, input + N, N, modulus);
}

const std::shared_ptr<const RLWEConversionContext> RLWEtoRLWEConverter::GetContext() const {
    return m_params;
}