//
// Created by leonard on 11/12/24.
//

#include <cassert>
#include <utility>

#include "utils/math_utils.h"
#include "utils/generic_utils.h"

#include "operators/automorphism_evaluation.h"

AutomorphismContext::AutomorphismContext(KeyDistribution source_key_distribution, uint64_t modulus, uint64_t N,
                                         uint64_t basis, uint64_t digits, double std,
                                         uint32_t automorphism_index) :  m_automorphism_index(automorphism_index) {
    m_rlwe_conversion = std::make_shared<RLWEConversionContext>(source_key_distribution, modulus, N, basis, digits, std);
}

AutomorphismContext::AutomorphismContext(KeyDistribution source_key_distribution,
                                         std::shared_ptr<MathWorker> ntt, uint64_t basis, uint64_t digits,
                                         double std, uint32_t automorphism_index) : m_automorphism_index(automorphism_index) {
    m_rlwe_conversion = std::make_shared<RLWEConversionContext>(source_key_distribution, ntt, basis, digits, std);
}

AutomorphismContext::AutomorphismContext(const AutomorphismContext &other)  : enable_shared_from_this(other) {
    auto ptr = other.m_rlwe_conversion;
    m_rlwe_conversion = std::make_shared<RLWEConversionContext>(*ptr);
    m_automorphism_index = other.GetAutomorphismIndex();
}

long double AutomorphismContext::ComputeOutputVariance(long double input_variance) const {
    return m_rlwe_conversion->ComputeOutputVariance(input_variance);
}

Container AutomorphismContext::GetInputContainer() const {
    return std::make_shared<RLWEContainerImpl>(m_rlwe_conversion->GetDimension(), m_rlwe_conversion->GetModulus(), 0.0);
}

Container AutomorphismContext::GetOutputContainer(Container input) const {
    return m_rlwe_conversion->GetOutputContainer(input);
}

std::unique_ptr<AutomorphismEvaluator> AutomorphismContext::ConstructOperator(const std::vector<GenericKey> &keys) const {
    auto N = m_rlwe_conversion->GetDimension();
    auto Q = m_rlwe_conversion->GetModulus();
    auto& key = keys[0].GetKey();

    std::vector<uint64_t> key_default(N, 0);
    std::vector<uint64_t> key_auto(N, 0);

    std::copy(key.begin(), key.end(), key_default.data());
    EvalNegacyclicAutomorphism(key_auto.data(), key_default.data(), m_automorphism_index, N, Q);

    std::vector<GenericKey> bundle = {
            {"AUTO_KEY", key_auto.data(), N },
            {"DEFAULT_KEY", key_default.data(), N}
    };

    auto rlwe_converter = m_rlwe_conversion->ConstructOperator(bundle);
    auto op = std::unique_ptr<AutomorphismEvaluator>(new AutomorphismEvaluator(shared_from_this(), std::move(rlwe_converter)));

    return std::move(op);
}

OperatorID AutomorphismContext::GetOperatorID() const {
    return EVAL_AUTO;
}

void AutomorphismContext::SetAutomorphismIndex(uint32_t idx) {
    m_automorphism_index = idx;
}

uint32_t AutomorphismContext::GetAutomorphismIndex() const {
    return m_automorphism_index;
}

KeyDistribution AutomorphismContext::GetSourceKeyDistribution() const {
    return m_rlwe_conversion->GetSourceKeyDistribution();
}

uint64_t AutomorphismContext::GetModulus() const {
    return m_rlwe_conversion->GetModulus();
}

uint64_t AutomorphismContext::GetDimension() const {
    return m_rlwe_conversion->GetDimension();
}

uint64_t AutomorphismContext::GetGadgetBasis() const {
    return m_rlwe_conversion->GetGadgetBasis();
}

uint64_t AutomorphismContext::GetGadgetBasisLog2() const {
    return m_rlwe_conversion->GetGadgetBasisLog2();
}

uint64_t AutomorphismContext::GetGadgetDigits() const {
    return m_rlwe_conversion->GetGadgetDigits();
}

double AutomorphismContext::GetStd() const {
    return m_rlwe_conversion->GetStd();
}

std::shared_ptr<MathWorker> AutomorphismContext::GetNTT() const {
    return m_rlwe_conversion->GetNTT();
}

void AutomorphismContext::SetSourceKeyDistribution(KeyDistribution distribution) {
    m_rlwe_conversion->SetSourceKeyDistribution(distribution);
}

void AutomorphismContext::SetModulus(uint64_t modulus) {
    m_rlwe_conversion->SetModulus(modulus);
}

void AutomorphismContext::SetDimension(uint64_t input_dimension) {
    m_rlwe_conversion->SetDimension(input_dimension);
}

void AutomorphismContext::SetGadgetBasis(uint64_t basis) {
    m_rlwe_conversion->SetGadgetBasis(basis);
}

void AutomorphismContext::SetGadgetDigits(uint64_t digits) {
    m_rlwe_conversion->SetGadgetDigits(digits);
}

void AutomorphismContext::SetStd(double std) {
    m_rlwe_conversion->SetStd(std);
}

void AutomorphismContext::SetNTT(std::shared_ptr<MathWorker> ntt) {
    m_rlwe_conversion->SetNTT(std::move(ntt));
}





AutomorphismEvaluator::AutomorphismEvaluator(std::shared_ptr<const AutomorphismContext> params, std::unique_ptr<RLWEtoRLWEConverter> conv)
   : m_automorphism_params(std::move(params)), m_converter(std::move(conv)) {
    m_params_set = true;
    m_keys_generated = true;
    m_auto_buffer.resize(2 * m_automorphism_params->GetDimension());
}



void AutomorphismEvaluator::Eval(uint64_t *output, const uint64_t *const input) {
    assert(m_keys_generated);
    auto m_ntt = m_automorphism_params->GetNTT();
    auto N = m_ntt->GetDimension();
    ApplyAutomorphism(m_auto_buffer.data(), input);
    ApplyAutomorphism(m_auto_buffer.data() + N , input + N);
    // TODO Remove me
    m_ntt->ForwardNTT(m_auto_buffer.data() + N, m_auto_buffer.data() + N);
    m_converter->Convert(output, m_auto_buffer.data());
}

void AutomorphismEvaluator::Eval(std::vector<uint64_t> &output, const std::vector<uint64_t> &input) {
    assert(m_keys_generated);
    auto m_ntt = m_automorphism_params->GetNTT();
    auto N = m_ntt->GetDimension();
    ApplyAutomorphism(m_auto_buffer.data(), input.data());
    ApplyAutomorphism(m_auto_buffer.data() + N , input.data() + N);
    // TODO Remove me
    m_ntt->ForwardNTT(m_auto_buffer.data() + N, m_auto_buffer.data() + N);
    m_converter->Convert(output.data(), m_auto_buffer.data());
}

void AutomorphismEvaluator::EvalSwitchOnly(uint64_t *output, const uint64_t *const input) {
    assert(m_keys_generated);
    ZERO_UINT64_ARR(m_auto_buffer.data(), m_auto_buffer.size());
    ApplyAutomorphism(m_auto_buffer.data(), input);
    m_converter->Convert(output, m_auto_buffer.data());
}

void AutomorphismEvaluator::ApplyAutomorphism(std::vector<uint64_t> &output, const std::vector<uint64_t> &input) {
    ApplyAutomorphism(output.data(), input.data());
}

void AutomorphismEvaluator::ApplyAutomorphism(uint64_t *output, const uint64_t *const input) {
    // TODO benchmark
    auto N = m_automorphism_params->GetDimension();
    auto dim2_mask = 2 * N - 1;
    auto dim_mask = N - 1;
    auto Q = m_automorphism_params->GetModulus();
    auto m_auto_idx = m_automorphism_params->GetAutomorphismIndex();
    const uint32_t stride = 4;


    for(uint64_t i = 0; i < N; i+=stride) {
        for(uint32_t j = 0; j < stride; j++) {
            auto idx = i + j;
            auto new_idx = m_auto_idx * idx;
            auto flip_sign = (new_idx & dim2_mask) >= N;
            auto poly_idx = new_idx & dim_mask;
            output[poly_idx] = input[idx];

            if (flip_sign) {
                output[poly_idx] = output[poly_idx] == 0 ? 0 : Q - output[poly_idx]; // intel::hexl::SubUIntMod(0, output[poly_idx], Q);
            }
        }
    }

}
