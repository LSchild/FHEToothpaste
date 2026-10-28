//
// Created by leonard on 5/30/26.
//


#include <ranges>
#include <utility>

#include "operators/trace_evaluation.h"
#include "utils/generic_utils.h"
#include "utils/math_utils.h"

#include <iostream>

TraceEvaluationContext::TraceEvaluationContext(KeyDistribution source_key_distribution, uint64_t modulus, uint64_t N,
                                               uint64_t basis, uint64_t digits, double std) {
    m_auto_context = std::make_shared<AutomorphismContext>(source_key_distribution, modulus, N, basis, digits, std, 3);

}

TraceEvaluationContext::TraceEvaluationContext(KeyDistribution source_key_distribution,
                                               std::shared_ptr<MathWorker> ntt, uint64_t basis, uint64_t digits,
                                               double std) {
    m_auto_context = std::make_shared<AutomorphismContext>(source_key_distribution, ntt, basis, digits, std, 3);
}

TraceEvaluationContext::TraceEvaluationContext(const TraceEvaluationContext &other) {
    m_auto_context = std::make_shared<AutomorphismContext>(other.GetSourceKeyDistribution(), other.GetNTT(), other.GetGadgetBasis(), other.GetGadgetDigits(), other.GetStd(), 3);
}

Container TraceEvaluationContext::GetInputContainer() const {
    return std::make_shared<RLWEContainerImpl>(m_auto_context->GetDimension(), m_auto_context->GetModulus(), 0.0);
}

long double TraceEvaluationContext::ComputeOutputVariance(long double input_variance) const {
    long double N = m_auto_context->GetDimension();
    auto auto_result_rlwe = std::dynamic_pointer_cast<RLWEContainerImpl>(m_auto_context->GetInputContainer());
    return auto_result_rlwe->GetVariance() * (std::pow(N, 2) / 3);

}

Container TraceEvaluationContext::GetOutputContainer(Container input) const {
    long double N = m_auto_context->GetDimension();
    auto input_rlwe =  std::dynamic_pointer_cast<RLWEContainerImpl>(m_auto_context->GetOutputContainer(input));
    // todo: implement the better one
    auto o_var = ComputeOutputVariance(input_rlwe->GetVariance());
    return std::make_shared<RLWEContainerImpl>(N, m_auto_context->GetModulus(), o_var);
}

std::unique_ptr<TraceEvaluator> TraceEvaluationContext::ConstructOperator(const std::vector<GenericKey> &keys) const {
    std::vector<std::unique_ptr<AutomorphismEvaluator>> evaluators;
    auto N = GetDimension();
    for(uint32_t i = 2; i <= N; i *= 2) {
        auto tmp_context = std::make_shared<AutomorphismContext>(*m_auto_context);
        tmp_context->SetAutomorphismIndex(i + 1);
        evaluators.push_back(tmp_context->ConstructOperator(keys));
    }

    auto op = std::unique_ptr<TraceEvaluator>(new TraceEvaluator(shared_from_this(), std::move(evaluators)));

    return op;
}

OperatorID TraceEvaluationContext::GetOperatorID() const {
    return EVAL_TRACE;
}

KeyDistribution TraceEvaluationContext::GetSourceKeyDistribution() const {
    return m_auto_context->GetSourceKeyDistribution();
}

uint64_t TraceEvaluationContext::GetModulus() const {
    return m_auto_context->GetModulus();
}

uint64_t TraceEvaluationContext::GetDimension() const {
    return m_auto_context->GetDimension();
}

uint64_t TraceEvaluationContext::GetGadgetBasis() const {
    return m_auto_context->GetGadgetBasis();
}

uint64_t TraceEvaluationContext::GetGadgetBasisLog2() const {
    return m_auto_context->GetGadgetBasisLog2();
}

uint64_t TraceEvaluationContext::GetGadgetDigits() const {
    return m_auto_context->GetGadgetDigits();
}

double TraceEvaluationContext::GetStd() const {
    return m_auto_context->GetStd();
}

std::shared_ptr<MathWorker> TraceEvaluationContext::GetNTT() const {
    return m_auto_context->GetNTT();
}

void TraceEvaluationContext::SetSourceKeyDistribution(KeyDistribution distribution) {
    m_auto_context->SetSourceKeyDistribution(distribution);
}

void TraceEvaluationContext::SetModulus(uint64_t modulus) {
    m_auto_context->SetModulus(modulus);
}

void TraceEvaluationContext::SetDimension(uint64_t input_dimension) {
    m_auto_context->SetDimension(input_dimension);
}

void TraceEvaluationContext::SetGadgetBasis(uint64_t basis) {
    m_auto_context->SetGadgetBasis(basis);
}

void TraceEvaluationContext::SetGadgetDigits(uint64_t digits) {
    m_auto_context->SetGadgetDigits(digits);
}

void TraceEvaluationContext::SetStd(double std) {
    m_auto_context->SetStd(std);
}

void TraceEvaluationContext::SetNTT(std::shared_ptr<MathWorker> ntt) {
    m_auto_context->SetNTT(std::move(ntt));
}





TraceEvaluator::TraceEvaluator(std::shared_ptr<const TraceEvaluationContext> context,
                               std::vector<std::unique_ptr<AutomorphismEvaluator>> evaluators) : m_context(std::move(context)), m_trace_evaluators(std::move(evaluators)) {}

void TraceEvaluator::Eval(uint64_t *output, const uint64_t *input) {

    auto ntt = m_context->GetNTT();
    auto N = m_context->GetDimension();
    auto Q = m_context->GetModulus();
    auto N_inverse = ntt->InvMod(N);

    // TODO: Preallocate ?
    AlignedBuffer buffer(2 * N, 0);
    const auto buf_p = buffer.data();
    ntt->MulSModElw(output, input, N_inverse, 2 * N);

    for(auto & m_trace_evaluator : std::views::reverse(m_trace_evaluators)) {

        std::cerr << std::endl;
        m_trace_evaluator->Eval(buf_p, output);
        ntt->BackwardNTT(buf_p, buf_p);
        ntt->BackwardNTT(buf_p + N, buf_p + N);

        ntt->AddEqModElw(output, buf_p, 2 * N);

        ZERO_UINT64_ARR(buf_p, 2 * N);
    }

    ntt->ForwardNTT(output, output);
    ntt->ForwardNTT(output + N, output + N);

}

void TraceEvaluator::Eval(std::vector<uint64_t> &output, const std::vector<uint64_t> &input) {
    Eval(output.data(), input.data());
}

std::shared_ptr<const TraceEvaluationContext> TraceEvaluator::GetContext() const {
    return m_context;
}

void TraceEvaluator::EvalAuto(uint64_t *output, const uint64_t *input, uint64_t auto_idx) {
    // checks that the automorphism index is power of 2 + 1
    auto index = IntLog2(auto_idx - 1);
    assert(((1 << index) + 1) == auto_idx);
    assert(index <= m_trace_evaluators.size());

    m_trace_evaluators[index-1]->Eval(output, input);
}

void TraceEvaluator::EvalAuto(std::vector<uint64_t> &output, const std::vector<uint64_t> &input, uint64_t auto_idx) {
    // checks that the automorphism index is power of 2 + 1
    auto index = IntLog2(auto_idx - 1);
    assert(((1 << index) + 1) == auto_idx);
    assert(index < m_trace_evaluators.size());

    m_trace_evaluators[index]->Eval(output, input);
}