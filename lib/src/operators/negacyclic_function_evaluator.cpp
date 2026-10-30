//
// Created by leonard on 7/27/26.
//

#include <iostream>
#include "operators/negacyclic_function_evaluator.h"
#include "utils/generic_utils.h"
#include "static/sample_extraction.h"

NegacyclicFunctionEvaluationContext::NegacyclicFunctionEvaluationContext(
        std::shared_ptr<OperatorContext<BlindRotator>> blind_rotation_context) : m_rotation_context(std::move(blind_rotation_context)) {

}

std::shared_ptr<OperatorContext<BlindRotator>> NegacyclicFunctionEvaluationContext::GetBlindRotationContext() const {
    return m_rotation_context;
}



OperatorID NegacyclicFunctionEvaluationContext::GetOperatorID() const {
    return NEGACYCLIC_FUNC_BOOT;
}

Container NegacyclicFunctionEvaluationContext::GetInputContainer() const {

    // todo: revisit
    auto br_input_gen = m_rotation_context->GetInputContainer();
    auto br_input = std::dynamic_pointer_cast<TupleContainerImpl>(br_input_gen);

    return br_input->GetElem(0);
}

Container NegacyclicFunctionEvaluationContext::GetOutputContainer(Container input) const {

    auto br_input_gen = m_rotation_context->GetInputContainer();

    // TODO

    return br_input_gen;
}

long double NegacyclicFunctionEvaluationContext::ComputeOutputVariance(long double input_variance) const {
    auto post_br_variance = m_rotation_context->ComputeOutputVariance(input_variance);

    return post_br_variance;
}

std::unique_ptr<NegacyclicFunctionEvaluator>
NegacyclicFunctionEvaluationContext::ConstructOperator(const std::vector<GenericKey> &keys) const {

    auto blind_rotator = m_rotation_context->ConstructOperator(keys);

    std::vector<GenericKey> switch_keys;
    switch_keys.push_back(keys[1]);
    switch_keys.push_back(keys[0]);

    return std::unique_ptr<NegacyclicFunctionEvaluator>(new NegacyclicFunctionEvaluator(shared_from_this(), std::move(blind_rotator)));
}

NegacyclicFunctionEvaluator::NegacyclicFunctionEvaluator(std::shared_ptr<const NegacyclicFunctionEvaluationContext> ctx,
                                                         std::unique_ptr<BlindRotator> rotator) :
                                         m_context(std::move(ctx)), m_rotator(std::move(rotator))    {
    auto cont_in = m_context->GetBlindRotationContext()->GetInputContainer();
    auto cont_rlwe = std::dynamic_pointer_cast<RLWEContainerImpl>(cont_in);
    m_acc_buffer.resize(4 * cont_rlwe->GetN());
    std::fill(m_acc_buffer.begin(), m_acc_buffer.end(), 0);
}

void NegacyclicFunctionEvaluator::EvalFunc(uint64_t *__restrict result, const uint64_t *__restrict input_lwe,
                                           std::function<uint64_t(uint64_t)> &function) {
    auto rot_ctx = m_rotator->GetContext();
    auto br_output = std::dynamic_pointer_cast<RLWEContainerImpl>(rot_ctx->GetOutputContainer(rot_ctx->GetInputContainer()));
    auto N = br_output->GetN();
    auto Q = br_output->GetQ();

    auto ntt = SelectWorker(Q, N);

    AlignedBuffer tmp_rlwe(2 * N);
    std::fill(tmp_rlwe.begin(), tmp_rlwe.end(), 0);

    // initial [f(0), f(N - 1), ..., f(1)]
    tmp_rlwe[0] = function(0);
    for(uint64_t i = 1; i < N; i++) {
        tmp_rlwe[i + N] = ntt->SubMod(0, function(N-i));
    }

    EvalFunc(result, input_lwe, tmp_rlwe.data());
}


void NegacyclicFunctionEvaluator::EvalFunc(uint64_t * RESTRICTED result, const uint64_t* RESTRICTED input, uint64_t *RESTRICTED lut_rlwe) {
    auto rot_ctx = m_rotator->GetContext();
    auto br_output = std::dynamic_pointer_cast<RLWEContainerImpl>(rot_ctx->GetOutputContainer(rot_ctx->GetInputContainer()));
    auto N = br_output->GetN();
    auto Q = br_output->GetQ();

    m_rotator->BlindRotate(m_acc_buffer.data(), input, lut_rlwe, true);
    SampleExtract(result, m_acc_buffer.data(), 0, N, Q);

}

const std::shared_ptr<const NegacyclicFunctionEvaluationContext> NegacyclicFunctionEvaluator::GetContext() const {
    return m_context;
}

void NegacyclicFunctionEvaluator::EvalFunc(std::vector<uint64_t> &result, const std::vector<uint64_t>& input, std::vector<uint64_t> &lut_rlwe) {
    EvalFunc(result.data(), input.data(), lut_rlwe.data());
}

void NegacyclicFunctionEvaluator::EvalFunc(std::vector<uint64_t> &result, const std::vector<uint64_t>& input, std::function<uint64_t(uint64_t)> &function) {
    EvalFunc(result.data(), input.data(), function);
}



