//
// Created by leonard on 7/31/26.
//

#include <numeric>
#include <utility>
#include "operators/ks_function_evaluator.h"

KSFunctionEvaluationContext::KSFunctionEvaluationContext(std::shared_ptr<OperatorContext<BlindRotator>> blind_rotation_context,
        std::shared_ptr<LWEConversionContext> converter, uint64_t padding_width) : m_blind_rotation_context(std::move(blind_rotation_context)),
        m_conversion_context(std::move(converter)), m_estimate_padding(false), m_padding_width(padding_width) {

    auto br_input = blind_rotation_context->GetInputContainer();
    auto br_input_tuple = std::dynamic_pointer_cast<TupleContainerImpl>(br_input);
    auto acc_input = br_input_tuple->GetElem(1);
    auto rlwe_input = std::dynamic_pointer_cast<RLWEContainerImpl>(acc_input);

    auto N = rlwe_input->GetN();
    auto Q = rlwe_input->GetQ();

    m_worker = SelectWorker(Q, N);

}

KSFunctionEvaluationContext::KSFunctionEvaluationContext(
        std::shared_ptr<OperatorContext<BlindRotator>> blind_rotation_context,
        std::shared_ptr<LWEConversionContext> converter, bool estimate_padding) : m_blind_rotation_context(std::move(blind_rotation_context))
        , m_conversion_context(std::move(converter)), m_estimate_padding(estimate_padding), m_padding_width(1) {

    auto br_input = blind_rotation_context->GetInputContainer();
    auto br_input_tuple = std::dynamic_pointer_cast<TupleContainerImpl>(br_input);
    auto acc_input = br_input_tuple->GetElem(1);
    auto rlwe_input = std::dynamic_pointer_cast<RLWEContainerImpl>(acc_input);

    auto N = rlwe_input->GetN();
    auto Q = rlwe_input->GetQ();

    m_worker = SelectWorker(Q, N);

}

std::unique_ptr<KSFunctionEvaluator>
KSFunctionEvaluationContext::ConstructOperator(const std::vector<GenericKey> &keys) const {

    auto rotator = m_blind_rotation_context->ConstructOperator(keys);

    std::vector<GenericKey> conversion_keys = {keys[1], keys[0]};

    auto converter = m_conversion_context->ConstructOperator(conversion_keys);

    return std::unique_ptr<KSFunctionEvaluator>(new KSFunctionEvaluator(shared_from_this(), std::move(rotator), std::move(converter), m_estimate_padding, m_padding_width));
}

OperatorID KSFunctionEvaluationContext::GetOperatorID() const {
    return FUNC_BOOT_KS;
}

std::shared_ptr<MathWorker> KSFunctionEvaluationContext::GetWorker() const {
    return m_worker;
}

std::shared_ptr<OperatorContext<BlindRotator>> KSFunctionEvaluationContext::GetBlindRotationContext() const {
    return m_blind_rotation_context;
}

std::shared_ptr<LWEConversionContext> KSFunctionEvaluationContext::GetLWEConversionContext() const {
    return m_conversion_context;
}

long double KSFunctionEvaluationContext::ComputeOutputVariance(long double input_variance) const {
    // todo
    return 0;
}

Container KSFunctionEvaluationContext::GetInputContainer() const {
    // todo
    return {};
}

Container KSFunctionEvaluationContext::GetOutputContainer(Container input) const {
    // todo
    return {};
}


KSFunctionEvaluator::KSFunctionEvaluator(std::shared_ptr<const KSFunctionEvaluationContext> context,
                                         std::unique_ptr<BlindRotator> rotator,
                                         std::unique_ptr<LWEtoLWEConverter> converter, bool estimate_padding,
                                         uint64_t padding_width) : m_context(std::move(context)), m_rotator(std::move(rotator)), m_converter(std::move(converter)), m_estimate_padding(estimate_padding), m_padding_width(padding_width) {}

uint64_t KSFunctionEvaluator::EstimatePadding(const uint64_t *poly) const {

   auto in = m_context->GetBlindRotationContext()->GetInputContainer();
   auto out = std::dynamic_pointer_cast<RLWEContainerImpl>(m_context->GetOutputContainer(in));

   auto N = out->GetN();

   std::vector<uint64_t> block_sizes;
   block_sizes.reserve(N);
   block_sizes.emplace_back(1);

   auto current_block_index = 0;

   for(uint64_t i = 1; i < N; i++) {
       // tee-hee I like to be """clever"""
       if (poly[i] == poly[i - 1])
           block_sizes[current_block_index]++;
       else {
           ++current_block_index;
           block_sizes.emplace_back(1);
       }
   }

   auto padding_width = block_sizes[0];
   for(uint64_t i = 1; i < block_sizes.size(); i++) {
       padding_width = std::gcd(padding_width, block_sizes[i]);
   }

    return padding_width;
}

uint64_t KSFunctionEvaluator::EstimatePadding(std::vector<uint64_t> &poly) const {
    return EstimatePadding(poly.data());
}

const std::shared_ptr<const KSFunctionEvaluationContext> KSFunctionEvaluator::GetContext() const {
    return m_context;
}

uint64_t KSFunctionEvaluator::GetEstimatePadding() const {
    return m_estimate_padding;
}

uint64_t KSFunctionEvaluator::GetPaddingWidth() const {
    return m_padding_width;
}


void KSFunctionEvaluator::EvalFunc(std::vector<uint64_t> &result, const std::vector<uint64_t> &input,
                                   std::vector<uint64_t> &lut_rlwe) {

}

void
KSFunctionEvaluator::EvalFunc(uint64_t *result, const uint64_t *input_lwe, uint64_t *lut_rlwe) {

}

void KSFunctionEvaluator::EvalFunc(std::vector<uint64_t> &result, const std::vector<uint64_t> &input,
                                   std::function<uint64_t(uint64_t)> &function) {

}

void
KSFunctionEvaluator::EvalFunc(uint64_t *result, const uint64_t *input_lwe, std::function<uint64_t(uint64_t)> &function) {

}
