//
// Created by leonard on 10/27/26.
//

#include "operators/lmp_decomposition.h"
#include "static/modulus_switching.h"

LMPDecompositionContext::LMPDecompositionContext(std::shared_ptr<OperatorContext<BlindRotator>> blind_rotation_context,
                                                 std::shared_ptr<LWEConversionContext> lwe_conversion_context,
                                                 uint64_t default_radix) :
                                                 m_rotation_context(std::move(blind_rotation_context)),
                                                 m_conversion_context(std::move(lwe_conversion_context)),
                                                 m_default_radix(default_radix) {

}

std::shared_ptr<OperatorContext<BlindRotator>> LMPDecompositionContext::GetBlindRotationContext() const {
    return m_rotation_context;
}

std::shared_ptr<LWEConversionContext> LMPDecompositionContext::GetLWEConversionContext() const {
    return m_conversion_context;
}

OperatorID LMPDecompositionContext::GetOperatorID() const {
    return DECOMP_LMP;
}

const uint64_t LMPDecompositionContext::GetDefaultRadix() const {
    return m_default_radix;
}

Container LMPDecompositionContext::GetInputContainer() const {
    auto conv_in = m_conversion_context->GetInputContainer();
    auto conv_out = m_conversion_context->GetOutputContainer(conv_in);

    return conv_out;
}

Container LMPDecompositionContext::GetOutputContainer(Container container) const {
    // the converter switches from dim N to dim n, 2 * N, being the output modulus too :)
    auto q_out = 2 * m_conversion_context->GetSourceDimension();

    auto conv_in = m_conversion_context->GetInputContainer();
    auto conv_out = m_conversion_context->GetOutputContainer(conv_in);
    auto lwe_out = std::dynamic_pointer_cast<LWEContainerImpl>(conv_in);

    // todo key distr
    auto var_out = EstimateModulusSwitchingVariance(lwe_out->GetVariance(), lwe_out->GetQ(), q_out, BINARY, lwe_out->GetN() / 2);
    return std::make_shared<LWEContainerImpl>(lwe_out->GetN(), q_out, var_out);
}