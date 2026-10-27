//
// Created by leonard on 10/27/26.
//

#include "operators/lmp_decomposition.h"
#include "static/modulus_switching.h"
#include "utils/math_utils.h"

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

void LMPDecompositionContext::SetLWEConversionContext(std::shared_ptr<LWEConversionContext> new_lwe_context) {
    m_conversion_context = std::move(new_lwe_context);
}

void LMPDecompositionContext::SetBlindRotationContext(std::shared_ptr<OperatorContext<BlindRotator>> new_rot_context) {
    m_rotation_context = std::move(new_rot_context);
}

long double LMPDecompositionContext::ComputeOutputVariance(long double input_variance) const {
    auto input_trivial = GetInputContainer();
    auto output_trivial = GetOutputContainer(input_trivial);
    auto lwe_out = std::dynamic_pointer_cast<LWEContainerImpl>(output_trivial);

    auto var_out = std::max(input_variance, lwe_out->GetVariance());

    return var_out;
}

void LMPDecompositionContext::SetDefaultRadix(uint64_t new_radix) {
    m_default_radix = new_radix;
}

std::unique_ptr<LMPDecomposer> LMPDecompositionContext::ConstructOperator(const std::vector<GenericKey> &keys) const {
    auto N = m_conversion_context->GetSourceDimension();
    auto alpha_lmp = (2 * N) / m_default_radix;
    auto beta = alpha_lmp / 4;

    std::vector<GenericKey> conv_keys = {keys[1], keys[0]};

    auto rot = m_rotation_context->ConstructOperator(keys);
    auto conv = m_conversion_context->ConstructOperator(conv_keys);

    return std::unique_ptr<LMPDecomposer>(new LMPDecomposer(shared_from_this(), std::move(rot), std::move(conv), beta));
}

LMPDecomposer::LMPDecomposer(std::shared_ptr<const LMPDecompositionContext> ctx, std::unique_ptr<BlindRotator> rotator,
                             std::unique_ptr<LWEtoLWEConverter> converter, uint64_t beta) :
                             m_context(std::move(ctx)), m_rotator(std::move(rotator)), m_lwe_converter(std::move(converter)),
                             m_beta(beta) {

}

const std::shared_ptr<const LMPDecompositionContext> LMPDecomposer::GetContext() const {
    return m_context;
}

void LMPDecomposer::Decompose(std::vector<uint64_t> &output, const std::vector<uint64_t> &input, uint64_t radix) {
    Decompose(output.data(), input.data(), radix);
}

void LMPDecomposer::Decompose(uint64_t *output, const uint64_t *const input, uint64_t radix) {

    auto conversion_context = m_lwe_converter->GetContext();
    auto N = conversion_context->GetSourceDimension();
    auto n = conversion_context->GetTargetDimension();
    auto Q = conversion_context->GetModulus();

    auto alpha = (2 * N / radix);
    auto alpha_bits = IntLog2(alpha);
    if (alpha  <= 4 * m_beta) {
        throw std::invalid_argument("Input radix is not compatible with parameters !");
    }

    AlignedBuffer input_copy(n + 1, 0);
    AlignedBuffer bootstrap_input_buffer(n + 1, 0);
    AlignedBuffer bootstrap_output_buffer(n + 1, 0);


    /* dynamic input modulus estimation */
    uint64_t q_in = 1ull << IntLog2(input[0]);
    for(uint64_t i = 0; i < n + 1; i++) {
        input_copy[i] = input[i];
        if (input_copy[i] >= q_in)
            q_in *= 2;
    }
    auto current_modulus = q_in;
    auto current_modulus_mask = q_in - 1;

    auto current_digit_index = 0;

    /* get least significant digit */
    auto mask_2N = (2 * N) - 1;


    while (current_modulus > (2 * N)) {

        for(uint64_t i = 0; i < n + 1; i++) {
            output[current_digit_index * (n + 1) + i] = input_copy[i] & mask_2N;
        }
        current_digit_index += 1;

        // HomFloor
        input_copy[n] = (input_copy[n] + m_beta) & current_modulus_mask;
        for(uint64_t i = 0; i < n+1; i++) {
            bootstrap_input_buffer[i] = input_copy[i] & mask_2N;
        }

        EvalF0(bootstrap_output_buffer.data(), bootstrap_input_buffer.data());
        ModulusSwitch(bootstrap_output_buffer.data(), n + 1, Q, current_modulus, ModulusSwitchType::ROUND);

        for(uint64_t i = 0; i < n + 1; i++) {
            input_copy[i] = (input_copy[i] - bootstrap_output_buffer[i]) & current_modulus_mask;
        }
        input_copy[n] = (input_copy[n] + m_beta - (N / 2)) & current_modulus_mask;

        for(uint64_t i = 0; i < n+1; i++) {
            bootstrap_input_buffer[i] = input_copy[i] & mask_2N;
        }

        EvalF1(bootstrap_output_buffer.data(), bootstrap_input_buffer.data());
        ModulusSwitch(bootstrap_output_buffer.data(), n + 1, Q, current_modulus, ModulusSwitchType::ROUND);
        for(uint64_t i = 0; i < n + 1; i++) {
            input_copy[i] = (input_copy[i] - bootstrap_output_buffer[i]) & current_modulus_mask;
        }

        ModulusSwitch(input_copy.data(), n + 1, current_modulus, (current_modulus_mask >> alpha_bits), ModulusSwitchType::ROUND);

        current_modulus >>= alpha_bits;
        current_modulus_mask = current_modulus - 1;

    }



}

void LMPDecomposer::EvalF0(uint64_t *output, uint64_t *input) {

    auto conversion_context = m_lwe_converter->GetContext();
    auto N = conversion_context->GetSourceDimension();

    AlignedBuffer acc(2 * N);
    for(uint64_t i = 0; i < N; i++) {
        acc[i]
    }
}