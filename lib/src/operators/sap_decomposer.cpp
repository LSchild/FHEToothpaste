//
// Created by leonard on 6/3/26.
//

#include <utility>
#include <iostream>

#include "operators/sap_decomposition.h"
#include "utils/generic_utils.h"
#include "utils/math_utils.h"
#include "static/modulus_switching.h"
#include "static/sample_extraction.h"

SAPDecompositionContext::SAPDecompositionContext(std::shared_ptr<OperatorContext<BlindRotator>> blind_rotation_context,
std::shared_ptr<TraceEvaluationContext> trace_context,
        std::shared_ptr<LWEConversionContext> lwe_conversion_context,
uint64_t default_radix) :
m_rotation_context(std::move(blind_rotation_context)),
m_packing_context(std::move(trace_context)),
m_conversion_context(std::move(lwe_conversion_context)), m_default_radix(default_radix) {

    m_restart_iter = ComputeRestartIterationForRadix(default_radix, IMPLICIT_SK_L0);
}

std::shared_ptr<LWEConversionContext> SAPDecompositionContext::GetLWEConversionContext() const {
    return m_conversion_context;
}

std::shared_ptr<OperatorContext<BlindRotator>> SAPDecompositionContext::GetBlindRotationContext() const {
    return m_rotation_context;
}

std::shared_ptr<TraceEvaluationContext> SAPDecompositionContext::GetTraceContext() const {
    return m_packing_context;
}

Container SAPDecompositionContext::GetInputContainer() const {
    auto dummy = m_conversion_context->GetInputContainer();
    auto out = m_conversion_context->GetOutputContainer(dummy);
    auto out_lwe = std::dynamic_pointer_cast<LWEContainerImpl>(out);

    auto in_n = out_lwe->GetN();
    auto in_q = UINT64_MAX;

    return std::make_shared<LWEContainerImpl>(in_n, in_q, 0.0);
}

long double SAPDecompositionContext::ComputeOutputVariance(long double input_variance) const {

    auto N = m_packing_context->GetDimension();

    long double out_var = 0.0;
    long double cluster_var = 0.0;

    for(uint32_t i = 0; i < m_restart_iter; i++) {
        auto tmp_var = m_rotation_context->ComputeOutputVariance(out_var);
        cluster_var = tmp_var * N * std::pow(m_default_radix, 2.0) / 12.0;
        out_var = m_conversion_context->ComputeOutputVariance(cluster_var);
    }

    return cluster_var;
}

Container SAPDecompositionContext::GetOutputContainer(Container container) const {
    auto dummy_rlwe_in = m_rotation_context->GetInputContainer();
    auto dummy_rlwe_out = std::dynamic_pointer_cast<RLWEContainerImpl>(m_rotation_context->GetOutputContainer(dummy_rlwe_in));

    auto o_var = ComputeOutputVariance(0.0);
    return std::make_shared<LWEContainerImpl>(dummy_rlwe_out->GetN(), dummy_rlwe_out->GetQ(), o_var);
}

const uint64_t SAPDecompositionContext::GetDefaultRadix() const {
    return m_default_radix;
}

const uint64_t SAPDecompositionContext::GetRestartIteration() const {
    return m_restart_iter;
}

const uint64_t SAPDecompositionContext::ComputeRestartIterationForRadix(uint64_t radix, uint64_t sk_hamming) const {

    long double log_max_p_fail = -60;
    long double blind_rotation_var = 0;
    long double last_state_extraction_var = 0;
    long double current_var = 0;

    auto Q = m_packing_context->GetModulus();
    auto N = m_packing_context->GetDimension();

    auto iter = 0;
    auto worst_trunc_factor = radix;
    auto state_elem_bound = N / radix;

    long double trigger_var;
    uint64_t P;

    while (true) {

        // compute iteration variance
        blind_rotation_var = m_rotation_context->ComputeOutputVariance(current_var);
        auto trunc_addition_var = worst_trunc_factor * blind_rotation_var;
        current_var = m_packing_context->ComputeOutputVariance(trunc_addition_var);

        if (current_var == 0.0) {
            return UINT64_MAX;
        }

        auto factor = std::max(radix, state_elem_bound);
        auto extraction_var = (blind_rotation_var * N * factor * factor) / 12.0;

        auto p_fail = EstimateFailureProbability(extraction_var, Q, factor);

        // we apply log next so we have to make sure we're within range
        // if we're at or below 0, we're fine
        if (p_fail <= 0.0) {
            iter++;
            continue;
        }

        auto log_p_fail= std::log2l(p_fail);
        if (log_p_fail <= log_max_p_fail) {
            iter++;
            last_state_extraction_var = extraction_var;
            continue;
        } else {
            break;
        }
    }

    auto lwe_conversion_var = m_conversion_context->ComputeOutputVariance(last_state_extraction_var);
    auto modulus_switch_var = EstimateModulusSwitchingVariance(lwe_conversion_var, Q, 2 * N, BINARY, sk_hamming);
    auto lwe_p_fail = EstimateFailureProbability(modulus_switch_var, 2 * N, state_elem_bound);
    if (lwe_p_fail <= 0.0)
        return iter;

    auto log_lwe_p_fail = std::log2(lwe_p_fail);
    if (log_lwe_p_fail <= log_max_p_fail) {
        return iter;
    } else {
        std::cerr << __FUNCTION__ << " Warning: LWE KeySwitch for restart induces large variance. Check your hamming weight & other params." << std::endl;
        std::cerr << __FUNCTION__ << " Warning: This is probably fine though." << std::endl;
        return iter - 1;
     }
}

OperatorID SAPDecompositionContext::GetOperatorID() const {
    return DECOMP_SAP;
}

void SAPDecompositionContext::SetDefaultRadix(uint64_t new_radix) {
    m_default_radix = new_radix;
    m_restart_iter = ComputeRestartIterationForRadix(new_radix, IMPLICIT_SK_L0);
}

void SAPDecompositionContext::SetBlindRotationContext(std::shared_ptr<OperatorContext<BlindRotator>> new_rot_context) {
    m_rotation_context = std::move(new_rot_context);
}

void SAPDecompositionContext::SetTraceContext(std::shared_ptr<TraceEvaluationContext> new_packing_context) {
    m_packing_context = new_packing_context;
}

void SAPDecompositionContext::SetLWEConversionContext(std::shared_ptr<LWEConversionContext> new_lwe_context) {
    m_conversion_context = new_lwe_context;
}

std::unique_ptr<SAPDecomposer> SAPDecompositionContext::ConstructOperator(const std::vector<GenericKey> &keys) const {

    std::vector<GenericKey> conv_keys;
    conv_keys.emplace_back(keys[1].GetKey());
    conv_keys.emplace_back(keys[0].GetKey());

    std::vector<GenericKey> pack_keys = {keys[1].GetKey()};

    auto lwe_l0 = 0;
    for(auto& v : keys[0].GetKey()) {
        if (v != 0)
            lwe_l0++;
    }

    auto rotator = m_rotation_context->ConstructOperator(keys);
    auto packer = m_packing_context->ConstructOperator(pack_keys);
    auto conv = m_conversion_context->ConstructOperator(conv_keys);

    return std::unique_ptr<SAPDecomposer>(new SAPDecomposer(shared_from_this(),std::move(rotator),std::move(packer),std::move(conv),m_default_radix, lwe_l0, m_restart_iter));
}



void SAPDecomposer::HomTrunc(uint64_t *output, uint64_t *input, uint64_t radix) {
    // check that radix is a power of 2
    assert((radix & (radix - 1)) == 0);
    auto trace_context = m_context->GetTraceContext();
    auto N = trace_context->GetDimension();
    auto Q = trace_context->GetModulus();

    auto log_N = IntLog2(N);
    auto log_radix = IntLog2(radix);
    auto b_rev_shift = (sizeof(N) * 8) - (log_N - log_radix);


    auto precon = intel::hexl::InverseMod(N, Q);

    auto packing_p = m_packing_buffer.data();
    auto poly_p = m_trunc_pad_poly.data();
    auto ntt = trace_context->GetNTT();
    ZERO_UINT64_ARR(packing_p, m_packing_buffer.size());
    ZERO_UINT64_ARR(output, 2 * N);

    AlignedBuffer automorphism_result(2 * N);
    if (m_last_radix != radix) {
        ZERO_UINT64_ARR(poly_p, N);
        // m_trunc_pad_poly = sum_{i = 0}^{radix - 1} X^{i - (radix - 1)}
        std::fill(poly_p + N - (radix - 1), poly_p + N, Q - 1);
        poly_p[0] = 1;
        ntt->ForwardNTT(poly_p, poly_p);
        m_last_radix = radix;
    }

    intel::hexl::EltwiseFMAMod(packing_p, input, precon, nullptr, 2 * N, Q, 1);

    // multiply by poly to move into correct coefs
    intel::hexl::EltwiseMultMod(packing_p, packing_p, poly_p, N, Q, 1);
    intel::hexl::EltwiseMultMod(packing_p + N, packing_p + N, poly_p, N, Q, 1);

    // partial trace nextlog_N - log_radix
    for(uint64_t i = 0; i < log_N - log_radix - 1; i++) {
        auto auto_index = (1 << (log_N - i)) + 1;
        ntt->BackwardNTT(automorphism_result.data(), packing_p);
        ntt->BackwardNTT(automorphism_result.data() + N, packing_p + N);

        m_packer->EvalAuto(packing_p + 2 * N, automorphism_result.data(), auto_index);
        intel::hexl::EltwiseAddMod(packing_p, packing_p, packing_p + 2 * N, 2 * N, Q);
        ZERO_UINT64_ARR(packing_p + 2 * N, 2 * N);
    }
    ZERO_UINT64_ARR(packing_p + 2 * N, m_packing_buffer.size() - 2 * N);


    // clean stack
    std::vector<std::pair<uint64_t, uint64_t>> stack;
    stack.reserve(2 * log_N);
    uint64_t stack_head = 1;
    stack.emplace_back(log_radix + 1, 0);

    // TODO optimize me away
    // TODO: monomial NTTs can be computed in O(N) instead of O(log(N) * N)
    AlignedBuffer shift_poly(N);

    while (stack_head != 0) {

        // pop elem from stack
        auto [level, exponent] = stack[--stack_head];
        auto current_elem = packing_p + stack_head * 2 * N;
        ZERO_UINT64_ARR(shift_poly.data(), N);

        if (level == 0) {
            auto exponent_brev = ReverseBitsGeneric(exponent) >> b_rev_shift;
            shift_poly[exponent_brev] = 1;
            ntt->ForwardNTT(shift_poly.data(), shift_poly.data());
            intel::hexl::EltwiseMultMod(current_elem, current_elem, shift_poly.data(), N, Q, 1);
            intel::hexl::EltwiseMultMod(current_elem + N, current_elem + N, shift_poly.data(), N, Q, 1);
            intel::hexl::EltwiseAddMod(output, current_elem, output, 2 * N, Q);
            ZERO_UINT64_ARR(current_elem, 2 * N);
            stack.pop_back();
            continue;
        } else {
            // unless for leaves, p_shift = X^{-N/2^level} = X^{2 * N - N/2^level} = -X^{N - N/2^{level}}
            auto shift_exponent_abs = N >> level;
            auto shift_exponent = N - shift_exponent_abs;
            shift_poly[shift_exponent] = Q - 1;
            ntt->ForwardNTT(shift_poly.data(), shift_poly.data());
        }

        ntt->BackwardNTT(current_elem + 2 * N, current_elem);
        ntt->BackwardNTT(current_elem + 2 * N + N, current_elem + N);
        // auto(current_elem)
        m_packer->EvalAuto(current_elem + 4 * N, current_elem + 2 * N, (1 << level) + 1);
        // current_elem - auto(current_elem)
        intel::hexl::EltwiseSubMod(current_elem + 2 * N, current_elem, current_elem + 4 * N, 2 * N, Q);
        // (current_elem - auto(current_elem)) * X^{-N/2^level}
        intel::hexl::EltwiseMultMod(current_elem + 2 * N, current_elem + 2 * N, shift_poly.data(), N, Q, 1);
        intel::hexl::EltwiseMultMod(current_elem + 2 * N + N, current_elem + 2 * N + N, shift_poly.data(), N, Q, 1);
        // current_elem + auto(current_elem)
        intel::hexl::EltwiseAddMod(current_elem, current_elem + 4 * N, current_elem, 2 * N, Q);
        ZERO_UINT64_ARR(current_elem + 4 * N, 2 * N);
        stack.pop_back();
        stack.emplace_back(level - 1, 2 * exponent);
        stack.emplace_back(level - 1, 2 * exponent + 1);
        stack_head += 2;
    }
}

void SAPDecomposer::ResetAccumulator(uint64_t *acc, uint64_t radix) {
    // resetting the accumulator consists of
    // - extracting the bits above \log_2(\alpha)
    // - key-switch
    // blind-rotate + truncate

    auto packing_context = m_packer->GetContext();
    auto Q = packing_context->GetModulus();
    auto N = packing_context->GetDimension();

    auto radix_log2 = IntLog2(radix);
    auto worker = packing_context->GetNTT();
    // ks mod
    auto Qks = m_lwe_converter->GetContext()->GetModulus();
    auto output_lwe_dimension = m_lwe_converter->GetContext()->GetTargetDimension();

    // define scratch space
    uint64_t* const m_buffer = m_packing_buffer.data();
    ZERO_UINT64_ARR(m_buffer, m_packing_buffer.size());

    uint64_t* const overflow_extraction_poly = m_buffer;
    uint64_t* const extracted_sample_N = overflow_extraction_poly + N;
    uint64_t* const extracted_sample_n = extracted_sample_N + 2 * N;
    uint64_t* const refreshed_accumulator = extracted_sample_n + 4 * N;


    // set up extraction poly
    overflow_extraction_poly[0] = 0;
    for(uint64_t i = 1; i < N; i++) {
        overflow_extraction_poly[i] = intel::hexl::SubUIntMod(0, ((N - i) >> radix_log2) << radix_log2, Q);
    }

    // perform extraction of current "overflow"
    worker->ForwardNTT(overflow_extraction_poly, overflow_extraction_poly);
    intel::hexl::EltwiseMultMod(acc, acc, overflow_extraction_poly, N, Q, 1);
    intel::hexl::EltwiseMultMod(acc + N, acc + N, overflow_extraction_poly, N, Q, 1);
    worker->BackwardNTT(acc, acc);
    worker->BackwardNTT(acc + N, acc + N);

    // sample extract
    SampleExtract(extracted_sample_N, acc, 0, N, Q);

    // optional mod switch
    if (Qks != Q) {
        ModulusSwitch(extracted_sample_N, N + 1, Q, Qks, ModulusSwitchType::ROUND);
    }

    m_lwe_converter->Convert(extracted_sample_n, extracted_sample_N);

    // mod switch to LWE Q
    ModulusSwitch(extracted_sample_n, output_lwe_dimension + 1, Qks, 2 * N, ModulusSwitchType::FLOOR);

    // input lwe now can be re-bootstrapped
    refreshed_accumulator[N] = Q / (2 * N);
    m_rotator->BlindRotate(refreshed_accumulator, extracted_sample_n, refreshed_accumulator);

    std::copy(refreshed_accumulator, refreshed_accumulator + 2 * N, acc);

}

void SAPDecomposer::Decompose(uint64_t *output, const uint64_t *const input, uint64_t radix) {

    // We extract expected input shapes
    auto br_input_container = m_rotator->GetContext()->GetInputContainer();
    auto br_input_tuple = std::dynamic_pointer_cast<TupleContainerImpl>(br_input_container);

    // LWE input
    auto lwe_container = std::dynamic_pointer_cast<LWEContainerImpl>(br_input_tuple->GetElem(0));
    auto lwe_dim_in = lwe_container->GetN();
    auto lwe_mod_in = lwe_container->GetQ();

    // RLWE input
    auto rlwe_container = std::dynamic_pointer_cast<RLWEContainerImpl>(br_input_tuple->GetElem(1));
    auto Q = rlwe_container->GetQ();
    auto N = rlwe_container->GetN();
    auto N_log2 = IntLog2(N);

    // Math worker
    auto worker = m_context->GetTraceContext()->GetNTT();

    // set up constants
    auto perturb_lo = radix * m_beta;
    auto perturb_hi = m_beta;
    auto radix_log2 = IntLog2(radix);
    // assumes radix is power of two
    auto radix_mask = radix - 1;

    // buffer for accumulator and extraction poly
    AlignedBuffer rlwe_scratch(5 * N, 0);
    auto accumulator_buffer = rlwe_scratch.data();
    // acc = Q / (2 * N) * 1
    // TODO: the scale can be set to Q / radix IF we do not restart
    accumulator_buffer[N] = Q / (2 * N);
    worker->ForwardNTT(accumulator_buffer + N, accumulator_buffer + N);
    // extraction poly encodes map x -> x mod radix, encoded in reverse
    auto extract_poly = accumulator_buffer + 2 * N;
    extract_poly[0] = 0;
    for(uint64_t i = 1; i < N; i++) {
        extract_poly[N - i] = Q - (i & radix_mask);
    }
    worker->ForwardNTT(extract_poly, extract_poly);
    // extraction output buffer
    auto extraction_buffer = accumulator_buffer + 3 * N;

    // buffers for current phase digit and updated input
    std::vector<uint64_t> lwe_scratch(2 * (lwe_dim_in + 1), 0);
    auto input_copy = lwe_scratch.data();
    auto current_sub_phase = lwe_scratch.data() + lwe_dim_in + 1;
    std::copy(input, input + lwe_dim_in + 1, input_copy);

    // dynamically determine input Q
    // assumes input Q is power of 2
    auto current_modulus = 1ull << IntLog2(input[0]);
    for(uint64_t i = 1; i <= lwe_dim_in; i++) {
        if (current_modulus <= input[i]) {
            current_modulus <<= 1;
        }
    }

    uint64_t restart_iteration = 0;
    if (radix == m_max_radix and m_beta <= SAPDecompositionContext::IMPLICIT_SK_L0) {
        restart_iteration = m_restart_iteration;
    } else {
        restart_iteration = m_context->ComputeRestartIterationForRadix(radix, m_beta);
    }

    auto current_output_digit = 0;

    // todo get rid of me
    AlignedBuffer tmp(2 * N);
    while (current_modulus > 1) {
        /* start by obtaining current phase digit */
        for(uint64_t i = 0; i < lwe_dim_in + 1; i++) {
            current_sub_phase[i] = input_copy[i] & radix_mask;
            input_copy[i] >>= radix_log2;
        }
        current_modulus >>= radix_log2;

        // perturb as required
        current_sub_phase[lwe_dim_in] = intel::hexl::AddUIntMod(current_sub_phase[lwe_dim_in], perturb_lo, lwe_mod_in);
        input_copy[lwe_dim_in] = intel::hexl::SubUIntMod(input_copy[lwe_dim_in], perturb_hi, current_modulus);

        worker->BackwardNTT(accumulator_buffer, accumulator_buffer);
        worker->BackwardNTT(accumulator_buffer + N, accumulator_buffer + N);
        m_rotator->BlindRotate(accumulator_buffer, current_sub_phase, accumulator_buffer);

        // multiply with extraction poly
        intel::hexl::EltwiseMultMod(extraction_buffer, accumulator_buffer, extract_poly, N, Q, 1);
        intel::hexl::EltwiseMultMod(extraction_buffer + N, accumulator_buffer + N, extract_poly, N, Q, 1);
        worker->BackwardNTT(extraction_buffer, extraction_buffer);
        worker->BackwardNTT(extraction_buffer + N, extraction_buffer + N);

        // actually extract
        auto current_digit_buf = output + current_output_digit * (N + 1);
        SampleExtract(current_digit_buf, extraction_buffer, 0, N, Q);
        std::cerr << current_digit_buf[N] << std::endl;


        if ((current_output_digit + 1) % restart_iteration == 0 and current_modulus != 0) {
            ResetAccumulator(accumulator_buffer, 0);
        }

        // truncate: accumulator_buffer: Ntt->Ntt
        worker->BackwardNTT(accumulator_buffer + N, accumulator_buffer + N);
        for(uint64_t  i = 0; i < N; i++) {
            std::cerr << accumulator_buffer[N + i] << " [" << i << "], ";
        }
        std::cerr << std::endl;
        worker->ForwardNTT(accumulator_buffer + N, accumulator_buffer + N);
        HomTrunc(tmp.data(), accumulator_buffer, radix);
        std::copy(tmp.begin(), tmp.end(), accumulator_buffer);
        worker->BackwardNTT(accumulator_buffer + N, accumulator_buffer + N);
        for(uint64_t  i = 0; i < N; i++) {
            std::cerr << accumulator_buffer[N + i] << " [" << i << "], ";
        }
        std::cerr << std::endl;
        worker->ForwardNTT(accumulator_buffer + N, accumulator_buffer + N);

        current_output_digit++;
    }
}

void SAPDecomposer::Decompose(std::vector<uint64_t> &output, const std::vector<uint64_t> &input, uint64_t radix) {
    Decompose(output.data(), input.data(), radix);
}

const std::shared_ptr<const SAPDecompositionContext> SAPDecomposer::GetContext() const {
    return m_context;
}

SAPDecomposer::SAPDecomposer(std::shared_ptr<const SAPDecompositionContext> ctx, std::unique_ptr<BlindRotator> rotator,
                             std::unique_ptr<TraceEvaluator> packer, std::unique_ptr<LWEtoLWEConverter> conv, uint64_t max_radix,
                             uint64_t lwe_sk_hamming_weight, uint64_t reset_period) : m_context(std::move(ctx)),
                             m_rotator(std::move(rotator)), m_packer(std::move(packer)), m_lwe_converter(std::move(conv)),
                             m_max_radix(max_radix), m_beta(lwe_sk_hamming_weight), m_restart_iteration(reset_period) {
    auto N = m_context->GetTraceContext()->GetDimension();
    auto log_N = IntLog2(N);
    m_packing_buffer.resize(2 * log_N * N);
    m_trunc_pad_poly.resize(N);
}