//
// Created by leonard on 7/16/26.
//

#include "operators/sap_decomposition.h"
#include "operators/lwe_to_rgsw_conversion.h"
#include "operators/cggi_blind_rotation.h"
#include "operators/endo_glwe_conversion.h"
#include "static/modulus_switching.h"

#include <gtest/gtest.h>
#include <random>

class SAPDecompositionTests : public testing::Test {

protected:

    std::vector<GenericKey> test_keys;

    std::shared_ptr<SAPDecompositionContext> sap_context;

    void SetUp() override {

        // declare params
        uint64_t Q = -(1ull << 32) + 1;
        uint32_t N = 1 << 11;
        uint64_t n = 100;
        double std = 0;
        uint64_t basebits = 16;
        uint64_t base = 1 << basebits;
        uint64_t digits = 4;

        auto params_bin = std::make_shared<CGGIBlindRotationContext>(KeyDistribution::BINARY, Q, N, n, base, digits, std);
        auto ntt = params_bin->GetNTT();
        auto trace_params = std::make_shared<TraceEvaluationContext>(BINARY, ntt, base, digits, std);
        auto lwe_conversion_params = std::make_shared<LWEConversionContext>(BINARY, Q, N, n, base, digits, std);

        std::shared_ptr<OperatorContext<BlindRotator>> params_rot = std::reinterpret_pointer_cast<OperatorContext<BlindRotator>>(params_bin);
        sap_context = std::make_shared<SAPDecompositionContext>(params_rot, trace_params, lwe_conversion_params, 32);

        std::vector<uint64_t> tmp_key(std::max<uint64_t>(N, n), 0);


        // gen sparse LWE key
        for(uint32_t i = 0; i < 32; i++) {
            tmp_key[i] = rand() % 2;
        }
        test_keys.emplace_back("LWE_KEY", tmp_key.data(), n);

        // gen RLWE key
        for (uint32_t i = 0; i< N; i++) {
            tmp_key[i] = rand() % 2;
        }
        test_keys.emplace_back("RLWE_KEY", tmp_key.data(), N);

        // ntt'ed RLWE key
        ntt->ForwardNTT(tmp_key.data(), tmp_key.data());
        test_keys.emplace_back("RLWE_KEY_NTT", tmp_key.data(), N);
    }

};

TEST_F(SAPDecompositionTests, TestHomTrunc) {

    auto decomp = sap_context->ConstructOperator(test_keys);
    auto ntt = sap_context->GetTraceContext()->GetNTT();

    auto N = ntt->GetDimension();
    auto Q = ntt->GetModulus();
    auto radix = sap_context->GetDefaultRadix();
    auto quot_space = N / radix;

    std::random_device random_device;
    std::mt19937 gen(random_device());

    std::uniform_int_distribution<uint64_t> hi_dist(0, quot_space);
    std::uniform_int_distribution<uint64_t> lo_dist(0, radix);

    auto truncated_part = hi_dist(gen);
    auto exponent = truncated_part * radix + lo_dist(gen);

    AlignedBuffer rlwe_sample(2 * N);
    AlignedBuffer rlwe_sample_out(2 * N, 0);
    rlwe_sample[N + exponent] = 1;
    ntt->ForwardNTT(rlwe_sample.data() + N, rlwe_sample.data() + N);

    decomp->HomTrunc(rlwe_sample_out.data(), rlwe_sample.data(), radix);

    auto ntt_key = test_keys[2];
    intel::hexl::EltwiseMultMod(rlwe_sample_out.data(), rlwe_sample_out.data(), ntt_key.GetKeyPtr(), N, Q, 1);
    intel::hexl::EltwiseSubMod(rlwe_sample_out.data() + N, rlwe_sample_out.data() + N, rlwe_sample_out.data(), N, Q);
    ntt->BackwardNTT(rlwe_sample_out.data(), rlwe_sample_out.data() + N);


    for(uint64_t i = 0; i < N; i++) {
        if (i != truncated_part) {
            EXPECT_EQ(rlwe_sample_out[i], 0);
        } else {
            EXPECT_EQ(rlwe_sample_out[i], 1);
        }
    }
}

TEST_F(SAPDecompositionTests, TestResetAccumulator) {

    auto decomp = sap_context->ConstructOperator(test_keys);
    auto ntt = sap_context->GetTraceContext()->GetNTT();

    auto N = ntt->GetDimension();
    auto Q = ntt->GetModulus();
    auto radix = sap_context->GetDefaultRadix();
    auto quot_space = N / radix;

    std::random_device random_device;
    std::mt19937 gen(random_device());

    std::uniform_int_distribution<uint64_t> hi_dist(0, quot_space - 1);
    std::uniform_int_distribution<uint64_t> lo_dist(0, radix - 1);

    auto truncated_part = hi_dist(gen);
    auto exponent = truncated_part * radix + lo_dist(gen);

    AlignedBuffer rlwe_sample(2 * N, 0);
    rlwe_sample[N + exponent] = Q / (2 * N);
    ntt->ForwardNTT(rlwe_sample.data() + N, rlwe_sample.data() + N);

    decomp->ResetAccumulator(rlwe_sample.data(), radix);

    auto ntt_key = test_keys[2];
    intel::hexl::EltwiseMultMod(rlwe_sample.data(), rlwe_sample.data(), ntt_key.GetKeyPtr(), N, Q, 1);
    intel::hexl::EltwiseSubMod(rlwe_sample.data() + N, rlwe_sample.data() + N, rlwe_sample.data(), N, Q);
    ntt->BackwardNTT(rlwe_sample.data(), rlwe_sample.data() + N);

    ModulusSwitch(rlwe_sample.data(), N, Q, 2 * N, ModulusSwitchType::ROUND);

    auto first_non_zero = N + 1;

    auto test = Q / (2 * N);
    ModulusSwitch(&test, 1, Q, 2 * N, ModulusSwitchType::ROUND);

    for(uint64_t i = 0; i < N; i++) {
        if (rlwe_sample[i] != 0) {
            first_non_zero = i / radix;
            std::cerr << i << " " << exponent << " " <<  (exponent > i ? (exponent - i) : (i - exponent)) <<  std::endl;
            break;
        }
    }

    EXPECT_EQ(truncated_part, first_non_zero);
}

TEST_F(SAPDecompositionTests, TestDecomp) {

    auto decomp = sap_context->ConstructOperator(test_keys);
    auto ntt = sap_context->GetTraceContext()->GetNTT();

    auto N = ntt->GetDimension();
    auto Q = ntt->GetModulus();
    auto radix = sap_context->GetDefaultRadix();
    auto quot_space = N / radix;

    auto lwe_Q= radix;
    // ugly
    auto lwe_n = decomp->GetContext()->GetLWEConversionContext()->GetTargetDimension();
    for(uint64_t i = 0; i < 3; i++) {
        lwe_Q = lwe_Q * lwe_Q;
    }
    auto encryptor = LWEEncryptor(lwe_Q,lwe_n,0.0);
    auto decryptor = LWEEncryptor(Q, N, 0.0);

    std::vector<uint64_t> ct(lwe_n + 1);
    std::vector<uint64_t> cts_out((N + 1) * 8);
    encryptor.MakeLWE(ct.data(), (1ull << 10) | (1ull << 5) | 1, test_keys[0].GetKeyPtr());

    decomp->Decompose(cts_out.data(), ct.data(), radix);

    std::vector<uint64_t> phase(2 * N, 0);
    for(uint64_t i = 0; i < 8; i++) {

        auto res = decryptor.PhaseLWE(cts_out.data() + i * (N + 1), test_keys[1].GetKeyPtr());
        std::cerr << res << std::endl;

    }


}