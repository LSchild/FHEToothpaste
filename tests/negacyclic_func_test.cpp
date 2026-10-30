//
// Created by Leonard on 10/30/26.
//
#include "operators/negacyclic_function_evaluator.h"

#include <gtest/gtest.h>
#include <random>

#include "operators/cggi_blind_rotation.h"

class NegacyclicFunctionTest : public testing::Test {

protected:

    std::vector<GenericKey> test_keys;

    std::shared_ptr<NegacyclicFunctionEvaluationContext> sap_context;
    std::unique_ptr<NegacyclicFunctionEvaluator> m_eval;

    void SetUp() override {

        // declare params
        uint64_t Q = -(1ull << 32) + 1;  uint64_t basebits = 16;
        //uint64_t Q = 36028797018972161; uint64_t basebits = 14;
        uint32_t N = 1 << 12;
        uint64_t n = 10;
        double std = 0;
        uint64_t base = 1 << basebits;
        uint64_t digits = 4;

        uint64_t auto_base = 1 << 16;
        uint64_t auto_digits = 4;

        auto params_bin = std::make_shared<CGGIBlindRotationContext>(KeyDistribution::BINARY, Q, N, n, base, digits, std);
        auto ntt = params_bin->GetNTT();

        auto params_rot = std::reinterpret_pointer_cast<OperatorContext<BlindRotator>>(params_bin);
        auto eval_context = std::make_shared<NegacyclicFunctionEvaluationContext>(params_rot);

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
        /*
        uint64_t toggle = 2;
        for(uint64_t i =0; i < N; i++) {
            tmp_key[i] = 0;// 1- (i & 1);
        }
        std::fill(tmp_key.begin(), tmp_key.begin() + toggle, 1);
*/
        test_keys.emplace_back("RLWE_KEY", tmp_key.data(), N);

        // ntt'ed RLWE key
        ntt->ForwardNTT(tmp_key.data(), tmp_key.data());
        test_keys.emplace_back("RLWE_KEY_NTT", tmp_key.data(), N);

        m_eval = eval_context->ConstructOperator(test_keys);
    }

};
