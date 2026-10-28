//
// Created by leonard on 6/1/26.
//

#include "operators/trace_evaluation.h"
#include "base_crypto.h"

#include <gtest/gtest.h>

class TraceTestGroup : public testing::Test {

protected:

    AlignedBuffer secret_nonoise, secret_ntt;
    std::shared_ptr<TraceEvaluationContext> m_auto_params;

    void SetUp() override {
        uint64_t Q = -(1ull << 32) + 1;
        uint32_t N = 1 << 12;
        uint32_t basis = 1 << 16;
        uint32_t digits = 4;
        double std = 0;

        m_auto_params = std::make_shared<TraceEvaluationContext>(BINARY, Q, N, basis, digits, std);
        secret_nonoise = AlignedBuffer(N);
        secret_ntt = AlignedBuffer(N);

        std::srand(4325243);

        for (uint32_t i = 0; i< N; i++) {
            secret_nonoise[i] = rand() % 2;
        }

        m_auto_params->GetNTT()->ForwardNTT(secret_ntt.data(),secret_nonoise.data());
    }

};

TEST_F(TraceTestGroup, TestAutomorphism) {

    std::vector<GenericKey> key = {
            {"RLWE", secret_nonoise.data(), secret_nonoise.size()}
    };

    auto eval = m_auto_params->ConstructOperator(key);
    auto N = m_auto_params->GetDimension();
    auto Q = m_auto_params->GetModulus();
    auto e_ntt = m_auto_params->GetNTT();
    auto encryptor = RLWEEncryptor(e_ntt, m_auto_params->GetStd());

    AlignedBuffer msg(N, 0);
    AlignedBuffer ct(2 * N, 0);
    AlignedBuffer ct_out(2 * N, 0);

    srand(time(nullptr));
    for(uint32_t i = 0; i < N; i++) {
        msg[i] = rand() % Q;
    }


    encryptor.MakeRLWE(ct.data(),msg.data(), secret_ntt.data(), false);
    e_ntt->BackwardNTT(ct.data(),ct.data());
    e_ntt->BackwardNTT(ct.data() + N,ct.data() + N);

    eval->Eval(ct_out.data(), ct.data());

    std::fill(ct.begin(),ct.end(), 0);
    encryptor.PhaseRLWE(ct.data(), ct_out.data(), secret_ntt.data());


    EXPECT_EQ(ct[0], msg[0]);
    for(uint32_t i = 1; i < N; i++) {
        EXPECT_EQ(ct[i], 0);
    }

}
