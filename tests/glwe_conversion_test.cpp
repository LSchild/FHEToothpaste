#include <gtest/gtest.h>
#include "base_crypto.h"
#include "backend/aligned_vector.h"
#include "operators/endo_glwe_conversion.h"

//
// Created by leonard on 4/26/26.
//
class GLWEConversionTestGroup : public testing::Test {

protected:

    std::shared_ptr<RLWEEncryptor> enc_rlwe;
    std::shared_ptr<LWEEncryptor> enc_lwe_in;
    std::shared_ptr<LWEEncryptor> enc_lwe_out;
    std::shared_ptr<LWEtoLWEConverter> lwe_conv;
    std::shared_ptr<RLWEtoRLWEConverter> rlwe_conv;

    AlignedBuffer secret_rlwe_src, secret_rlwe_target, secret_rlwe_src_ntt, secret_rlwe_target_ntt;
    AlignedBuffer secret_lwe_src, secret_lwe_target;

    void SetUp() override {
        uint64_t Q = 36028797018972161; -(1ull << 32) + 1;
        uint32_t N = 1 << 11;
        uint64_t n_in = N;
        uint64_t n_out = 2048;
        uint64_t L = 1 << 14;
        uint64_t digits = 4;
        double std = 0;

        enc_rlwe = std::make_shared<RLWEEncryptor>(Q, N, std);
        enc_lwe_in = std::make_shared<LWEEncryptor>(Q, n_in, std);
        enc_lwe_out = std::make_shared<LWEEncryptor>(Q, n_out, std);

        secret_lwe_src.resize(n_in);
        secret_lwe_target.resize(n_out);
        secret_rlwe_src.resize(N);
        secret_rlwe_target.resize(N);
        secret_rlwe_src_ntt.resize(N);
        secret_rlwe_target_ntt.resize(N);


        std::srand(time(nullptr));

        for (uint32_t i = 0; i < N; i++) {
            secret_rlwe_src[i] = random() % 2;
            secret_rlwe_target[i] = random() % 2;
        }

        for(uint64_t i = 0; i < n_in; i++) {
            secret_lwe_src[i] = random() % 2;
        }

        for(uint64_t i = 0; i < n_out; i++) {
            secret_lwe_target[i] = random() % 2;
        }

        enc_rlwe->GetNTT()->ForwardNTT(secret_rlwe_target_ntt.data(), secret_rlwe_target.data());
        enc_rlwe->GetNTT()->ForwardNTT(secret_rlwe_src_ntt.data(), secret_rlwe_src.data());

        auto lwe_params = std::make_shared<LWEConversionContext>(BINARY, Q, n_in, n_out, L, digits, std);
        auto rlwe_params = std::make_shared<RLWEConversionContext>(BINARY, Q, N, L, digits, std);

        std::vector<GenericKey> rlwe_conv_keys = {
                {"RLWE_SOURCE", secret_rlwe_src.data(), N},
                {"RLWE_TARGET", secret_rlwe_target.data(), N}
        };

        std::vector<GenericKey> lwe_conv_keys = {
                {"LWE_SOURCE", secret_lwe_src.data(), secret_lwe_src.size()},
                {"LWE_TARGET", secret_lwe_target.data(), secret_lwe_target.size()}
        };

        rlwe_conv = rlwe_params->ConstructOperator(rlwe_conv_keys);

        lwe_conv = lwe_params->ConstructOperator(lwe_conv_keys);

    }

};

TEST_F(GLWEConversionTestGroup, TestLWEtoLWE) {

    auto Q = enc_lwe_in->GetModulus();
    auto t = 32;
    auto scal = Q / t;

    AlignedBuffer lwe_in(enc_lwe_in->GetDimension() + 1);
    AlignedBuffer lwe_out(enc_lwe_out->GetDimension() + 1);

    auto msg = 1 + (random() % (t - 1));
    auto encoded_msg = scal * msg;
    enc_lwe_in->MakeLWE(lwe_in.data(), encoded_msg, secret_lwe_src.data());
    auto start = std::chrono::high_resolution_clock::now();

    lwe_conv->Convert(lwe_out.data(), lwe_in.data());


    auto stop = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(stop-start).count();
    std::cerr << elapsed<< std::endl;

    auto val = enc_lwe_out->PhaseLWE(lwe_out.data(), secret_lwe_target.data());
    uint64_t corrected_val = std::round(double(t * val) / double(Q));

    EXPECT_EQ(corrected_val, msg);
}

TEST_F(GLWEConversionTestGroup, TestRLWEtoRLWE) {

    auto N = enc_rlwe->GetDimension();
    auto Q = enc_rlwe->GetDimension();

    AlignedBuffer rlwe_in(2 * N);
    AlignedBuffer rlwe_out(2 * N, 0);
    AlignedBuffer msg(N, 0);
    AlignedBuffer res(N, 0);

    enc_rlwe->MakeRLWE(rlwe_in.data(), msg.data(), secret_rlwe_src_ntt.data());
    enc_rlwe->GetNTT()->BackwardNTT(rlwe_in.data(), rlwe_in.data());

    rlwe_conv->Convert(rlwe_out.data(), rlwe_in.data());
    enc_rlwe->PhaseRLWE(res.data(), rlwe_out.data(), secret_rlwe_target_ntt.data());

    for(uint64_t i = 0; i < N; i++) {
        EXPECT_EQ(res[i], msg[i]);
    }

}