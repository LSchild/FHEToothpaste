//
// Created by leonard on 3/24/26.
//

#include "base_crypto.h"
#include "backend/aligned_vector.h"
#include "gtest/gtest.h"

class RLWEBaseTestGroup : public testing::Test {

protected:

    std::shared_ptr<RLWEEncryptor> enc_nonoise;
    std::shared_ptr<RLWEEncryptor> enc_noise;

    AlignedBuffer secret_nonoise;
    AlignedBuffer secret_noise;

    void SetUp() override {
        uint64_t Q = -(1ull << 32) + 1; 2251799813773313;
        uint32_t N = 1024;
        double std = 0;
        double std2 = 3.19;

        enc_nonoise = std::make_shared<RLWEEncryptor>(Q, N, std);
        enc_noise = std::make_shared<RLWEEncryptor>(Q, N, std2);

        secret_nonoise = AlignedBuffer(N);
        secret_noise = AlignedBuffer(N);

        std::srand(432534263);

        for (uint32_t i = 0; i< N; i++) {
            secret_nonoise[i] = rand() % 2;
            secret_noise[i] = rand() % 2;
        }

        enc_noise->GetNTT()->ForwardNTT(secret_noise.data(), secret_noise.data());
        enc_nonoise->GetNTT()->ForwardNTT(secret_nonoise.data(), secret_nonoise.data());
    }

};

TEST_F(RLWEBaseTestGroup, TestRLWEEncDec) {

    auto N = enc_noise->GetDimension();
    auto ntt_noise = enc_noise->GetNTT();
    auto ntt_nonoise = enc_nonoise->GetNTT();

    AlignedBuffer msg(N);
    AlignedBuffer msg_scale(N);

    AlignedBuffer rlwe(2 * N);
    AlignedBuffer rlwe_scale(2 * N);

    uint64_t scale = 1ull << 32;

    for (uint32_t i = 0; i < N; i++) {
        msg[i] = (i + 1);
        msg_scale[i] = scale * (i + 1);
    }

    enc_nonoise->MakeRLWE(rlwe.data(), msg.data(), secret_nonoise.data(), false);
    enc_noise->MakeRLWE(rlwe_scale.data(), msg_scale.data(), secret_noise.data(), false);

    enc_nonoise->PhaseRLWE(msg.data(), rlwe.data(), secret_nonoise.data());
    enc_noise->PhaseRLWE(msg_scale.data(), rlwe_scale.data(), secret_noise.data());

    for (uint32_t i = 0; i < N; i++) {
        EXPECT_EQ(msg[i], i+1);
        EXPECT_EQ((msg_scale[i] + 64) >> 32, i + 1);
    }
}

class LWEBaseTestGroup : public testing::Test {

protected:

    std::shared_ptr<LWEEncryptor> enc_nonoise;
    std::shared_ptr<LWEEncryptor> enc_noise;

    AlignedBuffer secret_nonoise;
    AlignedBuffer secret_noise;

    void SetUp() override {
        uint64_t Q = -(1ull << 32) + 1;// 2251799813773313;
        uint32_t n = 734;
        double std = 0;
        double std2 = 3.19;

        enc_nonoise = std::make_shared<LWEEncryptor>(Q, n, std);
        enc_noise = std::make_shared<LWEEncryptor>(Q, n, std2);

        secret_nonoise = AlignedBuffer(n);
        secret_noise = AlignedBuffer(n);

        std::srand(432534263);

        for (uint32_t i = 0; i< n; i++) {
            secret_nonoise[i] = 0;//rand() % 2;
            secret_noise[i] = rand() % 2;
        }

    }


};

TEST_F(LWEBaseTestGroup, TestLWEEncDecNoiseless) {

    auto dim = enc_nonoise->GetDimension();
    auto mod = enc_nonoise->GetModulus();

    auto value = random() % mod;
    AlignedBuffer ct(dim + 1);

    enc_nonoise->MakeLWE(ct.data(), value, secret_nonoise.data());

    auto res = enc_nonoise->PhaseLWE(ct.data(), secret_nonoise.data());

    EXPECT_EQ(res, value);
}

TEST_F(LWEBaseTestGroup, TestLWEEncDecNoise) {

    auto dim = enc_noise->GetDimension();
    auto mod = enc_noise->GetModulus();

    auto t = 128;

    auto value = random() % t;
    auto delta = mod / t;
    auto delta2 = mod / (2 * t);
    AlignedBuffer ct(dim + 1);

    enc_noise->MakeLWE(ct.data(), value * delta, secret_noise.data());

    auto res = enc_noise->PhaseLWE(ct.data(), secret_noise.data());
    auto res_decoded = (res + delta2) % mod;
    res_decoded = res_decoded / delta;

    EXPECT_EQ(res_decoded, value);
}