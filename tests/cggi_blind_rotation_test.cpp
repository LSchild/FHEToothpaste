//
// Created by leonard on 4/8/26.
//

#include <gtest/gtest.h>
#include "base_crypto.h"
#include "operators/mux_operator.h"
#include "operators/cggi_blind_rotation.h"

class CGGIBlindRotTestGroup : public testing::Test {

protected:

    std::shared_ptr<CGGIBlindRotationContext> ctx_bin;
    std::shared_ptr<CGGIBlindRotationContext> ctx_ter;

    std::shared_ptr<CGGIBlindRotator> rotatorBinary;
    std::shared_ptr<CGGIBlindRotator> rotatorTernary;


    AlignedBuffer rlwe_secret;
    AlignedBuffer rlwe_secret_ntt;
    AlignedBuffer lwe_secret_binary;
    AlignedBuffer lwe_secret_ternary;

    void SetUp() override {

        uint64_t Q = -(1ull << 32) + 1;
        uint32_t N = 1 << 11;
        uint64_t n = 4096;
        double std = 0;
        uint64_t basebits = 11;
        uint64_t base = 1 << basebits;
        uint64_t digits = 6;

        ctx_bin = std::make_shared<CGGIBlindRotationContext>(KeyDistribution::BINARY, Q, N, n, base, digits, std);
        ctx_ter = std::make_shared<CGGIBlindRotationContext>(KeyDistribution::TERNARY, Q, N, n, base, digits, std);

        rlwe_secret = AlignedBuffer(N);
        rlwe_secret_ntt = AlignedBuffer(N);
        lwe_secret_binary = AlignedBuffer(n);
        lwe_secret_ternary = AlignedBuffer(n);

        std::srand(time(nullptr));

        for (uint32_t i = 0; i< N; i++) {
            rlwe_secret[i] = rand() % Q;
        }

        for(uint32_t i = 0; i < n; i++) {
            lwe_secret_binary[i] = rand() % 2;
            lwe_secret_ternary[i] = (2 * N + (rand() % 3) - 1) % (2 * N);
        }

        std::vector<GenericKey> bundle_binary = {{"LOL",  lwe_secret_binary.data(), lwe_secret_binary.size()},
                                                 {"LOL2", rlwe_secret.data(),       rlwe_secret.size()}};

        std::vector<GenericKey> bundle_ternary = {{"LOL",  lwe_secret_ternary.data(), lwe_secret_binary.size()},
                                                  {"LOL2", rlwe_secret.data(),        rlwe_secret.size()}};


        rotatorBinary = ctx_bin->ConstructOperator(bundle_binary);
        rotatorTernary = ctx_ter->ConstructOperator(bundle_ternary);

    }

};

TEST_F(CGGIBlindRotTestGroup, TestBinaryBlindRotate) {

    auto encryptor =  rotatorBinary->GetEncryptor();
    auto ntt= encryptor->GetNTT();
    auto Q = ctx_bin->GetModulus();
    auto N = ctx_bin->GetRingDimension();
    auto n = ctx_bin->GetLWEDimension();

    ntt->ForwardNTT(rlwe_secret_ntt.data(), rlwe_secret.data());

    AlignedBuffer data(2 * N);
    AlignedBuffer phase(N);
    AlignedBuffer lwe(n + 1);

    // create bogus lwe sample
    auto lwe_enc = LWEEncryptor(2 * N, n, 0.0);
    uint64_t msg = N + (rand() % (N));
    lwe_enc.MakeLWE(lwe.data(), msg, lwe_secret_binary.data());
    /*
    lwe[n] = msg;
    for(uint64_t  i = 0; i < n; i++) {
        lwe[i] = rand() % (2 * N);
        lwe[n] += lwe[i] * lwe_secret_binary[i];
    }
    lwe[n] %= (2 * N); */

    // Set up acc message
    for(uint64_t  i = 0; i < N; i++) {
        data[i] = i;
        data[i + N] = (Q - i) % Q;
    }

    // Set up acc
    AlignedBuffer rlwe_acc(2 * N);
    encryptor->MakeRLWE(rlwe_acc.data(), data.data(), rlwe_secret_ntt.data(), false);
    ntt->BackwardNTT(rlwe_acc.data(), rlwe_acc.data());
    ntt->BackwardNTT(rlwe_acc.data() + N, rlwe_acc.data() + N);



    AlignedBuffer result(2 * N, 0);
    auto start = std::chrono::high_resolution_clock::now();

    rotatorBinary->BlindRotate(result.data(), lwe.data(), rlwe_acc.data(), false);

    auto stop = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(stop-start).count();
    std::cerr << "Took " << elapsed << std::endl;

    encryptor->PhaseRLWE(phase.data(), result.data(), rlwe_secret_ntt.data());

    auto start_idx = (2 * N - msg) % (2 * N);
    for(uint64_t i = 0; i < N; i++) {
        EXPECT_EQ(phase[i], data[(i + start_idx) % (2 * N)]);
    }

}

TEST_F(CGGIBlindRotTestGroup, TestTernaryBlindRotate) {
    auto encryptor =  rotatorTernary->GetEncryptor();
    auto ntt= encryptor->GetNTT();
    auto Q = ctx_ter->GetModulus();
    auto N = ctx_ter->GetRingDimension();
    auto n = ctx_ter->GetLWEDimension();

    ntt->ForwardNTT(rlwe_secret_ntt.data(), rlwe_secret.data());

    AlignedBuffer data(2 * N);
    AlignedBuffer phase(N);
    AlignedBuffer lwe(n + 1);

    // create bogus lwe sample
    uint64_t msg = N + (rand() % (N));
    lwe[n] = msg;
    for(uint64_t  i = 0; i < n; i++) {
        lwe[i] = rand() % (2 * N);
        lwe[n] += (lwe[i] * lwe_secret_ternary[i]) % (2 * N);
    }
    lwe[n] %= (2 * N);

    // Set up acc message
    for(uint64_t  i = 0; i < N; i++) {
        data[i] = i;
        data[i + N] = (Q - i) % Q;
    }

    // Set up acc
    AlignedBuffer rlwe_acc(2 * N);
    encryptor->MakeRLWE(rlwe_acc.data(), data.data(), rlwe_secret_ntt.data(), false);
    ntt->BackwardNTT(rlwe_acc.data(), rlwe_acc.data());
    ntt->BackwardNTT(rlwe_acc.data() + N, rlwe_acc.data() + N);

    AlignedBuffer result(2 * N, 0);

    auto start = std::chrono::high_resolution_clock::now();
    rotatorTernary->BlindRotate(result.data(), lwe.data(), rlwe_acc.data(), false);
    auto stop = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(stop-start).count();
    std::cerr << "Took " << elapsed << std::endl;

    encryptor->PhaseRLWE(phase.data(), result.data(), rlwe_secret_ntt.data());

    auto start_idx = (2 * N - msg) % (2 * N);
    for(uint64_t i = 0; i < N; i++) {
        EXPECT_EQ(phase[i], data[(i + start_idx) % (2 * N)]);
    }
}