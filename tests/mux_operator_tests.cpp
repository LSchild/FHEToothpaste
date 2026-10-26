//
// Created by leonard on 3/24/26.
//

#include "base_crypto.h"
#include "operators/mux_operator.h"
#include "utils/generic_utils.h"
#include "gtest/gtest.h"

class MuxOperatorTestGroup : public testing::Test {

protected:

    std::shared_ptr<RLWEEncryptor> encryptor;
    std::shared_ptr<MuxOperator> mux_operator;
    AlignedBuffer secret;

    void SetUp() override {

        uint64_t Q = -(1ull << 32) + 1;
        uint32_t N = 1 << 11;
        double std = 0;
        uint64_t basebits = 16;
        uint64_t digits = 4;

        encryptor = std::make_shared<RLWEEncryptor>(Q, N, std);
        mux_operator = std::make_shared<MuxOperator>(encryptor->GetNTT(), basebits, digits);

        secret = AlignedBuffer(N);

        std::srand(time(nullptr));

        for (uint32_t i = 0; i< N; i++) {
            secret[i] = rand() % 2;
        }

        encryptor->GetNTT()->ForwardNTT(secret.data(), secret.data());

    }

};

TEST_F(MuxOperatorTestGroup, TestRLWEPrimeProduct) {
    auto N = encryptor->GetDimension();

    AlignedBuffer rgsw = AlignedBuffer(4 * N * mux_operator->GetRGSWDigits());
    AlignedBuffer poly = AlignedBuffer(N);

    AlignedBuffer data = AlignedBuffer(3 * N, 0);
    data[0] = 1;

    for (uint64_t i = 0; i < N; i++) {
        data[i + 2 * N] = i;
    }

    encryptor->MakeRGSW(rgsw.data(), data.data(), secret.data(), mux_operator->GetRGSWBasis(), mux_operator->GetRGSWDigits());
    std::fill(data.begin(), data.begin() + 2 * N, 0);
    mux_operator->RLWEPrimeProduct(data.data(), rgsw.data() + 2 * N * mux_operator->GetRGSWDigits(),
        data.data() + 2 * N, false);
    encryptor->PhaseRLWE(data.data() + 2 * N, data.data(), secret.data());

    for (uint32_t i = 0; i < N; i++) {
        EXPECT_EQ(data[i + 2 * N], i);
    }
    std::cerr << std::endl;
}

TEST_F(MuxOperatorTestGroup, TestExternalProduct) {
    auto N = encryptor->GetDimension();

    AlignedBuffer rgsw = AlignedBuffer(4 * N * mux_operator->GetRGSWDigits());
    AlignedBuffer rlwe0 = AlignedBuffer(2 * N);
    AlignedBuffer rlwe1 = AlignedBuffer(2 * N);

    AlignedBuffer data = AlignedBuffer(3 * N, 0);
    uint64_t bit = rand() % 2;
    data[0] = bit;

    for (uint64_t i = 0; i < N; i++) {
        data[i + N] = i;
        data[i + 2 * N] = i + N;
    }

    encryptor->MakeRGSW(rgsw.data(), data.data(), secret.data(), mux_operator->GetRGSWBasis(), mux_operator->GetRGSWDigits());
    encryptor->MakeRLWE(rlwe0.data(), data.data() + N, secret.data(), false);
    encryptor->MakeRLWE(rlwe1.data(), data.data() + 2 * N, secret.data(), false);

    mux_operator->BinaryMux(data.data(), rgsw.data(), rlwe0.data(), rlwe1.data());
    encryptor->PhaseRLWE(data.data() + 2 * N, data.data(), secret.data());

    for (uint32_t i = 0; i < N; i++) {
        EXPECT_EQ(data[2 * N + i], i + bit * N);
    }
}

TEST_F(MuxOperatorTestGroup, TestBinaryCmux) {
    auto N = encryptor->GetDimension();
    auto Q = encryptor->GetModulus();

    AlignedBuffer rgsw = AlignedBuffer(4 * N * mux_operator->GetRGSWDigits());
    AlignedBuffer rlwe0 = AlignedBuffer(2 * N);

    AlignedBuffer data = AlignedBuffer(3 * N, 0);
    uint64_t bit = rand() % 2;
    data[0] = bit;
    uint64_t delta = rand() % N;

    for (uint64_t i = 0; i < N; i++) {
        data[i + N] = i;
    }
    data[2 * N] = Q - 1;
    data[2 * N + delta] = (data[2 * N + delta] + 1) % Q;

    encryptor->MakeRGSW(rgsw.data(), data.data(), secret.data(), mux_operator->GetRGSWBasis(), mux_operator->GetRGSWDigits());
    encryptor->MakeRLWE(rlwe0.data(), data.data() + N, secret.data(), false);
    encryptor->GetNTT()->ForwardNTT(data.data() + 2 * N, data.data() + 2 * N);

    mux_operator->BinaryCMux(rlwe0.data(), rgsw.data(), data.data() + 2 * N);
    encryptor->PhaseRLWE(data.data() + 2 * N,   rlwe0.data(), secret.data());


    if (bit == 1) {
        for (uint64_t i = 0; i < delta; i++) {
            uint64_t e_val = Q - (N - delta + i);
            EXPECT_EQ(data[2 * N + i], e_val);
        }
        for (uint64_t i = delta; i < N; i++) {
            EXPECT_EQ(data[2 * N + i], i - delta);
        }
    } else {
        for (uint64_t i = 0; i < N; i++) {
            EXPECT_EQ(data[2 * N + i], i);
        }
    }
}

TEST_F(MuxOperatorTestGroup, TestTernaryCMUX) {
    auto N = encryptor->GetDimension();
    auto Q = encryptor->GetModulus();

    AlignedBuffer rgsw = AlignedBuffer(8 * N * mux_operator->GetRGSWDigits());
    AlignedBuffer rlwe0 = AlignedBuffer(2 * N);

    AlignedBuffer data = AlignedBuffer(3 * N, 0);
    uint64_t twit = rand() % 3;
    if (twit == 0) {
        data[0] = 0;
        encryptor->MakeRGSW(rgsw.data(), data.data(), secret.data(), mux_operator->GetRGSWBasis(), mux_operator->GetRGSWDigits());
        encryptor->MakeRGSW(rgsw.data() + 4 * N * mux_operator->GetRGSWDigits(), data.data(), secret.data(), mux_operator->GetRGSWBasis(), mux_operator->GetRGSWDigits());
    }
    if (twit == 1) {
        data[0] = 1;
        encryptor->MakeRGSW(rgsw.data(), data.data(), secret.data(), mux_operator->GetRGSWBasis(), mux_operator->GetRGSWDigits());
        data[0] = 0;
        encryptor->MakeRGSW(rgsw.data() + 4 * N * mux_operator->GetRGSWDigits(), data.data(), secret.data(), mux_operator->GetRGSWBasis(), mux_operator->GetRGSWDigits());
    }
    if (twit == 2) {
        data[0] = 0;
        encryptor->MakeRGSW(rgsw.data(), data.data(), secret.data(), mux_operator->GetRGSWBasis(), mux_operator->GetRGSWDigits());
        data[0] = 1;
        encryptor->MakeRGSW(rgsw.data() + 4 * N * mux_operator->GetRGSWDigits(), data.data(), secret.data(), mux_operator->GetRGSWBasis(), mux_operator->GetRGSWDigits());
    }
    for (uint64_t i = 0; i < N; i++) {
        data[i] = i;
    }
    encryptor->MakeRLWE(rlwe0.data(), data.data(), secret.data(), false);
    std::fill(data.begin(), data.end(), 0);

    uint64_t a = random() % (2 * N);
    uint64_t a_neg = (2 * N - a) % (2*N);
    if (a >= N) {
        uint64_t a_idx = a - N;
        data[a_idx] = Q - 1;
        data[N + a_neg] = 1;
    } else {
        uint64_t a_idx = a_neg - N;
        data[a] = 1;
        data[N + a_idx] = Q - 1;
    }
    data[0] = (data[0] + Q - 1) % Q;
    data[N] = (data[N] + Q - 1) % Q;

    encryptor->GetNTT()->ForwardNTT(data.data(), data.data());
    encryptor->GetNTT()->ForwardNTT(data.data() + N, data.data() + N);

    mux_operator->TernaryCMux(rlwe0.data(), rgsw.data(), data.data());
    std::fill(data.begin(), data.end(), 0);
    encryptor->PhaseRLWE(data.data(),   rlwe0.data(), secret.data());

    // Set up expected output
    for (uint64_t i = 0; i < N; i++) {
        data[i + N] = i;
        data[i + 2 * N] = (Q - i) % Q;
    }

    if (twit == 0) {
        for (uint64_t i = 0; i < N; i++) {
            EXPECT_EQ(data[i], i);
        }
    } else {
        if (twit == 1) {
            for (uint64_t i = 0; i < N; i++) {
                EXPECT_EQ(data[i], data[(i + a_neg) % (2 * N) + N]);
            }
        } else {
            for (uint64_t i = 0; i < N; i++) {
                EXPECT_EQ(data[i], data[(i + a) % (2 * N) + N]);
            }
        }
    }
}

TEST_F(MuxOperatorTestGroup, TestMultiMUX) {
    auto N = encryptor->GetDimension();
    auto Q = encryptor->GetModulus();
    auto L = mux_operator->GetRGSWBasis();
    auto l = mux_operator->GetRGSWDigits();
    auto rgsw_size = 4 * N * l;

    // random number of RGSW controls to generate
    auto k = random() % 16;

    AlignedBuffer weight_polys(k * N, 0);
    AlignedBuffer rgsw = AlignedBuffer(k * 4 * N * l);
    AlignedBuffer rlwe_in = AlignedBuffer(2 * N);
    AlignedBuffer data(N, 0);

    std::vector<uint64_t> bits(k);
    std::vector<uint64_t> weights(k);

    auto rgsw_pointer = rgsw.data();
    uint64_t expected_val = 0;
    uint64_t hm = 0;
    for(uint32_t i = 0; i < k; i++ ) {
        bits[i] = random() % 2;
        hm += bits[i];
        data[0] = bits[i];
        // rgsw
        encryptor->MakeRGSW(rgsw_pointer + rgsw_size * i, data.data(), secret.data(), L, l);
        // weight poly
        weights[i] = random() % N;
        auto w_idx = (N - weights[i]) % N;
        weight_polys[i * N + w_idx] = Q - 1;
        mux_operator->GetNTT()->ForwardNTT(weight_polys.data() + i * N, weight_polys.data() + i * N);
        expected_val += bits[i] * weights[i];
    }

    for (uint64_t i = 0; i < N; i++) {
        data[i] = i;
    }
    encryptor->MakeRLWE(rlwe_in.data(), data.data(), secret.data(), false);

    mux_operator->MultiMux(rlwe_in.data(), k, rgsw_pointer, weight_polys.data());

    // Set up expected output
    encryptor->PhaseRLWE(data.data(), rlwe_in.data(), secret.data());

    EXPECT_EQ(data[0], expected_val);


}