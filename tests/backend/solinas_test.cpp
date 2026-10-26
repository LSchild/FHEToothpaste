//
// Created by leonard on 8/19/26.
//

#include "backend/solinas.h"
#include <gtest/gtest.h>
#include <random>

class SolinasTestGroup : public testing::Test {

protected:

    enum TestOps {
        Add = 0, Sub = 1, Mul = 2
    };

    std::random_device m_random_device;
    std::mt19937 m_engine;
    std::uniform_int_distribution<uint64_t> m_mod_distr;
    std::uniform_int_distribution<uint64_t> m_op_distr;

    std::shared_ptr<SolinasBE::SolinasWorker> m_worker;

    void SetUp() override {
        m_engine = std::mt19937(m_random_device());
        m_mod_distr = std::uniform_int_distribution<uint64_t>(0, SolinasBE::SolinasWorker::kModulus - 1);
        m_op_distr = std::uniform_int_distribution<uint64_t>(0, 3 - 1);
        m_worker = std::make_shared<SolinasBE::SolinasWorker>(1u << 11);
    }

    uint64_t GenElem() {
        return m_mod_distr(m_engine);
    }

    uint64_t ComputeExpected(TestOps op, uint64_t op0, uint64_t op1) {
        auto Q = SolinasBE::SolinasWorker::kModulus;
        switch (op) {
            case Add: {
                __uint128_t sum = __uint128_t(op0) + __uint128_t(op1);
                if (sum >= Q) {
                    sum -= Q;
                }
                return sum;
            }
            case Sub: {
                if (op1 == 0) {
                    return op0;
                }
                auto inv = __uint128_t(Q) - op1;
                return ComputeExpected(Add, op0, inv);
            }
            case Mul: {
                auto prod = __uint128_t(op0) * __uint128_t(op1);
                uint64_t res = (prod % __uint128_t(Q));
                return res;
            }
        }
    }

};

TEST_F(SolinasTestGroup, TestAdd) {
    const auto rounds = 100000;
    for(uint64_t i = 0 ; i < rounds; i++) {
        auto op0 = GenElem();
        auto op1 = GenElem();

        auto expected = ComputeExpected(Add, op0, op1);
        auto computed = m_worker->AddMod(op0, op1);

        EXPECT_EQ(expected, computed);
    }
}

TEST_F(SolinasTestGroup, TestSub) {
    const auto rounds = 100000;
    for(uint64_t i = 0 ; i < rounds; i++) {
        auto op0 = GenElem();
        auto op1 = GenElem();

        auto expected = ComputeExpected(Sub, op0, op1);
        auto computed = m_worker->SubMod(op0, op1);

        EXPECT_EQ(expected, computed);
    }
}

TEST_F(SolinasTestGroup, TestMul) {
    const auto rounds = 100000;
    for(uint64_t i = 0 ; i < rounds; i++) {
        auto op0 = GenElem();
        auto op1 = GenElem();

        auto expected = ComputeExpected(Mul, op0, op1);
        auto computed = m_worker->MulMod(op0, op1);

        if (expected != computed) {
            std::cerr << op0 << " " << op1 << " " << expected << " " << computed << std::endl;
        }

        EXPECT_EQ(expected, computed);
    }
}

TEST_F(SolinasTestGroup, TestFMA) {
    const auto rounds = 100000;
    for(uint64_t i = 0 ; i < rounds; i++) {
        auto op0 = GenElem();
        auto op1 = GenElem();
        auto op2 = GenElem();

        auto expected_prod = ComputeExpected(Mul, op0, op1);
        auto expected = ComputeExpected(Add, expected_prod, op2);
        auto computed = m_worker->FMAMod(op0, op1, op2);

        if (expected != computed) {
            std::cerr << op0 << " " << op1 << " " << op2 << " " << expected << " " << computed << std::endl;
        }

        EXPECT_EQ(expected, computed);
    }
}

// if inversion works, so does exponentiation
TEST_F(SolinasTestGroup, TestInv) {
    const auto rounds = 100000;
    for(uint64_t i = 0 ; i < rounds; i++) {
        auto op0 = GenElem();
        auto computed = m_worker->InvMod(op0);
        auto expected_one = ComputeExpected(Mul, op0, computed);

        if (expected_one != 1) {
            std::cerr << op0 << " " << computed << std::endl;
        }

        EXPECT_EQ(expected_one, 1);
    }
}

TEST_F(SolinasTestGroup, TestNTTTrivial) {

    std::vector<uint64_t> vector_in(m_worker->GetDimension(), 0);
    std::vector<uint64_t> vector_out(m_worker->GetDimension(), 0);

    m_worker->ForwardNTT(vector_out, vector_in);

    for(auto& v: vector_out) {
        EXPECT_EQ(v, 0);
    }

}

TEST_F(SolinasTestGroup, TestForwardBackward) {

    std::vector<uint64_t> vector_ref(m_worker->GetDimension(), 0);

    AlignedBuffer vector_in(m_worker->GetDimension());

    for(uint64_t i = 0; i < vector_in.size(); i++) {
        vector_in[i] = GenElem();
        vector_ref[i] = vector_in[i];
    }

    AlignedBuffer vector_out(m_worker->GetDimension());

    m_worker->ForwardNTT(vector_in.data(), vector_in.data());
    m_worker->BackwardNTT(vector_in.data(), vector_in.data());

    for(uint64_t i = 0; i < vector_in.size(); i++) {
        EXPECT_EQ(vector_in[i], vector_ref[i]);
    }
}