//
// Created by leonard on 4/27/26.
//
#include <cstdint>
#include <vector>
#include <algorithm>

#include "utils/math_utils.h"
#include "hexl/hexl.hpp"

uint64_t IntLog2(uint64_t input) {

    if (input <= 1) {
        return 0;
    }

    uint64_t res = 0;
    input -= 1;
    while (input > 0) {
        res++;
        input >>= 1;
    }

    return res;
}

std::vector<long double> EstimateSampleMean(uint64_t* seq, uint64_t n_samples, uint64_t len, uint64_t modulus) {
    std::vector<long double> means(len, 0.0);
    std::vector<long double> comp(len, 0.0);

    // numerically stable sum, do not optimize
    // see https://en.wikipedia.org/wiki/Kahan_summation_algorithm
    for(uint64_t i = 0; i < n_samples; i++) {
        for(uint64_t j = 0; j < len; j++) {
            long double v_ij = UnsignedToSignedRepr(seq[i * len + j], modulus);
            auto y = v_ij - comp[j];
            auto t = means[j] + y;
            comp[j] = (t - means[j]) - y;
            means[j] = t;
        }
    }

    std::transform(means.begin(), means.end(), means.begin(), [n_samples](long double x) {return x / ((long double) n_samples - 1);});
    return means;
}

std::vector<long double> EstimateSampleVariance(uint64_t* seq, uint64_t n_samples, uint64_t len, uint64_t modulus) {

    auto means = EstimateSampleMean(seq, n_samples, len, modulus);
    std::vector<long double> variances(len, 0.0);
    std::vector<long double> comp(len, 0.0);

    // numerically stable sum, do not optimize
    // see https://en.wikipedia.org/wiki/Kahan_summation_algorithm
    for(uint64_t i = 0; i < n_samples; i++) {
        for(uint64_t j = 0; j < len; j++) {
            long double v_ij = UnsignedToSignedRepr(seq[i * len + j], modulus);
            auto diff_square = (v_ij - means[j]) * (v_ij - means[j]);
            auto y = diff_square - comp[j];
            auto t = variances[j] + y;
            comp[j] = (t - variances[j]) - y;
            variances[j] = t;
        }
    }

    std::transform(variances.begin(), variances.end(), variances.begin(), [n_samples](long double x) {return x / ((long double) (n_samples - 1));});
    return variances;
}

void EvalNegacyclicAutomorphism(uint64_t* __restrict output, const uint64_t* __restrict input, uint64_t index, uint64_t n, uint64_t Q) {

    // todo in-place algo
    auto n_bits = IntLog2(n);
    auto n_mod_mask = n - 1;

    for(uint64_t buffer_idx = 0; buffer_idx < n; buffer_idx++) {
        auto permuted_index = buffer_idx * index;
        auto new_index = permuted_index & n_mod_mask;
        bool negate = ((permuted_index >> n_bits) & 1) == 1;
        auto read_val = input[buffer_idx];
        auto write_val = negate ? intel::hexl::SubUIntMod(0, read_val, Q) : read_val;
        output[new_index] = write_val;
    }

}

long double EstimateFailureProbability(long double variance, uint64_t modulus, uint64_t plaintext_modulus) {

    long double lQ = modulus;
    long double lP = plaintext_modulus;

    auto x = lQ / (2 * lP * std::sqrt(variance * 2));

    return std::erfcl(x);

}

