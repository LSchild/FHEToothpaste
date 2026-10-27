//
// Goldilocks negacyclic MathWorker, AVX2 implementation.
//
// NTT layout (radix-2, roots of X^N + 1 merged into the twiddles):
//   Forward  = Cooley-Tukey, natural order in, bit-reversed order out.
//   Backward = Gentleman-Sande, bit-reversed order in, natural order out.
// A level with m blocks uses primitive 4m-th roots of unity as twiddles. For
// m <= 16 these are 64th roots of unity, i.e. powers of two (8 = 2^3 has order
// 64), so those 5 levels multiply with shifts instead of full products. On the
// inverse side the 1/N = 2^-k scaling is also a power of two and is merged into
// the last level, which therefore needs no general multiplication at all.
// The two levels with t = 2 and t = 1 (twiddle changes inside a vector) are
// fused into one pass that permutes data in registers.
//

#include "backend/solinas_worker.h"
#include "backend/solinas_internal.h"

#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>

namespace SolinasBE {

namespace {

using namespace avx2;

uint32_t BitReverse(uint32_t x, unsigned bits) {
    uint32_t r = 0;
    for (unsigned i = 0; i < bits; ++i) {
        r = (r << 1) | (x & 1);
        x >>= 1;
    }
    return r;
}

// ---- Level kernels --------------------------------------------------------

// One CT block: (u, v) -> (u + v*w, u - v*w). dsum/ddiff may alias up/vp.
template <class Mul>
inline void CtBlock(const uint64_t* up, const uint64_t* vp, uint64_t* dsum, uint64_t* ddiff,
                    size_t t, const Mul& mul) {
    for (size_t j = 0; j < t; j += 4) {
        const __m256i u = load(up + j);
        __m256i ts = mul.apply(load(vp + j));
        if constexpr (!Mul::kCanonical) ts = canon_s(ts);
        __m256i s, d;
        ct_butterfly(u, ts, s, d);
        store(dsum + j, s);
        store(ddiff + j, d);
    }
}

// One GS block: (u, v) -> (u + v, (u - v)*w). dsum/ddiff may alias up/vp.
template <class Mul>
inline void GsBlock(const uint64_t* up, const uint64_t* vp, uint64_t* dsum, uint64_t* ddiff,
                    size_t t, const Mul& mul) {
    for (size_t j = 0; j < t; j += 4) {
        __m256i s, d;
        gs_butterfly(load(up + j), load(vp + j), s, d);
        store(dsum + j, s);
        store(ddiff + j, flip(mul.apply(d)));
    }
}

// Forward levels t = 2 and t = 1, fused. Outputs are canonical.
//   w2: expanded twiddles of the t = 2 level, w1: twiddles of the t = 1 level.
void ForwardLastTwo(uint64_t* a, size_t n, const uint64_t* w2, const uint64_t* w1) {
    for (size_t c = 0; c < n / 8; ++c) {
        uint64_t* p = a + 8 * c;
        const __m256i x = load(p);      // a0 a1 a2 a3
        const __m256i y = load(p + 4);  // a4 a5 a6 a7
        // t = 2: pairs (a0 a1 | a2 a3) and (a4 a5 | a6 a7)
        __m256i u = _mm256_permute2x128_si256(x, y, 0x20);  // a0 a1 a4 a5
        __m256i v = _mm256_permute2x128_si256(x, y, 0x31);  // a2 a3 a6 a7
        {
            const __m256i w = load(w2 + 4 * c);
            const __m256i ts = canon_s(mul_s(v, w, _mm256_srli_epi64(w, 32)));
            ct_butterfly(u, ts, u, v);
        }
        // t = 1: pairs (a0 a1), (a2 a3), (a4 a5), (a6 a7)
        __m256i u1 = _mm256_unpacklo_epi64(u, v);  // a0 a2 a4 a6
        __m256i v1 = _mm256_unpackhi_epi64(u, v);  // a1 a3 a5 a7
        {
            const __m256i w = load(w1 + 4 * c);
            const __m256i ts = canon_s(mul_s(v1, w, _mm256_srli_epi64(w, 32)));
            ct_butterfly(u1, ts, u1, v1);
            u1 = canon(u1);
            v1 = canon(v1);
        }
        const __m256i lo = _mm256_unpacklo_epi64(u1, v1);  // a0 a1 a4 a5
        const __m256i hi = _mm256_unpackhi_epi64(u1, v1);  // a2 a3 a6 a7
        store(p, _mm256_permute2x128_si256(lo, hi, 0x20));
        store(p + 4, _mm256_permute2x128_si256(lo, hi, 0x31));
    }
}

// Inverse levels t = 1 and t = 2, fused. Reads in, writes out (may alias).
void InverseFirstTwo(uint64_t* out, const uint64_t* in, size_t n, const uint64_t* w1,
                     const uint64_t* w2) {
    for (size_t c = 0; c < n / 8; ++c) {
        const __m256i x = load(in + 8 * c);      // a0 a1 a2 a3
        const __m256i y = load(in + 8 * c + 4);  // a4 a5 a6 a7
        const __m256i a = _mm256_permute2x128_si256(x, y, 0x20);  // a0 a1 a4 a5
        const __m256i b = _mm256_permute2x128_si256(x, y, 0x31);  // a2 a3 a6 a7
        __m256i u = _mm256_unpacklo_epi64(a, b);                  // a0 a2 a4 a6
        __m256i v = _mm256_unpackhi_epi64(a, b);                  // a1 a3 a5 a7
        {
            const __m256i w = load(w1 + 4 * c);
            __m256i s, d;
            gs_butterfly(u, v, s, d);
            u = s;
            v = flip(mul_s(d, w, _mm256_srli_epi64(w, 32)));
        }
        __m256i u2 = _mm256_unpacklo_epi64(u, v);  // a0 a1 a4 a5
        __m256i v2 = _mm256_unpackhi_epi64(u, v);  // a2 a3 a6 a7
        {
            const __m256i w = load(w2 + 4 * c);
            __m256i s, d;
            gs_butterfly(u2, v2, s, d);
            u2 = s;
            v2 = flip(mul_s(d, w, _mm256_srli_epi64(w, 32)));
        }
        store(out + 8 * c, _mm256_permute2x128_si256(u2, v2, 0x20));
        store(out + 8 * c + 4, _mm256_permute2x128_si256(u2, v2, 0x31));
    }
}

// ---- Elementwise helpers --------------------------------------------------

template <class V, class S>
inline void Map2(uint64_t* out, const uint64_t* a, const uint64_t* b, uint64_t len, V vec, S sc) {
    uint64_t i = 0;
    for (; i + 4 <= len; i += 4) store(out + i, vec(load(a + i), load(b + i)));
    for (; i < len; ++i) out[i] = sc(a[i], b[i]);
}

template <class V, class S>
inline void Map3(uint64_t* out, const uint64_t* a, const uint64_t* b, const uint64_t* c,
                 uint64_t len, V vec, S sc) {
    uint64_t i = 0;
    for (; i + 4 <= len; i += 4) store(out + i, vec(load(a + i), load(b + i), load(c + i)));
    for (; i < len; ++i) out[i] = sc(a[i], b[i], c[i]);
}

template <class V, class S>
inline void Map1(uint64_t* out, const uint64_t* a, uint64_t len, V vec, S sc) {
    uint64_t i = 0;
    for (; i + 4 <= len; i += 4) store(out + i, vec(load(a + i)));
    for (; i < len; ++i) out[i] = sc(a[i]);
}

const auto kVecAdd = [](__m256i a, __m256i b) { return add_canon(a, b); };
const auto kVecSub = [](__m256i a, __m256i b) { return sub_canon(a, b); };
const auto kVecMul = [](__m256i a, __m256i b) { return mul_canon(a, b, hi32(b)); };
const auto kVecFma = [](__m256i a, __m256i b, __m256i c) {
    return mul_add_canon(a, b, hi32(b), c);
};

}  // namespace

// ---------------------------------------------------------------------------
// Construction
// ---------------------------------------------------------------------------

SolinasWorker::SolinasWorker(uint64_t dimension) : m_N(dimension), m_log_N(0), m_psi(0) {
    if (dimension < 128 || dimension > 65536 || (dimension & (dimension - 1)) != 0) {
        throw std::invalid_argument("GoldilocksWorker: dimension must be a power of two in [2^7, 2^16], got " +
                                    std::to_string(dimension));
    }
    while ((uint64_t{1} << m_log_N) < m_N) ++m_log_N;

    // 7 generates the multiplicative group, so psi has order exactly 2N.
    m_psi = elemops::pow_mod(7, (kP - 1) >> (m_log_N + 1));
    if (elemops::pow_mod(m_psi, m_N) != kP - 1) throw std::logic_error("GoldilocksWorker: bad root of unity");
    const uint64_t psi_inv = elemops::inv_mod(m_psi);

    std::vector<uint64_t> pw(m_N), pw_inv(m_N);
    pw[0] = pw_inv[0] = 1;
    for (uint64_t i = 1; i < m_N; ++i) {
        pw[i] = elemops::mul(pw[i - 1], m_psi);
        pw_inv[i] = elemops::mul(pw_inv[i - 1], psi_inv);
    }
    m_psi_rev.resize(m_N);
    m_psi_inv_rev.resize(m_N);
    for (uint64_t i = 0; i < m_N; ++i) {
        const uint32_t r = BitReverse(static_cast<uint32_t>(i), m_log_N);
        m_psi_rev[i] = pw[r];
        m_psi_inv_rev[i] = pw_inv[r];
    }

    // Level with N/4 blocks: vector c covers blocks 2c and 2c + 1.
    m_forward_twiddle_n4_rep2.resize(m_N / 2);
    m_backward_twiddle_n4_rep2.resize(m_N / 2);
    for (uint64_t c = 0; c < m_N / 8; ++c) {
        for (uint64_t lane = 0; lane < 4; ++lane) {
            m_forward_twiddle_n4_rep2[4 * c + lane] = m_psi_rev[m_N / 4 + 2 * c + lane / 2];
            m_backward_twiddle_n4_rep2[4 * c + lane] = m_psi_inv_rev[m_N / 4 + 2 * c + lane / 2];
        }
    }

    // Twiddles of the levels with <= 16 blocks are powers of two: find exponents.
    uint64_t pow2[192];
    pow2[0] = 1;
    for (int e = 1; e < 192; ++e) pow2[e] = elemops::add(pow2[e - 1], pow2[e - 1]);
    auto log2_of = [&](uint64_t w) -> unsigned {
        for (unsigned e = 0; e < 192; ++e)
            if (pow2[e] == w) return e;
        throw std::logic_error("GoldilocksWorker: twiddle is not a power of two");
    };
    for (unsigned k = 1; k < 32; ++k) {
        m_fwd_pow2_twiddle[k] = describe_pow2(log2_of(m_psi_rev[k]));
        m_backward_pow2_twiddle[k] = describe_pow2(log2_of(m_psi_inv_rev[k]));
    }
    // N^-1 = 2^-k = 2^(192 - k).
    m_inverse_dim = describe_pow2(192 - m_log_N);
    m_final_twiddle_scaled = describe_pow2((log2_of(m_psi_inv_rev[1]) + 192 - m_log_N) % 192);
    if (m_inverse_dim.kind != internals::Pow2Kind::kDiv || m_inverse_dim.negate)
        throw std::logic_error("GoldilocksWorker: unexpected scaling factor");
}

// ---------------------------------------------------------------------------
// NTT
// ---------------------------------------------------------------------------

void SolinasWorker::ForwardNTT(uint64_t* result, const uint64_t* input) {
    const size_t n = m_N;
    const uint64_t* src = input;
    size_t m = 1;
    size_t t = n / 2;

    // Levels m = 1..16: twiddles are +-2^e.
    for (; m <= 16; m <<= 1, t >>= 1) {
        for (size_t i = 0; i < m; ++i) {
            const internals::Pow2Twiddle& tw = m_fwd_pow2_twiddle[m + i];
            const size_t base = 2 * i * t;
            uint64_t* dsum = result + base;
            uint64_t* ddiff = result + base + t;
            // u + (-w)v = u - wv: swap where the two results go.
            if (tw.negate) std::swap(dsum, ddiff);
            with_pow2(tw, [&](const auto& mul) {
                CtBlock(src + base, src + base + t, dsum, ddiff, t, mul);
            });
        }
        src = result;
    }

    // Levels with t >= 4: general twiddle, constant per block.
    for (; t >= 4; m <<= 1, t >>= 1) {
        for (size_t i = 0; i < m; ++i) {
            const size_t base = 2 * i * t;
            uint64_t* p = result + base;
            CtBlock(p, p + t, p, p + t, t, MulGeneral(m_psi_rev[m + i]));
        }
    }

    // Levels t = 2 and t = 1.
    ForwardLastTwo(result, n, m_forward_twiddle_n4_rep2.data(), m_psi_rev.data() + n / 2);
}

void SolinasWorker::BackwardNTT(uint64_t* result, const uint64_t* input) {
    const size_t n = m_N;

    // Levels t = 1 and t = 2.
    InverseFirstTwo(result, input, n, m_psi_inv_rev.data() + n / 2, m_backward_twiddle_n4_rep2.data());

    size_t m = n / 8;
    size_t t = 4;

    // Levels with m >= 32: general twiddle, constant per block.
    for (; m >= 32; m >>= 1, t <<= 1) {
        for (size_t i = 0; i < m; ++i) {
            uint64_t* p = result + 2 * i * t;
            GsBlock(p, p + t, p, p + t, t, MulGeneral(m_psi_inv_rev[m + i]));
        }
    }

    // Levels m = 16..2: twiddles are +-2^e.
    for (; m >= 2; m >>= 1, t <<= 1) {
        for (size_t i = 0; i < m; ++i) {
            const internals::Pow2Twiddle& tw = m_backward_pow2_twiddle[m + i];
            uint64_t* p = result + 2 * i * t;
            const uint64_t* up = p;
            const uint64_t* vp = p + t;
            // (u - v)(-w) = (v - u)w: swap the operands (u + v is symmetric).
            if (tw.negate) std::swap(up, vp);
            with_pow2(tw, [&](const auto& mul) { GsBlock(up, vp, p, p + t, t, mul); });
        }
    }

    // Last level (m = 1) with the 1/N scaling merged in:
    //   a[j] = (u + v) / N,   a[j + N/2] = (u - v) * w / N.
    // Both factors are powers of two; outputs are canonical.
    {
        const uint64_t* up = result;
        const uint64_t* vp = result + t;
        if (m_final_twiddle_scaled.negate) std::swap(up, vp);
        const MulDiv scale(m_inverse_dim.amount);
        with_pow2(m_final_twiddle_scaled, [&](const auto& mul) {
            using M = std::decay_t<decltype(mul)>;
            for (size_t j = 0; j < t; j += 4) {
                __m256i s, d;
                gs_butterfly(load(up + j), load(vp + j), s, d);
                __m256i r = mul.apply(d);
                if constexpr (!M::kCanonical) r = canon_s(r);
                store(result + j, flip(scale.apply(s)));
                store(result + t + j, flip(r));
            }
        });
    }
}

void SolinasWorker::ForwardNTT(std::vector<uint64_t>& result, std::vector<uint64_t>& input) {
    if (input.size() != m_N) throw std::invalid_argument("GoldilocksWorker::ForwardNTT: input.size() != N");
    if (result.size() != m_N) result.resize(m_N);
    ForwardNTT(result.data(), input.data());
}

void SolinasWorker::BackwardNTT(std::vector<uint64_t>& result, std::vector<uint64_t>& input) {
    if (input.size() != m_N) throw std::invalid_argument("GoldilocksWorker::BackwardNTT: input.size() != N");
    if (result.size() != m_N) result.resize(m_N);
    BackwardNTT(result.data(), input.data());
}

// ---------------------------------------------------------------------------
// Scalar operations
// ---------------------------------------------------------------------------

uint64_t SolinasWorker::AddMod(uint64_t left, uint64_t right) { return elemops::add(left, right); }
uint64_t SolinasWorker::SubMod(uint64_t left, uint64_t right) { return elemops::sub(left, right); }
uint64_t SolinasWorker::MulMod(uint64_t left, uint64_t right) { return elemops::mul(left, right); }
uint64_t SolinasWorker::FMAMod(uint64_t left, uint64_t right, uint64_t summand) {
    return elemops::fma(left, right, summand);
}
uint64_t SolinasWorker::InvMod(uint64_t left) { return elemops::inv_mod(left); }
uint64_t SolinasWorker::PowMod(uint64_t base, uint64_t exponent) {
    return elemops::pow_mod(base, exponent);
}

// ---------------------------------------------------------------------------
// Elementwise operations
// ---------------------------------------------------------------------------

void SolinasWorker::AddModElw(uint64_t* __restrict result, const uint64_t* __restrict left,
                              const uint64_t* __restrict right, uint64_t len) {
    Map2(result, left, right, len, kVecAdd, elemops::add);
}

void SolinasWorker::AddEqModElw(uint64_t* __restrict left_out, const uint64_t* __restrict right,
                                uint64_t len) {
    Map2(left_out, left_out, right, len, kVecAdd, elemops::add);
}

void SolinasWorker::SubModElw(uint64_t* __restrict result, const uint64_t* __restrict left,
                              const uint64_t* __restrict right, uint64_t len) {
    Map2(result, left, right, len, kVecSub, elemops::sub);
}

void SolinasWorker::SubEqModElw(uint64_t* __restrict left_out, const uint64_t* __restrict right,
                                uint64_t len) {
    Map2(left_out, left_out, right, len, kVecSub, elemops::sub);
}

void SolinasWorker::MulModElw(uint64_t* __restrict result, const uint64_t* __restrict left,
                              const uint64_t* __restrict right, uint64_t len) {
    Map2(result, left, right, len, kVecMul, elemops::mul);
}

void SolinasWorker::MulEqModElw(uint64_t* __restrict left_out, const uint64_t* __restrict right,
                                uint64_t len) {
    Map2(left_out, left_out, right, len, kVecMul, elemops::mul);
}

void SolinasWorker::FMAModElw(uint64_t* __restrict result, const uint64_t* __restrict left,
                              const uint64_t* __restrict right, const uint64_t* __restrict summand,
                              uint64_t len) {
    Map3(result, left, right, summand, len, kVecFma, elemops::fma);
}

void SolinasWorker::FMAEqModElw(uint64_t* __restrict left_out, const uint64_t* __restrict right,
                                const uint64_t* __restrict summand, uint64_t len) {
    Map3(left_out, left_out, right, summand, len, kVecFma, elemops::fma);
}

void SolinasWorker::MulSModElw(uint64_t* __restrict result, const uint64_t* __restrict left,
                               uint64_t right, uint64_t len) {
    const MulGeneral mul(right);
    Map1(result, left, len, [&](__m256i a) { return flip(canon_s(mul.apply(a))); },
         [&](uint64_t a) { return elemops::mul(a, right); });
}

void SolinasWorker::MulSEqModElw(uint64_t* __restrict left_out, uint64_t right, uint64_t len) {
    const MulGeneral mul(right);
    Map1(left_out, left_out, len, [&](__m256i a) { return flip(canon_s(mul.apply(a))); },
         [&](uint64_t a) { return elemops::mul(a, right); });
}

void SolinasWorker::ReduceModElw(uint64_t* out, const uint64_t* in, uint64_t len) {
    Map1(out, in, len, [](__m256i a) { return canon(a); }, elemops::make_canonical);
}

}  // namespace goldilocks
