//
// Internal: scalar and AVX2 arithmetic modulo Q = 2^64 - 2^32 + 1.
//
// Useful identities (eps = 2^32 - 1):
//   2^64 = eps,   2^96 = -1,   2^192 = 1   (mod Q)
//   Q - 1 = 2^32 * eps, so Q - 1 = 0xFFFFFFFF00000000 = 2^64 - 2^32.
//
// AVX2 has no unsigned 64-bit compare. Where one is needed, values are kept in
// a "shifted" form x_s = x XOR 2^63, for which the signed compare gives the
// unsigned order. Functions whose name ends in _s take or return shifted values.
//
// Inside the NTT, values may be any 64-bit integer (not necessarily < Q);
// all functions below state which inputs must be canonical.
//

#ifndef LARGE_FUNCTIONS_GOLDILOCKS_ARITH_H
#define LARGE_FUNCTIONS_GOLDILOCKS_ARITH_H

#include "solinas_worker.h"

#include <cstdint>
#include <immintrin.h>

#if !defined(__AVX2__)
#error "The Goldilocks backend requires AVX2. Compile with -mavx2."
#endif

namespace SolinasBE {

    inline constexpr uint64_t kP = 0xFFFFFFFF00000001ULL;
    inline constexpr uint64_t kEps = 0x00000000FFFFFFFFULL;
    inline constexpr uint64_t kSign = 0x8000000000000000ULL;

    __extension__ using u128 = unsigned __int128;

// ---------------------------------------------------------------------------
// Scalar
// ---------------------------------------------------------------------------
    namespace elemops {

// Branch-free: the carries below are data dependent and would mispredict.
        inline uint64_t make_mux_mask(bool c) { return uint64_t{0} - static_cast<uint64_t>(c); }

        inline uint64_t make_canonical(uint64_t x) { return x - (kP & make_mux_mask(x >= kP)); }

// a, b canonical.
        inline uint64_t add(uint64_t a, uint64_t b) {
            uint64_t s;
            const bool carry = __builtin_add_overflow(a, b, &s);
            // On carry s + 2^64 - Q = s + eps; otherwise subtract Q if s >= Q. Both are s - Q mod 2^64.
            return s - (kP & make_mux_mask(carry | (s >= kP)));
        }

// a, b canonical.
        inline uint64_t sub(uint64_t a, uint64_t b) {
            uint64_t d;
            const bool borrow = __builtin_sub_overflow(a, b, &d);
            return d + (kP & make_mux_mask(borrow));
        }

// Any 128-bit x -> canonical.  x = lo + 2^64 (h0 + 2^32 h1) = lo + eps*h0 - h1.
        inline uint64_t reduce128(u128 x) {
            const uint64_t lo = static_cast<uint64_t>(x);
            const uint64_t hi = static_cast<uint64_t>(x >> 64);
            uint64_t t0;
            const bool borrow = __builtin_sub_overflow(lo, hi >> 32, &t0);
            t0 -= kEps & make_mux_mask(borrow);
            const uint64_t t1 = (hi << 32) - (hi & kEps);  // eps * h0
            uint64_t r;
            const bool carry = __builtin_add_overflow(t0, t1, &r);
            r += kEps & make_mux_mask(carry);
            return make_canonical(r);
        }

        inline uint64_t mul(uint64_t a, uint64_t b) { return reduce128(static_cast<u128>(a) * b); }

// a, b, c canonical (a*b + c < 2^128 then holds for any 64-bit a, b, c anyway).
        inline uint64_t fma(uint64_t a, uint64_t b, uint64_t c) {
            return reduce128(static_cast<u128>(a) * b + c);
        }

        inline uint64_t pow_mod(uint64_t base, uint64_t e) {
            uint64_t r = 1;
            base = make_canonical(base);
            while (e != 0) {
                if (e & 1) r = mul(r, base);
                base = mul(base, base);
                e >>= 1;
            }
            return r;
        }

        inline uint64_t square_n(uint64_t x, int n) {
            while (n-- > 0) x = mul(x, x);
            return x;
        }

// x^(Q-2) with a fixed addition chain (64 squarings, 9 multiplications).
// Q - 2 = (2^31 - 1) * 2^33 + (2^32 - 1). inv(0) = 0.
        inline uint64_t inv_mod(uint64_t x) {
            const uint64_t e1 = make_canonical(x);                     // x^(2^1 - 1)
            const uint64_t e2 = mul(square_n(e1, 1), e1);      // x^(2^2 - 1)
            const uint64_t e3 = mul(square_n(e2, 1), e1);      // x^(2^3 - 1)
            const uint64_t e6 = mul(square_n(e3, 3), e3);      // x^(2^6 - 1)
            const uint64_t e12 = mul(square_n(e6, 6), e6);     // x^(2^12 - 1)
            const uint64_t e24 = mul(square_n(e12, 12), e12);  // x^(2^24 - 1)
            const uint64_t e30 = mul(square_n(e24, 6), e6);    // x^(2^30 - 1)
            const uint64_t e31 = mul(square_n(e30, 1), e1);    // x^(2^31 - 1)
            const uint64_t e32 = mul(square_n(e31, 1), e1);    // x^(2^32 - 1)
            return mul(square_n(e31, 33), e32);
        }

    }  // namespace scalar

// ---------------------------------------------------------------------------
// AVX2, 4 lanes of 64 bit
// ---------------------------------------------------------------------------
    namespace avx2 {

        inline __m256i set1(uint64_t x) { return _mm256_set1_epi64x(static_cast<long long>(x)); }

        inline __m256i load(const uint64_t *p) {
            return _mm256_loadu_si256(reinterpret_cast<const __m256i *>(p));
        }

        inline void store(uint64_t *p, __m256i x) {
            _mm256_storeu_si256(reinterpret_cast<__m256i *>(p), x);
        }

// x <-> x_s
        inline __m256i flip(__m256i x) { return _mm256_xor_si256(x, set1(kSign)); }

// Copies the high 32 bits of each lane into the low 32 bits (runs on the shuffle
// port, leaving the shift/multiply ports free).
        inline __m256i hi32(__m256i x) {
            return _mm256_castps_si256(_mm256_movehdup_ps(_mm256_castsi256_ps(x)));
        }

// Any x_s -> canonical x_s (one conditional subtraction of Q suffices below 2^64).
        inline __m256i canon_s(__m256i xs) {
            const __m256i lt = _mm256_cmpgt_epi64(set1(kP ^ kSign), xs);  // x < Q
            return _mm256_add_epi64(xs, _mm256_andnot_si256(lt, set1(kEps)));
        }

        inline __m256i canon(__m256i x) {
            const __m256i lt = _mm256_cmpgt_epi64(set1(kP ^ kSign), flip(x));
            return _mm256_add_epi64(x, _mm256_andnot_si256(lt, set1(kEps)));
        }

// x_s + y (mod Q) for any x and y <= 2^64 - 2^32 (canonical values, or a
// product of two 32-bit values). Result is shifted, in [0, 2^64).
// The carry test only needs the high 32-bit halves: with that bound on y, a
// carry moves the high half by at least one, so a 32-bit compare (cheaper than
// the 64-bit one) gives the right answer in the high half of the mask.
        inline __m256i add_small_s(__m256i xs, __m256i y) {
            const __m256i r = _mm256_add_epi64(xs, y);
            const __m256i carry = _mm256_cmpgt_epi32(xs, r);
            return _mm256_add_epi64(r, _mm256_srli_epi64(carry, 32));  // +eps == -Q on wrap
        }

// x_s - y (mod Q) for any x and y <= 2^64 - 2^32. Result is shifted, in [0, 2^64).
        inline __m256i sub_small_s(__m256i xs, __m256i y) {
            const __m256i r = _mm256_sub_epi64(xs, y);
            const __m256i borrow = _mm256_cmpgt_epi32(r, xs);
            return _mm256_sub_epi64(r, _mm256_srli_epi64(borrow, 32));
        }

// Full 64x64 -> 128 bit product. y_hi must hold the high half of y in its low 32 bits.
        inline void mul64(__m256i x, __m256i y, __m256i y_hi, __m256i &hi, __m256i &lo) {
            const __m256i x_hi = hi32(x);
            const __m256i ll = _mm256_mul_epu32(x, y);
            const __m256i lh = _mm256_mul_epu32(x, y_hi);
            const __m256i hl = _mm256_mul_epu32(x_hi, y);
            const __m256i hh = _mm256_mul_epu32(x_hi, y_hi);
            const __m256i t0 = _mm256_add_epi64(hl, _mm256_srli_epi64(ll, 32));
            const __m256i t1 = _mm256_add_epi64(lh, _mm256_and_si256(t0, set1(kEps)));
            const __m256i t2 = _mm256_add_epi64(hh, _mm256_srli_epi64(t0, 32));
            hi = _mm256_add_epi64(t2, _mm256_srli_epi64(t1, 32));
            const __m256i t1_lo_up = _mm256_castps_si256(_mm256_moveldup_ps(_mm256_castsi256_ps(t1)));
            lo = _mm256_blend_epi32(ll, t1_lo_up, 0xAA);
        }

// Any 128-bit (hi, lo) -> shifted value in [0, 2^64).
        inline __m256i reduce128_s(__m256i hi, __m256i lo) {
            const __m256i t = sub_small_s(flip(lo), _mm256_srli_epi64(hi, 32));   // lo - h1
            return add_small_s(t, _mm256_mul_epu32(hi, set1(kEps)));              // + eps * h0
        }

// x * y for any x, y. Shifted result in [0, 2^64).
        inline __m256i mul_s(__m256i x, __m256i y, __m256i y_hi) {
            __m256i hi, lo;
            mul64(x, y, y_hi, hi, lo);
            return reduce128_s(hi, lo);
        }

// Canonical a, b -> canonical a + b.
        inline __m256i add_canon(__m256i a, __m256i b) {
            const __m256i nb = _mm256_sub_epi64(set1(kP), b);  // Q - b in [1, Q]
            const __m256i as = flip(a);
            const __m256i ds = _mm256_sub_epi64(as, nb);       // (a + b - Q), shifted
            const __m256i lt = _mm256_cmpgt_epi64(ds, as);     // borrow <=> a + b < Q
            return flip(_mm256_add_epi64(ds, _mm256_and_si256(lt, set1(kP))));
        }

// Canonical a, b -> canonical a - b.
        inline __m256i sub_canon(__m256i a, __m256i b) {
            const __m256i d = _mm256_sub_epi64(a, b);
            const __m256i lt = _mm256_cmpgt_epi64(flip(b), flip(a));
            return _mm256_add_epi64(d, _mm256_and_si256(lt, set1(kP)));
        }

// Any a, b -> canonical a * b.
        inline __m256i mul_canon(__m256i a, __m256i b, __m256i b_hi) {
            return flip(canon_s(mul_s(a, b, b_hi)));
        }

// Any a, b; canonical c -> canonical a * b + c.
        inline __m256i mul_add_canon(__m256i a, __m256i b, __m256i b_hi, __m256i c) {
            return flip(canon_s(add_small_s(mul_s(a, b, b_hi), c)));
        }

// ---------------------------------------------------------------------------
// Butterflies
// ---------------------------------------------------------------------------

// Cooley-Tukey: (u + t, u - t) for any u and canonical t (given shifted).
// Outputs are unshifted, in [0, 2^64).
        inline void ct_butterfly(__m256i u, __m256i ts, __m256i &sum, __m256i &diff) {
            const __m256i us = flip(u);
            __m256i s = _mm256_add_epi64(u, ts);                // (u + t) shifted
            const __m256i carry = _mm256_cmpgt_epi64(ts, s);    // u + t >= 2^64
            s = _mm256_add_epi64(s, _mm256_srli_epi64(carry, 32));
            __m256i d = _mm256_sub_epi64(us, ts);               // u - t, unshifted
            const __m256i borrow = _mm256_cmpgt_epi64(ts, us);  // t > u
            d = _mm256_sub_epi64(d, _mm256_srli_epi64(borrow, 32));
            sum = flip(s);
            diff = d;
        }

// Gentleman-Sande: (u + v, u - v) for any u, v. Outputs unshifted, in [0, 2^64).
        inline void gs_butterfly(__m256i u, __m256i v, __m256i &sum, __m256i &diff) {
            ct_butterfly(u, canon_s(flip(v)), sum, diff);
        }

// ---------------------------------------------------------------------------
// Multiplication by 2^e (Pow2Twiddle). apply() takes any x and returns x * 2^e
// shifted; kCanonical says whether that result is already < Q.
// ---------------------------------------------------------------------------

        inline __m128i shift_count(unsigned n) { return _mm_cvtsi32_si128(static_cast<int>(n)); }

        struct MulOne {
            static constexpr bool kCanonical = false;

            __m256i apply(__m256i x) const { return flip(x); }
        };

// e in [1, 32]: x * 2^e = lo + 2^64 * m with m < 2^32, so x * 2^e = lo + eps * m.
        struct MulLow {
            static constexpr bool kCanonical = false;
            __m128i left, right;

            explicit MulLow(unsigned e) : left(shift_count(e)), right(shift_count(64 - e)) {}

            __m256i apply(__m256i x) const {
                const __m256i lo = _mm256_sll_epi64(x, left);
                const __m256i m = _mm256_srl_epi64(x, right);
                return add_small_s(flip(lo), _mm256_mul_epu32(m, set1(kEps)));
            }
        };

// e in [33, 64]: 128-bit product formed with shifts, then the generic reduction.
        struct MulHigh {
            static constexpr bool kCanonical = false;
            __m128i left, right;

            explicit MulHigh(unsigned e) : left(shift_count(e)), right(shift_count(64 - e)) {}

            __m256i apply(__m256i x) const {
                // A shift count of 64 yields 0, which is what e = 64 needs.
                return reduce128_s(_mm256_srl_epi64(x, right), _mm256_sll_epi64(x, left));
            }
        };

// Write x = q * 2^r + rem with rem < 2^r, r in [1, 31]. Then
//   x / 2^r = q + rem * 2^-r = q - rem * 2^(96 - r) = q - c * eps,
// where c = rem * 2^(32 - r) < 2^32 is the low r bits of x moved to the top of a
// 32-bit word. Both q < 2^63 and c * eps < Q, so a single conditional
// correction yields a canonical result.
        struct MulDiv {  // x * 2^-r
            static constexpr bool kCanonical = true;
            __m128i r_count, l_count;

            explicit MulDiv(unsigned r) : r_count(shift_count(r)), l_count(shift_count(64 - r)) {}

            __m256i apply(__m256i x) const {
                const __m256i q = _mm256_srl_epi64(x, r_count);
                const __m256i c = _mm256_srli_epi64(_mm256_sll_epi64(x, l_count), 32);
                const __m256i ce = _mm256_mul_epu32(c, set1(kEps));
                return sub_small_s(flip(q), ce);  // q - c*eps: canonical in both branches
            }
        };

        struct MulNegDiv {  // x * 2^(96 - r) = -x * 2^-r = c * eps - q
            static constexpr bool kCanonical = true;
            __m128i r_count, l_count;

            explicit MulNegDiv(unsigned r) : r_count(shift_count(r)), l_count(shift_count(64 - r)) {}

            __m256i apply(__m256i x) const {
                const __m256i q = _mm256_srl_epi64(x, r_count);
                const __m256i c = _mm256_srli_epi64(_mm256_sll_epi64(x, l_count), 32);
                const __m256i ce = _mm256_mul_epu32(c, set1(kEps));
                return sub_small_s(flip(ce), q);  // q < 2^63 satisfies the bound
            }
        };

// Twiddle 2^e for e in [0, 192). A negative twiddle -2^e' is reported with
// negate = true; callers absorb the sign by swapping butterfly operands.
        inline internals::Pow2Twiddle describe_pow2(unsigned e) {
            using internals::Pow2Kind;
            internals::Pow2Twiddle d;
            if (e >= 161) {
                d.kind = Pow2Kind::kDiv;
                d.amount = static_cast<uint8_t>(192 - e);
                return d;
            }
            if (e >= 96) {
                d.negate = true;
                e -= 96;  // now in [0, 64]
            }
            if (e == 0) {
                d.kind = Pow2Kind::kOne;
            } else if (e <= 32) {
                d.kind = Pow2Kind::kLow;
                d.amount = static_cast<uint8_t>(e);
            } else if (e <= 64) {
                d.kind = Pow2Kind::kHigh;
                d.amount = static_cast<uint8_t>(e);
            } else {
                d.kind = Pow2Kind::kNegDiv;
                d.amount = static_cast<uint8_t>(96 - e);
            }
            return d;
        }

// Calls f(multiplier) with the multiplier type matching d (sign not included).
        template<class F>
        inline void with_pow2(const internals::Pow2Twiddle &d, F &&f) {
            using internals::Pow2Kind;
            switch (d.kind) {
                case Pow2Kind::kOne:
                    f(MulOne{});
                    break;
                case Pow2Kind::kLow:
                    f(MulLow(d.amount));
                    break;
                case Pow2Kind::kHigh:
                    f(MulHigh(d.amount));
                    break;
                case Pow2Kind::kNegDiv:
                    f(MulNegDiv(d.amount));
                    break;
                case Pow2Kind::kDiv:
                    f(MulDiv(d.amount));
                    break;
            }
        }

// Multiplication by an arbitrary twiddle w broadcast to all lanes.
        struct MulGeneral {
            static constexpr bool kCanonical = false;
            __m256i w, w_hi;

            explicit MulGeneral(uint64_t tw) : w(set1(tw)), w_hi(_mm256_srli_epi64(w, 32)) {}

            __m256i apply(__m256i x) const { return mul_s(x, w, w_hi); }
        };

    }  // namespace avx2
}  // namespace goldilocks

#endif  // LARGE_FUNCTIONS_GOLDILOCKS_ARITH_H
