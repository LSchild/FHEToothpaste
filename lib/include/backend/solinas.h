//
// MathWorker for the negacyclic ring Z_Q[X]/(X^N + 1) with the Goldilocks prime
// Q = 2^64 - 2^32 + 1 and N = 2^k (7 <= k <= 16; tuned for k in {10, 11, 12}).
//
// Contract:
//  - Inputs are expected canonical (in [0, Q)); every output is canonical.
//  - ForwardNTT maps coefficients (natural order) to evaluations at the roots of
//    X^N + 1 in bit-reversed order. The order is only meant for pointwise
//    operations followed by BackwardNTT.
//  - BackwardNTT includes the 1/N scaling: BackwardNTT(ForwardNTT(a)) == a.
//  - For both NTTs, result may equal input. The std::vector overloads require
//    input.size() == N and resize result to N if needed.
//  - InvMod(0) returns 0.
//  - All tables are read-only after construction, so one instance can be used
//    from several threads at once.
//  - Requires AVX2 (compile with -mavx2). Unaligned load/store instructions are
//    used throughout; they cost nothing extra on aligned data.
//

#ifndef LARGE_FUNCTIONS_GOLDILOCKS_WORKER_H
#define LARGE_FUNCTIONS_GOLDILOCKS_WORKER_H

#include "backend.h"
#include "aligned_vector.h"

#include <cstddef>
#include <cstdint>
#include <new>
#include <vector>

namespace SolinasBE {

namespace internals {

// Multiplication by a power of two 2^e mod Q (e in [0, 192), since 2^96 = -1 and 2^192 = 1).
// Which instruction sequence is used depends on e; see goldilocks_arith.h.
enum class Pow2Kind : uint8_t {
    kOne,     // e = 0                 : x
    kLow,     // e in [1, 32]          : x * 2^e, product fits in 96 bits
    kHigh,    // e in [33, 64]         : x * 2^e, product fits in 128 bits
    kNegDiv,  // e in [65, 95]         : x * 2^e = -x / 2^(96 - e)
    kDiv,     // e in [161, 191]       : x * 2^e =  x / 2^(192 - e)
};

struct Pow2Twiddle {
    Pow2Kind kind = Pow2Kind::kOne;
    uint8_t amount = 0;   // shift amount (kLow, kHigh) or divisor exponent (kNegDiv, kDiv)
    bool negate = false;  // the twiddle is -(2^e'); handled by swapping butterfly operands
};




}  // namespace detail

class SolinasWorker final : public MathWorker {
public:
    static constexpr uint64_t kModulus = 0xFFFFFFFF00000001ULL;

    // Throws std::invalid_argument unless dimension is a power of two in [2^7, 2^16].
    explicit SolinasWorker(uint64_t dimension);

    void ForwardNTT(uint64_t* result, const uint64_t* input) override;
    void ForwardNTT(std::vector<uint64_t>& result, std::vector<uint64_t>& input) override;
    void BackwardNTT(uint64_t* result, const uint64_t* input) override;
    void BackwardNTT(std::vector<uint64_t>& result, std::vector<uint64_t>& input) override;

    uint64_t AddMod(uint64_t left, uint64_t right) override;
    uint64_t SubMod(uint64_t left, uint64_t right) override;
    uint64_t MulMod(uint64_t left, uint64_t right) override;
    uint64_t FMAMod(uint64_t left, uint64_t right, uint64_t summand) override;
    uint64_t InvMod(uint64_t left) override;
    uint64_t PowMod(uint64_t base, uint64_t exponent) override;

    void AddModElw(uint64_t* __restrict result, const uint64_t* __restrict left,
                   const uint64_t* __restrict right, uint64_t len) override;
    void AddEqModElw(uint64_t* __restrict left_out, const uint64_t* __restrict right,
                     uint64_t len) override;
    void SubModElw(uint64_t* __restrict result, const uint64_t* __restrict left,
                   const uint64_t* __restrict right, uint64_t len) override;
    void SubEqModElw(uint64_t* __restrict left_out, const uint64_t* __restrict right,
                     uint64_t len) override;
    void MulModElw(uint64_t* __restrict result, const uint64_t* __restrict left,
                   const uint64_t* __restrict right, uint64_t len) override;
    void MulEqModElw(uint64_t* __restrict left_out, const uint64_t* __restrict right,
                     uint64_t len) override;
    void FMAModElw(uint64_t* __restrict result, const uint64_t* __restrict left,
                   const uint64_t* __restrict right, const uint64_t* __restrict summand,
                   uint64_t len) override;
    void FMAEqModElw(uint64_t* __restrict left_out, const uint64_t* __restrict right,
                     const uint64_t* __restrict summand, uint64_t len) override;
    void MulSModElw(uint64_t* __restrict result, const uint64_t* __restrict left, uint64_t right,
                    uint64_t len) override;
    void MulSEqModElw(uint64_t* __restrict left_out, uint64_t right, uint64_t len) override;
    void ReduceModElw(uint64_t* out, const uint64_t* in, uint64_t len) override;

    [[nodiscard]] uint64_t GetModulus() const override { return kModulus; }
    [[nodiscard]] uint64_t GetDimension() const override { return m_N; }

    // The primitive 2N-th root of unity psi used by the transform (psi = 7^((Q-1)/2N)).
    [[nodiscard]] uint64_t GetPsi() const { return m_psi; }

private:
    uint64_t m_N;
    unsigned m_log_N;
    uint64_t m_psi;

    // psi_rev_[k] = psi^brv(k), psi_inv_rev_[k] = psi^-brv(k) (brv over log_n_ bits).
    AlignedBuffer m_psi_rev;
    AlignedBuffer m_psi_inv_rev;
    // Twiddles of the level with N/4 blocks, each repeated twice so that one
    // 4-lane vector covers two blocks of 4 coefficients.
    AlignedBuffer m_forward_twiddle_n4_rep2;
    AlignedBuffer m_backward_twiddle_n4_rep2;

    // Levels with at most 16 blocks use 64th roots of unity, which are powers of
    // two. Index k = m + i as in psi_rev_ (k in [1, 32)).
    internals::Pow2Twiddle m_fwd_pow2_twiddle[32];
    internals::Pow2Twiddle m_backward_pow2_twiddle[32];
    internals::Pow2Twiddle m_final_twiddle_scaled;   // psi^-N/2 * N^-1, used on the last inverse level
    internals::Pow2Twiddle m_inverse_dim;  // N^-1 = 2^-log_n_
};

}  // namespace goldilocks

#endif  // LARGE_FUNCTIONS_GOLDILOCKS_WORKER_H
