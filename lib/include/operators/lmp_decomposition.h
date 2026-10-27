//
// Created by leonard on 10/27/26.
//

#ifndef TOOTHPASTE_LMP_DECOMPOSITION_H
#define TOOTHPASTE_LMP_DECOMPOSITION_H

#include "interfaces/decomposition_operator.h"
#include "interfaces/blindrotation_operator.h"
#include "operators/endo_glwe_conversion.h"

struct LMPDecomposer;

struct LMPDecompositionContext : OperatorContext<LMPDecomposer>, public std::enable_shared_from_this<LMPDecompositionContext> {

    LMPDecompositionContext(std::shared_ptr<OperatorContext<BlindRotator>> blind_rotation_context,
                            std::shared_ptr<LWEConversionContext> lwe_conversion_context,
                            uint64_t default_radix);


    /**
     * Computes output variance as a function of the input variance
     * In this case, the input is effectively ignored as the method is independent of the input variance
     * @param input_variance (ignored)
     * @return Output variance
     */
    long double ComputeOutputVariance(long double input_variance = 0) const override;

    /** Gets a pointer to the Rotation context
     *
     * @return BlindRotator
     */
    [[nodiscard]] std::shared_ptr<OperatorContext<BlindRotator>> GetBlindRotationContext() const;

    /** Gets a pointer to the Rotation context
     *
     * @return LWE to LWE conversion context
     */
    [[nodiscard]] std::shared_ptr<LWEConversionContext> GetLWEConversionContext() const;

    /** Returns default radix for decomposition
     *
     * @return value of default radix
     */
    [[nodiscard]] const uint64_t GetDefaultRadix() const;

    /**
     * Sets m_rotation_context to new rotation context
     * @param new_rot_context new rotation context
     */
    void SetBlindRotationContext(std::shared_ptr<OperatorContext<BlindRotator>> new_rot_context);


    /**
     * Sets m_conversion_context to new conversion context
     * @param new_packing_context new conversion context
     */
    void SetLWEConversionContext(std::shared_ptr<LWEConversionContext> new_lwe_context);

    /**
     * Sets new default radix
     * @param new_radix new radix value
     */
    void SetDefaultRadix(uint64_t new_radix);

    /**
     * Returns pointer to container describing the input format
     * @return pointer to expected input container
     */
    [[nodiscard]] Container GetInputContainer() const override;

    /** Returns pointer to container describing the output as a function of the input
     *
     * @param container input container (ignored here)
     * @return output container
     */
    [[nodiscard]] Container GetOutputContainer(Container container) const override;

    /**
     * Constructs a new SAP decomposition context
     * @param keys Key bundle including keys[0] = key compatible with the given blind-rotation context
     * @return unique_ptr to SAPDecomposer instance
     */
    std::unique_ptr<LMPDecomposer> ConstructOperator(const std::vector<GenericKey>& keys) const override;

    /**
     *
     * @return Operator id (DECOMP_LMP)
     */
    [[nodiscard]] OperatorID GetOperatorID() const override;

private:

    /* blind-rotation interface context */
    std::shared_ptr<OperatorContext<BlindRotator>> m_rotation_context;

    /* LWE to LWE conversion context */
    std::shared_ptr<LWEConversionContext> m_conversion_context;

    /* default radix */
    uint64_t m_default_radix;
};

struct LMPDecomposer : RadixDecomposer<LMPDecompositionContext> {

    friend struct LMPDecompositionContext;


    /** Performs the decomposition
     *
     * Note that the radix parameter usually needs to be specified for
     * both the associated context, and the method.
     * This is due to the fact that often, for any fixed basis B, the operator will also work
     * for bases B' either s.t. B' | B or B' < B
     *
     * @param output Pointer to output buffer
     * @param input Pointer to input buffer
     * @param radix Radix
     */
    void Decompose(uint64_t *output, const uint64_t *const input, uint64_t radix) override;

    /** Performs the decomposition
    *
    * Note that the radix parameter usually needs to be specified for
    * both the associated context, and the method.
    * This is due to the fact that often, for any fixed basis B, the operator will also work
    * for bases B' either s.t. B' | B or B' < B
    *
    * @param output Vector for output
    * @param input Vector for input
    * @param radix Radix
    */
    void Decompose(std::vector<uint64_t> &output, const std::vector<uint64_t> &input, uint64_t radix) override;

    /** Returns pointer to context that created the current operator.
    *
    * @return pointer to context
    */
    [[nodiscard]] const std::shared_ptr<const LMPDecompositionContext> GetContext() const override;

    ~LMPDecomposer() override = default;

private:

    LMPDecomposer(std::shared_ptr<const LMPDecompositionContext> ctx,
                  std::unique_ptr<BlindRotator> rotator,
                  std::unique_ptr<LWEtoLWEConverter> converter, uint64_t beta);

    void EvalF0(uint64_t* output, uint64_t* input);

    void EvalF1(uint64_t* output, uint64_t* input);

    /* context storing constants etc. */
    std::shared_ptr<const LMPDecompositionContext> m_context;
    /* pointer to interface instance for blind-rotation */
    std::unique_ptr<BlindRotator> m_rotator;
    /* pointer to lwe-to-lwe keyswitch for reset/restart */
    std::unique_ptr<LWEtoLWEConverter> m_lwe_converter;

    uint64_t m_beta;
};

#endif //TOOTHPASTE_LMP_DECOMPOSITION_H
