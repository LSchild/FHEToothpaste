//
// Created by leonard on 7/31/26.
//

#ifndef TOOTHPASTE_KS_FUNCTION_EVALUATOR_H
#define TOOTHPASTE_KS_FUNCTION_EVALUATOR_H

#include "interfaces/functional_bootstrap.h"
#include "interfaces/operator_context.h"
#include "interfaces/blindrotation_operator.h"
#include "endo_glwe_conversion.h"

struct KSFunctionEvaluator;

struct KSFunctionEvaluationContext : OperatorContext<KSFunctionEvaluator>, public std::enable_shared_from_this<KSFunctionEvaluationContext> {

    KSFunctionEvaluationContext(std::shared_ptr<OperatorContext<BlindRotator>> blind_rotation_context, std::shared_ptr<LWEConversionContext> converter, bool estimate_padding = false);

    KSFunctionEvaluationContext(std::shared_ptr<OperatorContext<BlindRotator>> blind_rotation_context, std::shared_ptr<LWEConversionContext> converter, uint64_t padding_width);

    [[nodiscard]] std::shared_ptr<OperatorContext<BlindRotator>> GetBlindRotationContext() const;

    [[nodiscard]] std::shared_ptr<LWEConversionContext> GetLWEConversionContext() const;

    [[nodiscard]] long double ComputeOutputVariance(long double input_variance = 0.0) const override;

    [[nodiscard]] std::unique_ptr<KSFunctionEvaluator> ConstructOperator(const std::vector<GenericKey>& keys) const override;

    /** Describes the input LWE sample
     *
     * @return A Container describing the expected input
     */
    [[nodiscard]] Container GetInputContainer() const override;

    /** Describes the output LWE sample. Note that several input containers may be present
     * - LWEContainer (ignored)
     * - RLWEContainer, if the bootstrapping accumulator is provided
     * - FlagContainer, to specify which output type is expected
     *
     * @param input Input container description.
     * @return LWEContainer
     */
    [[nodiscard]] Container GetOutputContainer(Container input) const override;

    /** Returns an ID describing the Operator
     *
     * @return The ID
     */
    [[nodiscard]] OperatorID GetOperatorID() const override;

    [[nodiscard]] std::shared_ptr<MathWorker> GetWorker() const;

private:

    bool m_estimate_padding;
    uint64_t m_padding_width;

    std::shared_ptr<MathWorker> m_worker;

    std::shared_ptr<OperatorContext<BlindRotator>> m_blind_rotation_context;
    std::shared_ptr<LWEConversionContext> m_conversion_context;

};

struct KSFunctionEvaluator : FunctionEvaluator<KSFunctionEvaluationContext> {

    friend struct KSFunctionEvaluationContext;

    /** Function evaluation with a provided negacyclic function mapping Z_{2N} to Z_Q
     * where N is the ring dimension and Q is the ring Q
     *
     * @param result Buffer for output LWE
     * @param function Function to be evaluated
     * @param output Stage at which to output, either right after blind-rotation or after key-switching.
     */
    void EvalFunc(uint64_t* RESTRICTED result, const uint64_t* RESTRICTED input_lwe, std::function<uint64_t(uint64_t)>& function) override;

    /** Function evaluation with a provided negacyclic function mapping Z_{2N} to Z_Q
     * where N is the ring dimension and Q is the ring Q
     *
     * @param result Buffer for output LWE
     * @param function Function to be evaluated
     * @param output Stage at which to output, either right after blind-rotation or after key-switching.
     */
    void EvalFunc(std::vector<uint64_t>& result, const std::vector<uint64_t>& input, std::function<uint64_t(uint64_t)>& function) override;

    /** Function evaluation with a (possibly encrypted) negacyclic LUT in lut_rlwe
     * where N is the ring dimension and Q is the ring Q
     *
     * @param result Buffer for output LWE
     * @param function LUT to perform the evaluation on
     * @param output Stage at which to output, either right after blind-rotation or after key-switching.
     */
    void EvalFunc(uint64_t* RESTRICTED result, const uint64_t* RESTRICTED input_lwe, uint64_t* RESTRICTED lut_rlwe) override;

    /** Function evaluation with a (possibly encrypted) negacyclic LUT in lut_rlwe
     * where N is the ring dimension and Q is the ring Q
     *
     * @param result Buffer for output LWE
     * @param function LUT to perform the evaluation on
     * @param output Stage at which to output, either right after blind-rotation or after key-switching.
     */
    void EvalFunc(std::vector<uint64_t>& result, const std::vector<uint64_t>& input,  std::vector<uint64_t>& lut_rlwe) override;

    [[nodiscard]] const std::shared_ptr<const KSFunctionEvaluationContext> GetContext() const override;

    uint64_t EstimatePadding(const uint64_t* poly) const;

    uint64_t EstimatePadding(std::vector<uint64_t>& poly) const;

    [[nodiscard]] uint64_t GetPaddingWidth() const;

    [[nodiscard]] uint64_t GetEstimatePadding() const;

    ~KSFunctionEvaluator() override = default;

protected:

    KSFunctionEvaluator(std::shared_ptr<const KSFunctionEvaluationContext> context, std::unique_ptr<BlindRotator> rotator, std::unique_ptr<LWEtoLWEConverter> converter, bool estimate_padding, uint64_t padding_width);

    std::shared_ptr<const KSFunctionEvaluationContext> m_context;

    bool m_estimate_padding;
    uint64_t m_padding_width;

    std::unique_ptr<BlindRotator> m_rotator;
    std::unique_ptr<LWEtoLWEConverter> m_converter;

    AlignedBuffer m_sign_poly_acc;
    AlignedBuffer m_padding_poly;


};

#endif //TOOTHPASTE_KS_FUNCTION_EVALUATOR_H
