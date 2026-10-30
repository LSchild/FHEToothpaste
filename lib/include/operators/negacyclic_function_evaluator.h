//
// Created by leonard on 7/27/26.
//

#ifndef TOOTHPASTE_NEGACYCLIC_FUNCTION_EVALUATOR_H
#define TOOTHPASTE_NEGACYCLIC_FUNCTION_EVALUATOR_H

#include "negacyclic_function_evaluator.h"
#include "interfaces/functional_bootstrap.h"
#include "interfaces/operator_context.h"
#include "interfaces/blindrotation_operator.h"

#include "operators/endo_glwe_conversion.h"

struct NegacyclicFunctionEvaluator;

struct NegacyclicFunctionEvaluationContext : OperatorContext<NegacyclicFunctionEvaluator>, public std::enable_shared_from_this<NegacyclicFunctionEvaluationContext> {

    NegacyclicFunctionEvaluationContext(std::shared_ptr<OperatorContext<BlindRotator>> blind_rotation_context);

    [[nodiscard]] std::shared_ptr<OperatorContext<BlindRotator>> GetBlindRotationContext() const;

    [[nodiscard]] std::shared_ptr<LWEConversionContext> GetLWEConversionContext() const;

    [[nodiscard]] long double ComputeOutputVariance(long double input_variance = 0.0) const override;

    [[nodiscard]] std::unique_ptr<NegacyclicFunctionEvaluator> ConstructOperator(const std::vector<GenericKey>& keys) const override;

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

private:

    std::shared_ptr<OperatorContext<BlindRotator>> m_rotation_context;

};

struct NegacyclicFunctionEvaluator : FunctionEvaluator<NegacyclicFunctionEvaluationContext> {

    friend struct NegacyclicFunctionEvaluationContext;

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
     * @param input Input LWE
     * @param lut_rlwe LUT as vector
     */
    void EvalFunc(std::vector<uint64_t>& result, const std::vector<uint64_t>& input,  std::vector<uint64_t>& lut_rlwe) override;

    [[nodiscard]] const std::shared_ptr<const NegacyclicFunctionEvaluationContext> GetContext() const override;

    ~NegacyclicFunctionEvaluator() override = default;

private:


    NegacyclicFunctionEvaluator(std::shared_ptr<const NegacyclicFunctionEvaluationContext> ctx,
                                std::unique_ptr<BlindRotator> rotator);

    /* pointer to context that constructed the operator */
    std::shared_ptr<const NegacyclicFunctionEvaluationContext> m_context;
    /* internal blind-rotator */
    std::unique_ptr<BlindRotator> m_rotator;
    /* temporary buffer for blind-rotation */
    AlignedBuffer m_acc_buffer;

};


#endif //TOOTHPASTE_NEGACYCLIC_FUNCTION_EVALUATOR_H
