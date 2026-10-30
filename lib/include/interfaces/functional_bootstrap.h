//
// Created by leonard on 7/27/26.
//

#ifndef TOOTHPASTE_FUNCTIONAL_BOOTSTRAP_H
#define TOOTHPASTE_FUNCTIONAL_BOOTSTRAP_H

#include <cstdint>
#include <functional>
#include <future>
#include <string>

#include "backend/definitions.h"


template<typename ContextType>
struct FunctionEvaluator {

    virtual void EvalFunc(uint64_t* RESTRICTED result, const uint64_t* RESTRICTED input_lwe, std::function<uint64_t(uint64_t)>& function) = 0;

    virtual void EvalFunc(std::vector<uint64_t>& result, const std::vector<uint64_t>& input_lwe, std::function<uint64_t(uint64_t)>& function) = 0;

    virtual void EvalFunc(uint64_t* RESTRICTED result, const uint64_t* RESTRICTED input_lwe, uint64_t* RESTRICTED lut_rlwe) = 0;

    virtual void EvalFunc(std::vector<uint64_t>& result, const std::vector<uint64_t>& input_lwe, std::vector<uint64_t>& lut_rlwe) = 0;

    virtual const std::shared_ptr<const ContextType> GetContext() const = 0;

    virtual ~FunctionEvaluator() = default;
};

#endif //TOOTHPASTE_FUNCTIONAL_BOOTSTRAP_H
