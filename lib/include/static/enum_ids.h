//
// Created by leonard on 5/27/26.
//

#ifndef LARGE_FUNCTIONS_ENUM_IDS_H
#define LARGE_FUNCTIONS_ENUM_IDS_H

enum PolyFormat {
    COEF,
    NTT
};

enum KeyDistribution {
    BINARY,
    TERNARY,
    GAUSSIAN
};

enum BlindRotationMethod {
    CGGI,
    BMMP
};

enum OperatorID {
    BR_CGGI,
    BR_BMMP,
    CONV_RLWE_RLWE,
    CONV_LWE_LWE,
    CONV_LWE_RGSW,
    DECOMP_SAP,
    DECOMP_LMP,
    EVAL_AUTO,
    EVAL_TRACE,
    PACK_LWE_TO_RLWE,
    FUNC_BOOT_KS,
    NEGACYCLIC_FUNC_BOOT
};

#endif //LARGE_FUNCTIONS_ENUM_IDS_H
