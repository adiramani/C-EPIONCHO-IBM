#ifndef ONCHO_PARAMS_HPP
#define ONCHO_PARAMS_HPP

#include "params.hpp"

struct OnchoSevereItchSequelaeParams : TimestepProbSequelaeParams {
    OnchoSevereItchSequelaeParams()
    : TimestepProbSequelaeParams(
        SequelaeType::SevereItch, SequelaeModelType::TimestepProb,
        0.1636701, SequelaeProbTimeUnit::Year, 2, 0, true, 3, false, false, 1
    )
    {}
};

struct OnchoRSDSequelaeParams : TimestepProbSequelaeParams {
    OnchoRSDSequelaeParams()
    : TimestepProbSequelaeParams(
        SequelaeType::ReactiveSkinDisease, SequelaeModelType::TimestepProb,
        0.04163095, SequelaeProbTimeUnit::Year, 2, 0, true, 3, false, false, 1
    )
    {}
};

struct OnchoAtrophySequelaeParams : TimestepProbSequelaeParams {
    OnchoAtrophySequelaeParams()
    : TimestepProbSequelaeParams(
        SequelaeType::Atrophy, SequelaeModelType::TimestepProb,
        0.036, SequelaeProbTimeUnit::Day, 0, 0, true, -1, false, false, 22.876
    )
    {}
};

struct OnchoHangingGroinSequelaeParams : TimestepProbSequelaeParams {
    OnchoHangingGroinSequelaeParams()
    : TimestepProbSequelaeParams(
        SequelaeType::HangingGroin, SequelaeModelType::TimestepProb,
        0.018, SequelaeProbTimeUnit::Day, 0, 0, true, -1, false, false, 22.876
    )
    {}
};

struct OnchoDepigmentationSequelaeParams : TimestepProbSequelaeParams {
    OnchoDepigmentationSequelaeParams()
    : TimestepProbSequelaeParams(
        SequelaeType::Depigmentation, SequelaeModelType::TimestepProb,
        0.059, SequelaeProbTimeUnit::Day, 0, 0, true, -1, false, false, 22.876
    )
    {}
};

struct OnchoBlindnessSequelaeParams : ExponentialProbSequelaeParams {
    OnchoBlindnessSequelaeParams()
    : ExponentialProbSequelaeParams(
        SequelaeType::Blindness, SequelaeModelType::ExponentialProb,
        0, SequelaeProbTimeUnit::Day, 0, 0, true, 730, true, true,
        0.003, 0.0099
    )
    {}
};

struct OnchoOAESequelaeParams : OAESequelaeParams {
    OnchoOAESequelaeParams()
    : OAESequelaeParams(
        SequelaeType::OAE, SequelaeModelType::OAE,
        0, SequelaeProbTimeUnit::NA, 3, 0, false, -1, false, true,
        -3.8934511926365, 0.42122568578622, 15
    )
    {}
};

static std::vector<std::unique_ptr<SequelaeParams>> get_all_oncho_sequelae_params() {
    std::vector<std::unique_ptr<SequelaeParams>> params;
    params.push_back(std::make_unique<OnchoSevereItchSequelaeParams>());
    params.push_back(std::make_unique<OnchoRSDSequelaeParams>());
    params.push_back(std::make_unique<OnchoAtrophySequelaeParams>());
    params.push_back(std::make_unique<OnchoHangingGroinSequelaeParams>());
    params.push_back(std::make_unique<OnchoDepigmentationSequelaeParams>());
    params.push_back(std::make_unique<OnchoBlindnessSequelaeParams>());
    params.push_back(std::make_unique<OnchoOAESequelaeParams>());
    return params;
}

struct DrugParamsIVM : DrugParams {
    double microfilaricidal_upsilon = 0.0096;
    double microfilaricidal_kappa = 1.25;
    double embryostatic_lambda_max = 32.4;
    double embryostatic_phi = 19.6;
    double permanent_infertility = 0.345;
};

struct DrugParamsMOX : DrugParams {
    double microfilaricidal_upsilon = 0.04;
    double microfilaricidal_kappa = 1.82;
    double embryostatic_lambda_max = 462;
    double embryostatic_phi = 4.83;
    double permanent_infertility = 0.345;
};

#endif