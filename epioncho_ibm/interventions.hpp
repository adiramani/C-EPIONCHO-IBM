#ifndef INTERVENTIONS_HPP
#define INTERVENTIONS_HPP

#include "params.hpp"
#include <iostream>
#include <string>

// -------------------- Base Class --------------------
class Intervention {
protected:
    InterventionParams params;

public:
    Intervention(const InterventionParams& params_) : params(params_) {}
    virtual ~Intervention() = default;

    virtual bool should_apply(int current_timestep) {
        return std::find(
            params.application_timesteps.begin(),
            params.application_timesteps.end(),
            current_timestep
        ) != params.application_timesteps.end();
    }
};

class Treatment : public Intervention {
protected:
    TreatmentParams treatment_params;

public:
    Treatment();
    Treatment(const TreatmentParams& params);
};

// -------------------- Derived Class: VectorControl --------------------
class VectorControl : public Intervention {
protected:
    VectorControlParams vc_params;
public:
    VectorControl();
    VectorControl(const VectorControlParams& params);
};

#endif
