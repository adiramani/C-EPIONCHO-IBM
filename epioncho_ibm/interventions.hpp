#ifndef INTERVENTIONS_HPP
#define INTERVENTIONS_HPP

#include "params.hpp"
#include <algorithm>
#include <iostream>
#include <string>

class Intervention {
protected:
    InterventionParams params;
    std::size_t index = 0;
    int curr_timesteps_per_year = -1;

public:
    std::vector<int> timestep_scaled_application_times;
    Intervention(const InterventionParams& params_) : params(params_) {}
    virtual ~Intervention() = default;


    virtual bool isMDA();
    virtual bool isVectorControl();
    virtual bool check_scale(int timesteps_per_year);
    virtual void scale_application_times(int timesteps_per_year);
    virtual bool during_intervention_period(int current_timestep);
    virtual bool should_apply(int current_timestep);
};

class Treatment : public Intervention {
protected:
    TreatmentParams treatment_params;

public:
    Treatment(const TreatmentParams& params);

    TreatmentParams apply();
};

class VectorControl : public Intervention {
protected:
    VectorControlParams vc_params;
public:
    VectorControl(const VectorControlParams& params);

    bool should_apply(int current_timestep);
    double apply(double initial_annual_biting_rate, int current_timestep);
};

#endif
