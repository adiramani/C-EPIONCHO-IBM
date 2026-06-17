# include "interventions.hpp"
#include <cmath>

bool Intervention::isMDA() {
    return params.intervention_type == InterventionType::MDA;
}

bool Intervention::isVectorControl() {
    return params.intervention_type == InterventionType::VectorControl;
}

bool Intervention::check_scale(int timesteps_per_year) {
    return curr_timesteps_per_year == timesteps_per_year;
}

void Intervention::scale_application_times(int timesteps_per_year) {
    curr_timesteps_per_year = timesteps_per_year;
    if (params.pre_converted_timesteps == true) {
        timesteps_per_year = 1;
    }
    timestep_scaled_application_times.resize(params.application_times.size());

    std::transform(
        params.application_times.begin(),
        params.application_times.end(),
        timestep_scaled_application_times.begin(),
        [timesteps_per_year](auto t) {
            return std::floor(t * timesteps_per_year);
        }
    );
}

bool Intervention::during_intervention_period(int current_timestep) {
    return (
        (current_timestep >= timestep_scaled_application_times.front()) && 
        (current_timestep <= timestep_scaled_application_times.back())
    );
}

bool Intervention::should_apply(int current_timestep) {
    return timestep_scaled_application_times[index] == current_timestep;
}

Treatment::Treatment(const TreatmentParams& params_)
    : Intervention(params_), treatment_params(params_) 
{}

TreatmentParams Treatment::apply() {
    index++;
    return treatment_params;
}

VectorControl::VectorControl(const VectorControlParams& params_)
    : Intervention(params_), vc_params(params_) 
{}

bool VectorControl::should_apply(int current_timestep) {
    return during_intervention_period(current_timestep);
}

double VectorControl::apply(double initial_annual_biting_rate, int current_time_step) {
    if (
        index + 1 < timestep_scaled_application_times.size() &&
        current_time_step >= timestep_scaled_application_times[index + 1]
    )
        index++;
    if (index >= vc_params.efficacies.size())
        throw std::out_of_range("Fatal: Model attemting to use VectorControl object longer than specified.");
    double updated_biting_rate = initial_annual_biting_rate * (1 - vc_params.efficacies[index]);
    return (
        updated_biting_rate
    );
}