#include "model.hpp"
#include <iostream>

Model::Model(const InputParams& input_params)
    : state(input_params) // state(std::make_unique<State>(input_params))
{
    std::cout << "Initialized Model\n";
    std::cout << "Timestep in years: " << (state.params.base.delta_time_days / state.params.base.year_length_days) << "\n";
    // TODO: Initialize model state variables here
    // state->people = ...
}

void Model::advance_timestep() {
    double current_year = state.current_timestep / state.params.base.year_length_days;
    
    // TODO: determine if treatment should be conducted this timestep
    
    double beta = state.params.blackfly.human_blood_index / state.params.blackfly.gonotrophic_cycle_length;
    // double m = state.params.blackfly.bite_rate_per_person_per_year / beta;
    auto exposure_vals = state.people.get_exposure(
        state.params.exposure,
        0.5
    );
    
    // Add new worms from L3 delay array (determine sex):
    std::vector<double> new_worms = state.people.get_new_worms();
    std::vector<double> new_male_worms(new_worms.size(), 0.0);
    std::vector<double> new_female_worms(new_worms.size(), 0.0);

    for (int i = 0; i < state.people.population_size; ++i) {
        for (int j = 0; j < new_worms[i]; ++i) {
            new_male_worms[i] += state.people.gender_dist(state.generator);
        }
        new_female_worms[i] = new_worms[i] - new_male_worms[i];
    }

    // Calculate rate of infection (new L3 to be added to delay array)
    state.people.new_established_l3(
        state.generator, state.params.blackfly.delta_h_zero, state.params.blackfly.delta_h_inf,
        state.params.blackfly.c_h, state.params.blackfly.bite_rate_per_person_per_year, exposure_vals,
        state.timestep_years
    );
    
    // Age worms
    // Age/Calculate new MF
    // Reset individuals who die
    state.people.age(state.generator, state.timestep_years, new_male_worms, new_female_worms);    
    auto& microfilariae = state.people.microfilariae;
    auto& blackflies = state.people.blackflies;

    // Update L1-L3 in blackfly (per indiv)
    for (int i = 0; i < state.people.population_size; ++i) {
        blackflies[i].calc_L1(
            state.timestep_years,
            beta,
            exposure_vals[i],
            microfilariae[i].get_mf_load()
        );
        blackflies[i].calc_L2(state.timestep_years);
        blackflies[i].calc_L3(
            state.params.blackfly.a_H,
            state.params.blackfly.gonotrophic_cycle_length,
            state.params.blackfly.mu_L3
        );

        blackflies[i].update_delay_index(
            blackflies[i].l1,
            microfilariae[i].get_mf_load(),
            exposure_vals[i]
        );
    }

    // Calculate morbidities
    // Update Ov16 serostatus

    // Process deaths
    state.people.process_deaths(state.generator);

    double mf_positives = 0.0;
    double mf_load = 0.0;
    double worm_load = 0.0;
    for (int i = 0; i < state.people.population_size; ++i) {
        mf_load += microfilariae[i].get_mf_load();
        mf_positives += microfilariae[i].get_mf_load() > 0;
        worm_load += state.people.worms["fertile_female"][i].get_worm_load();
    }

    if (current_year == std::floor(current_year)) {
        printf("Current Year: %f | ", current_year);
        printf("MF Prev: %f | ", mf_positives / state.people.population_size);
        printf("MF Load: %f | ", mf_load / state.people.population_size);
        printf("Fertile Female Worm Load: %f | ", worm_load / state.people.population_size);
        printf("Mean Blackfly L3: %f\n", state.people.mean_l3_per_blackfly());
    }

    // iterate to next timestep
    state.current_timestep += state.params.base.delta_time_days;
}
