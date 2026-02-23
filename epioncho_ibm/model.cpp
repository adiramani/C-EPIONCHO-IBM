#include "model.hpp"
#include <iostream>
#include <time.h>

Model::Model(const InputParams& input_params)
    : state(input_params) // state(std::make_unique<State>(input_params))
{
    std::cout << "Initialized Model\n";
    std::cout << "Timestep in years: " << (state.params.base.delta_time_days / state.params.base.year_length_days) << "\n";
    // TODO: Initialize model state variables here
    // state->people = ...
}

double get_elapsed_time(clock_t start_time) {
    return (double)(clock() - start_time)/CLOCKS_PER_SEC;
}

void Model::advance_timestep() {
    clock_t overall_start_time = clock();
    clock_t start_time = overall_start_time;
    double current_year = state.current_timestep / state.params.base.year_length_days;

    
    // TODO: determine if treatment should be conducted this timestep
    
    double beta = state.params.blackfly.human_blood_index / state.params.blackfly.gonotrophic_cycle_length;
    auto exposure_vals = state.people.get_exposure(
        state.params.exposure,
        0.5
    );
    
    // Add new worms from L3 delay array (determine sex):
    std::vector<double> new_worms = state.people.get_new_worms();
    std::vector<double> new_male_worms(new_worms.size(), 0.0);
    std::vector<double> new_female_worms(new_worms.size(), 0.0);

    get_new_worms_time += get_elapsed_time(start_time);
    
    start_time = clock();
    for (int i = 0; i < state.people.population_size; ++i) {
        std::binomial_distribution<int> worm_sex_gen(new_worms[i], state.params.worms.sex_ratio);
        new_male_worms[i] = worm_sex_gen(state.generator);
        new_female_worms[i] = new_worms[i] - new_male_worms[i];
    }
    worm_sex_time += get_elapsed_time(start_time);

    start_time = clock();
    // Calculate rate of infection (new L3 to be added to delay array)
    state.people.new_established_l3(
        state.generator, state.params.blackfly.delta_h_zero, state.params.blackfly.delta_h_inf,
        state.params.blackfly.c_h, state.params.blackfly.bite_rate_per_person_per_year, exposure_vals,
        state.timestep_years
    );
    new_l3_time += get_elapsed_time(start_time);

    start_time = clock();
    // Age worms
    // Age/Calculate new MF
    // Reset individuals who die
    std::vector<double> mf_loads_curr;
    mf_loads_curr.reserve(state.people.population_size);
    for (auto& mf : state.people.microfilariae) {
        mf_loads_curr.push_back(mf.get_mf_load());
    }
    // TODO: Bottleneck
    state.people.age(state.generator, state.timestep_years, new_male_worms, new_female_worms); 
    age_time += get_elapsed_time(start_time);   
    auto& microfilariae = state.people.microfilariae;
    auto& blackflies = state.people.blackflies;

    double mf_delay = 0.0;
    double expos_delay = 0.0;

    // Update L1-L3 in blackfly (per indiv)
    start_time = clock();
    for (int i = 0; i < state.people.population_size; ++i) {
        mf_delay += state.people.blackflies[i].delays[state.people.blackflies[i].delay_index].microfilariae;
        expos_delay += state.people.blackflies[i].delays[state.people.blackflies[i].delay_index].exposure;
        blackflies[i].calc_L1(
            state.timestep_years,
            beta,
            exposure_vals[i],
            mf_loads_curr[i]
        );
        double curr_l2 = blackflies[i].l2;
        blackflies[i].calc_L2(state.timestep_years);
        blackflies[i].calc_L3(
            curr_l2,
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
    update_blackfly_time += get_elapsed_time(start_time);

    // Calculate morbidities
    // Update Ov16 serostatus

    // Process deaths
    start_time = clock();
    state.people.process_deaths(state.generator);
    process_death_time += get_elapsed_time(start_time);

    start_time = clock();
    double mf_positives = 0.0;
    double mf_load = 0.0;
    double worm_load = 0.0;
    double infertile_female_worm_load = 0.0;
    double mean_expos = 0.0;
    double age_tot = std::accumulate(state.people.ages.begin(), state.people.ages.end(), 0.0);
    for (int i = 0; i < state.people.population_size; ++i) {
        mf_load += microfilariae[i].get_mf_load();
        mf_positives += microfilariae[i].get_mf_load() > 0;
        worm_load += state.people.worms["fertile_female"][i].get_worm_load();
        infertile_female_worm_load += state.people.worms["infertile_female"][i].get_worm_load();
        mean_expos += exposure_vals[i];
    }

    calc_outputs_time += get_elapsed_time(start_time);

    start_time = clock();
    overall_time += get_elapsed_time(overall_start_time); 
    // TODO: Store output instead of printing
    if (current_year == std::floor(current_year)) {
        printf("Current Year: %f | ", current_year);
        printf("Average Age: %f | ", age_tot / state.people.population_size);
        printf("MF Prev: %f\n", mf_positives / state.people.population_size);
        // printf("MF Load: %f | ", mf_load / state.people.population_size);
        // printf("Fertile Female Worm Load: %f | ", worm_load / state.people.population_size);
        // printf("IF Load: %f | ", infertile_female_worm_load / state.people.population_size);
        // printf("MF Delay: %f | ", mf_delay / state.people.population_size);
        // printf("Expos Delay: %f | ", expos_delay / state.people.population_size);
        // printf("Mean Expos Delay: %f | ", mean_expos / state.people.population_size);
        // printf("Mean Blackfly L1: %f | ", state.people.mean_l1_per_blackfly());
        // printf("Mean Blackfly L2: %f | ", state.people.mean_l2_per_blackfly());
        // printf("Mean Blackfly L3: %f\n", state.people.mean_l3_per_blackfly());
        printf(
            "Avg Clocktimes: %f | %f | %f | %f | %f | %f | %f | %f | %f\n",
            overall_time / state.current_timestep, get_new_worms_time / state.current_timestep,
            worm_sex_time / state.current_timestep, new_l3_time / state.current_timestep,
            age_time / state.current_timestep, update_blackfly_time / state.current_timestep,
            process_death_time / state.current_timestep, calc_outputs_time / state.current_timestep,
            print_time / state.current_timestep
        );
    }
    print_time += get_elapsed_time(start_time);

    // iterate to next timestep
    state.current_timestep += state.params.base.delta_time_days;
}
