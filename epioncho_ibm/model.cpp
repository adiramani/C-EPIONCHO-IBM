#include "model.hpp"
#include <iostream>
#include <numeric>
#include <ctime>
#include <fstream>

Model::Model(InputParams input_params, bool enable_timing)
    : state(std::move(input_params.params)),
    enable_timing(enable_timing)
{
    int timesteps_per_year = state.params.base.year_length_days / state.params.base.delta_time_days;
    for (auto& treat : input_params.treatments) {
        std::unique_ptr<Treatment> temp_treat = std::make_unique<Treatment>(treat);
        temp_treat->scale_application_times(timesteps_per_year);
        treatments.push_back(std::move(temp_treat));
    }
    for (auto& vc : input_params.vector_control) {
        std::unique_ptr<VectorControl> temp_vc = std::make_unique<VectorControl>(vc);
        temp_vc->scale_application_times(timesteps_per_year);
        vector_control.push_back(std::move(temp_vc));
    }
    // std::cout << "Initialized Model\n";
    // std::cout << "Timestep in years: "
    //           << (state.params.base.delta_time_days / state.params.base.year_length_days)
    //           << "\n";
    // std::cout << "Seed: " << state.params.base.seed << " | ABR: " << state.params.blackfly.bite_rate_per_person_per_year << " | kE: " << state.params.base.k_E;
    // std::cout << " | c_h: " << state.params.blackfly.c_h << "\n";
}

static double get_elapsed_time(clock_t start_time) {
    return (double)(clock() - start_time) / CLOCKS_PER_SEC;
}

std::vector<Intervention*> get_active_interventions(
    const std::vector<std::unique_ptr<Treatment>>& treatments,
    const std::vector<std::unique_ptr<VectorControl>>& vc,
    int current_timestep
) {
    std::vector<Intervention*> result;
    for (auto& i : treatments)
        if (i->during_intervention_period(current_timestep))
            result.push_back(i.get());
    for (auto& i : vc)
        if (i->during_intervention_period(current_timestep))
            result.push_back(i.get());
    return result;
}

void Model::advance_timestep(bool verbose) {
    clock_t start;
    clock_t overall_start;
    if (enable_timing) {
        overall_start = clock();
        start = overall_start;
    }
    state.people.generate_random_vals_for_diagnostic(state.generator);

    int current_timestep = state.current_timestep;

    double current_year = (double)current_timestep
                              / state.params.base.year_length_days;

    // Get all interventions that will be applied in the current year
    int num_vc_year = 0;
    int num_mda_year = 0;
    if (floor(state.current_year) != floor(current_year)) {
        state.current_interventions = get_active_interventions(
            treatments, vector_control, current_timestep
        );
    }
    state.current_year = current_year;

    if (verbose) {
        _print_yearly_stats(current_timestep, current_year, num_vc_year, num_mda_year, enable_timing);
    }


    double bite_rate = state.params.blackfly.bite_rate_per_person_per_year;

    // See which interventions are applied at the current timestep
    Treatment* applied_treatment = nullptr;
    int timesteps_per_year = state.params.base.year_length_days / state.params.base.delta_time_days;
    for (auto& intervention : state.current_interventions) {
        if (intervention->isVectorControl()) {
            num_vc_year += 1;
        } else if (intervention->isMDA()) {
            num_mda_year += 1;
        }
        if (!intervention->check_scale(timesteps_per_year)) {
            intervention->scale_application_times(timesteps_per_year);
        }
        if (intervention->should_apply(current_timestep)) {
            if (intervention->isVectorControl()) {
                bite_rate = static_cast<VectorControl*>(intervention)->apply(bite_rate, current_timestep);
            } else if (intervention->isMDA()) {
                applied_treatment = static_cast<Treatment*>(intervention);
            }
        }
    }

    // Apply MDA, determine who has been treated
    if (applied_treatment != nullptr) {
        TreatmentParams treatment_params = applied_treatment->apply();
        state.people.treatment_params = treatment_params;
        // Create or udpate compliance values
        // TODO: simplify compliance methods to just use the stored treatment_params object.
        if (state.people._current_rho == -1 && state.people._current_cov == -1) {
            state.people.calculate_population_compliance(
                state.generator, treatment_params.rho, treatment_params.total_population_coverage,
                treatment_params.proportion_never_treated
            );
        } else if (
                state.people._current_cov != treatment_params.total_population_coverage ||
                state.people._current_rho != treatment_params.rho
        ) {
            state.people.update_compliance(
                state.generator, treatment_params.rho, treatment_params.total_population_coverage
            );
        }

        state.people.apply_treatment_round(state.generator, treatment_params.min_age_of_treatment, current_timestep);
    }

    const auto exposure_vals = state.people.get_exposure(state.params.exposure, state.params.human.gender_ratio);

    std::vector<double> new_worms = state.people.get_new_worms();

    if (enable_timing) {
        get_new_worms_time += get_elapsed_time(start);
        start = clock();
    }

    std::vector<double> new_male_worms(state.people.population_size, 0.0);
    std::vector<double> new_female_worms(state.people.population_size, 0.0);
    for (int i = 0; i < state.people.population_size; ++i) {
        if ((int)new_worms[i] <= 0)
            continue;
        std::binomial_distribution<int> sex_dist(
            (int)new_worms[i], state.params.worms.sex_ratio
        );
        new_male_worms[i] = (double)sex_dist(state.generator);
        new_female_worms[i] = new_worms[i] - new_male_worms[i];
        if (new_male_worms[i] + new_female_worms[i] != new_worms[i])
            printf("We have an issue");
    }
    if (enable_timing) {
        worm_sex_time += get_elapsed_time(start);
        start = clock();
    }
    const double beta = state.params.blackfly.human_blood_index
                      / state.params.blackfly.gonotrophic_cycle_length;
    state.people.new_established_l3(
        state.generator,
        state.params.blackfly.delta_h_zero,
        state.params.blackfly.delta_h_inf,
        state.params.blackfly.c_h,
        bite_rate,
        exposure_vals,
        state.timestep_years
    );

    if (enable_timing) {
        new_l3_time += get_elapsed_time(start);
        start = clock();
    }

    state.people.microfilariae.get_all_raw_loads(state.people.temp_mf_loads);

    state.people.age(
        state.generator, current_timestep, state.timestep_years,
        new_male_worms, new_female_worms
    );

    if (enable_timing) {
        age_time += get_elapsed_time(start);
        start = clock();
    }

    state.people.blackflies.update_all(
        state.timestep_years,
        beta,
        exposure_vals,
        state.people.temp_mf_loads,
        state.params.blackfly.a_H,
        state.params.blackfly.gonotrophic_cycle_length,
        state.params.blackfly.mu_L3
    );

    if (enable_timing) {
        update_blackfly_time += get_elapsed_time(start);
        start = clock();
    }

    state.people.update_all_status(
        state.generator, state.timestep_years, state.params.base.year_length_days, 
        state.params.human.skin_snip_weight, state.params.human.skin_snip_number
    );

    if (enable_timing) {
        update_people_status_time += get_elapsed_time(start);
        start = clock();
    }

    state.people.process_deaths(state.generator);

    if (enable_timing) {
        process_death_time += get_elapsed_time(start);
        overall_time += get_elapsed_time(overall_start);
    }
    

    // if (verbose) {
    //     _print_yearly_stats(current_timestep, current_year, num_vc_year, num_mda_year, enable_timing);
    // }

    state.current_timestep += (int)state.params.base.delta_time_days;
}

void Model::_print_yearly_stats(
    int current_timestep, double current_year, int num_vc_year,
    int num_mda_year, bool enable_timing
) {
    clock_t start;
    if (enable_timing) {
        start = clock();
    }
    double mf_positives = 0.0;
    double mf_load = 0.0, male_worm_load = 0.0, fertile_worm_load = 0.0, infertile_worm_load = 0.0;
    const double age_tot = std::accumulate(
        state.people.ages.begin(), state.people.ages.end(), 0.0);

    for (int i = 0; i < state.people.population_size; ++i) {
        const double mf = state.people.temp_mf_loads[i];
        mf_load += mf;
        mf_positives += (mf > 0.0) ? 1.0 : 0.0;
        male_worm_load += state.people.male_worms.get_raw_load(i);
        fertile_worm_load += state.people.fertile_female_worms.get_raw_load(i);
        infertile_worm_load += state.people.infertile_female_worms.get_raw_load(i);
    }
    if (enable_timing) {
        calc_outputs_time += get_elapsed_time(start);
        start = clock();
    }

    if (current_year == std::floor(current_year) || true) {
        printf("Current Year: %f | ", current_year);
        printf("Avg MF, MW, WFF, WIF Load: %f, %f, %f, %f\n", mf_load / state.people.population_size, male_worm_load / state.people.population_size, fertile_worm_load / state.people.population_size, infertile_worm_load / state.people.population_size);
        if (enable_timing) {
            printf(
                "Avg Clocktimes: %.6f | %.6f | %.6f | %.6f | %.6f | %.6f | %.6f | %.6f | %.6f\n",
                overall_time / current_timestep,
                get_new_worms_time / current_timestep,
                worm_sex_time / current_timestep,
                new_l3_time / current_timestep,
                age_time / current_timestep,
                update_blackfly_time  / current_timestep,
                process_death_time / current_timestep,
                calc_outputs_time / current_timestep,
                print_time / current_timestep
            );
        }
    }
    if (enable_timing) {
        print_time += get_elapsed_time(start);
    }
}
