#include "state.hpp"

State::State(const Params& params)
    : params(params),
      timestep_years(
        params.base.delta_time_days
        / params.base.year_length_days
      ),
      people(params.base.n_people)
{
    generator = std::mt19937(params.base.seed);

    people.initialize_from_params(
        generator,
        timestep_years,
        params.human.mean_human_age,
        params.human.max_human_age,
        params.base.k_E,
        params.worms,
        params.mf,
        params.blackfly
    );
}

std::vector<int> State::get_sub_population(int age_start, int age_end) {
    return people.get_sub_population(age_start, age_end);
}

double State::mf_intensity(std::vector<int> filtered_individuals, int incubation_time_hours) {
    double total_intensity = 0.0;
    if (incubation_time_hours > 100) {
        total_intensity += 100;
    }
    for (auto& indiv : filtered_individuals) {
        total_intensity += people.microfilariae.get_skin_snip_load_person(generator, indiv, params.human.skin_snip_weight, params.human.skin_snip_number);
    }
    return total_intensity / filtered_individuals.size();
}

double State::mf_prevalence(std::vector<int> filtered_individuals, int incubation_time_hours) {
    double total_prevalence = 0.0;
    if (incubation_time_hours > 100) {
        return 0.0;
    }
    for (auto& indiv : filtered_individuals) {
        total_prevalence += people.microfilariae.get_skin_snip_load_person(generator, indiv, params.human.skin_snip_weight, params.human.skin_snip_number) > 0;
    }
    return total_prevalence / filtered_individuals.size();
}