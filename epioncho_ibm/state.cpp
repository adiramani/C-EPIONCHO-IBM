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
        params.base.k_E,
        params.human,
        params.worms,
        params.mf,
        params.blackfly,
        params.sequelae_probs,
        params.base.sequela_active
    );
}

std::vector<int> State::get_sub_population(int age_start, int age_end) {
    return people.get_sub_population(age_start, age_end);
}

double State::mf_intensity(const std::vector<int>& filtered_individuals, int incubation_time_hours) {
    double total_intensity = 0.0;
    if (incubation_time_hours > 100) {
        total_intensity += 100;
    }
    for (auto& indiv : filtered_individuals) {
        total_intensity += people.microfilariae.get_skin_snip_load_person(generator, indiv, params.human.skin_snip_weight, params.human.skin_snip_number);
    }
    return total_intensity / filtered_individuals.size();
}

double State::mf_prevalence(const std::vector<int>& filtered_individuals, int incubation_time_hours) {
    double total_prevalence = 0.0;
    if (incubation_time_hours > 100) {
        return 0.0;
    }
    for (auto& indiv : filtered_individuals) {
        total_prevalence += people.microfilariae.get_skin_snip_load_person(generator, indiv, params.human.skin_snip_weight, params.human.skin_snip_number) > 0;
    }
    return total_prevalence / filtered_individuals.size();
}

void State::update_state_summary(const std::vector<int>& filtered_individuals, StateSummary& summary, double sens, double spec) {
    summary.total_individuals = filtered_individuals.size();
    for (auto& indiv : filtered_individuals) {
        summary.pop_male_worm_load += people.male_worms.get_raw_load(indiv);
        summary.pop_fertile_female_worm_load += people.fertile_female_worms.get_raw_load(indiv);
        summary.pop_infertile_female_worm_load += people.infertile_female_worms.get_raw_load(indiv);
        summary.pop_perm_sterile_female_worm_load += people.sterilized_female_worms.get_raw_load(indiv);

        summary.pop_raw_microfilarial_load += people.microfilariae.get_raw_load(indiv);

        double skin_snip_mf_count = people.microfilariae.get_skin_snip_load_person(generator, indiv, params.human.skin_snip_weight, params.human.skin_snip_number);
        summary.pop_skin_snip_microfilarial_load += skin_snip_mf_count;
        summary.skin_snip_positives += skin_snip_mf_count > 0;

        summary.raw_ov16_seropositives += people.ov16_serostatus[indiv];
        summary.adjusted_ov16_seropositives += people.sample_serostatus_individual(indiv, sens, spec);

        for (size_t s = 0; s < params.base.sequela_active.size(); ++s) {
            summary.sequelae_positives[s] += people.sequelae.indiv_status(summary.active_sequelae[s], indiv);
        }
    }
}