#include "sequelae.hpp"

Sequelae::Sequelae(
    std::vector<SequelaeType> sequelae_types, int num_individuals,
    std::vector<double> probabilities, std::vector<SequelaeProbTimeUnit> timescales,
    std::vector<int> countdown_timesteps, std::vector<double> average_ages
) : sequelae_types(sequelae_types),
    num_individuals(num_individuals),
    probabilities(probabilities),
    timescales(timescales),
    countdown_timesteps(countdown_timesteps),
    average_ages(average_ages)
{
    num_sequelae = sequelae_types.size();
    positives.resize(num_sequelae * num_individuals, false);
    has_been_tested.resize(num_sequelae * num_individuals, false);
    countdowns.resize(num_sequelae * num_individuals, 0);
    for (int i = 0; i < num_individuals; ++i) {
        for (int s = 0 ; s < sequelae_types.size(); ++s)
            countdowns[(i * num_sequelae) + s] = countdown_timesteps[s];
    }
}

double Sequelae::calculate_timestep_probability(int sequelae_index, double timestep_years, int days_in_year) {
    SequelaeProbTimeUnit timescale = timescales[sequelae_index];
    double average_age = average_ages[sequelae_index];
    double probability = probabilities[sequelae_index];

    double exponent = timestep_years;
    if (timescale == SequelaeProbTimeUnit::Day)
        exponent *= days_in_year;
    if (average_age > 0)
        exponent *= 1/average_age;
    return 1 - pow(1 - probability, exponent);
}

void Sequelae::decrease_countdown_individual(int indiv_index) {
    for (int s = 0; s < num_sequelae; ++s) {
        int person_sequelae_index = (indiv_index * num_sequelae) + s;
        bool should_decrease = positives[person_sequelae_index] && countdown_timesteps[s] > 0;
        countdowns[person_sequelae_index] -= should_decrease ? 1 : 0;
        if (should_decrease && countdowns[person_sequelae_index] < 0)
            positives[person_sequelae_index] = false;
    }
}

bool Sequelae::should_test_individual(
    SequelaeType st, int person_sequelae_index, double age,
    double raw_mf_count, double observed_mf_count
) {
    bool neg_and_not_tested = !has_been_tested[person_sequelae_index] && !positives[person_sequelae_index];
    switch (st) {
        case SequelaeType::Blindness: return (raw_mf_count > 0) && neg_and_not_tested;
        case SequelaeType::SevereItch: return (observed_mf_count > 0 & age >= 2) && neg_and_not_tested;
        case SequelaeType::ReactiveSkinDisease: return (observed_mf_count > 0 & age >= 2) && neg_and_not_tested;
        case SequelaeType::Atrophy: return (observed_mf_count > 0) && neg_and_not_tested;
        case SequelaeType::HangingGroin: return (observed_mf_count > 0) && neg_and_not_tested;
        case SequelaeType::Depigmentation: return (observed_mf_count > 0) && neg_and_not_tested;
        default: return neg_and_not_tested;
    }

}

void Sequelae::update_sequelae(SequelaeType st, int person_sequelae_index, bool positive) {
    switch (st) {
        case SequelaeType::Blindness: has_been_tested[person_sequelae_index] = positive;
        case SequelaeType::SevereItch: positives[person_sequelae_index] = positive;
        case SequelaeType::ReactiveSkinDisease: positives[person_sequelae_index] = positive;
        case SequelaeType::Atrophy: positives[person_sequelae_index] = positive;
        case SequelaeType::HangingGroin: positives[person_sequelae_index] = positive;
        case SequelaeType::Depigmentation: positives[person_sequelae_index] = positive;
        default: positives[person_sequelae_index] = positives[person_sequelae_index];
    }
}

void Sequelae::update_all_sequelae_indiv(
    std::mt19937& generator, std::uniform_real_distribution<double>& uniform_dist,
    int indiv_index, double timestep_years, 
    int days_in_year, double age, double raw_mf_count,
    double observed_mf_count
) {
    for (int s = 0; s < num_sequelae; ++s) {
        int person_sequelae_index = (indiv_index * num_sequelae) + s;
        if (!should_test_individual(sequelae_types[s], person_sequelae_index, age, raw_mf_count, observed_mf_count))
            continue;
        double probability = calculate_timestep_probability(s, timestep_years, days_in_year);
        update_sequelae(sequelae_types[s], person_sequelae_index, uniform_dist(generator) < probability);
    }
}

void Sequelae::process_death(int indiv_index) {
    for (int s = 0; s < num_sequelae; ++s) {
        int person_sequelae_index = (indiv_index * num_sequelae) + s;
        positives[person_sequelae_index] = false;
        has_been_tested[person_sequelae_index] = false;
        countdowns[person_sequelae_index] = countdown_timesteps[s];
    }
}

int Sequelae::indiv_status(SequelaeType sequelae_type, int indiv_index) {
    for (int s = 0; s < num_sequelae; ++s)
        return positives[(indiv_index * num_sequelae) + s];
    return -1.0;
}
