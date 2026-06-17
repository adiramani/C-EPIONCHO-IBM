#include "sequelae.hpp"

Sequelae::Sequelae(
    SequelaeType sequelae_type, SequelaeModelType sequelae_model_type, int num_individuals,
    double base_probability, SequelaeProbTimeUnit prob_timescale,
    int min_age_test, int min_infection, bool retest_tested_indivs,
    int countdown_timesteps, bool status_end_countdown, bool use_raw_infection_for_test
) : sequelae_type(sequelae_type),
    smt(sequelae_model_type),
    num_individuals(num_individuals),
    base_probability(base_probability),
    prob_timescale(prob_timescale),
    min_age_test(min_age_test),
    min_infection(min_infection),
    retest_tested_indivs(retest_tested_indivs),
    countdown_timesteps(countdown_timesteps),
    status_end_countdown(status_end_countdown),
    use_raw_infection_for_test(use_raw_infection_for_test)
{
    positives.resize(num_individuals, false);
    has_been_tested.resize(num_individuals, 0);
    countdowns.resize(num_individuals, 0);
    for (int i = 0; i < num_individuals; ++i) {
        countdowns[i] = countdown_timesteps;
    }
}

SequelaeType Sequelae::get_type() {
    return sequelae_type;
}
SequelaeModelType Sequelae::get_model_type() {
    return smt;
}

bool Sequelae::use_raw_infection() {
    return use_raw_infection_for_test;
}

double Sequelae::get_probability(double timestep_years, int days_in_year, int infection_level) {
    double exponent = 1.0 / 365.0;//timestep_years;
    if (prob_timescale == SequelaeProbTimeUnit::Year)
        exponent *= 365.0;//days_in_year;
    return 1 - pow(1 - base_probability, exponent);
}

void Sequelae::decrease_countdown_individual(int indiv_index) {
    bool should_decrease = (
        (
            (positives[indiv_index] && !status_end_countdown) || 
            (!positives[indiv_index] && has_been_tested[indiv_index] == 2 && status_end_countdown)
        ) && 
        countdown_timesteps > 0
    );
    countdowns[indiv_index] -= should_decrease ? 1 : 0;
    if (should_decrease && countdowns[indiv_index] <= 0) {
        positives[indiv_index] = status_end_countdown;
    }
}

bool Sequelae::should_test_individual(
    int indiv_index, double age,
    int infection_level,
    bool mating_worm_pair, bool male_female_worm_pair
) {
    bool should_retest = retest_tested_indivs || has_been_tested[indiv_index] == 0;
    bool not_tested_positive = has_been_tested[indiv_index] != 2;
    bool negative = !positives[indiv_index];
    bool above_min_age = age >= min_age_test;
    bool infection_test = infection_level > min_infection;
    return should_retest && not_tested_positive && negative && above_min_age && infection_test;
}

void Sequelae::update_all_individuals(
    std::mt19937& generator, std::uniform_real_distribution<double>& uniform_dist,
    double timestep_years, int days_in_year, int num_individuals,
    const std::vector<double>& ages, const std::vector<int>& infection_levels_raw,
    const std::vector<int>& infection_levels_ss, bool mating_worm_pair, bool male_female_worm_pair
) {
    for (int indiv_index = 0; indiv_index < num_individuals; ++indiv_index) {
        if (
            !should_test_individual(
                indiv_index, ages[indiv_index], 
                use_raw_infection_for_test ? infection_levels_raw[indiv_index] : infection_levels_ss[indiv_index],
                mating_worm_pair, male_female_worm_pair
            )
        ) {
            decrease_countdown_individual(indiv_index);
            continue;
        }
        double probability = get_probability(timestep_years, days_in_year, infection_levels_ss[indiv_index]);
        bool positive = uniform_dist(generator) < probability;
        if (positive)
            countdowns[indiv_index] = countdown_timesteps;

        if (!retest_tested_indivs)
            has_been_tested[indiv_index] = 1;

        if ((countdown_timesteps > 0 && status_end_countdown == true) && positive)
            has_been_tested[indiv_index] = 2;
        else
            positives[indiv_index] = positive;
    }
}

void Sequelae::process_death(int indiv_index) {
    positives[indiv_index] = false;
    has_been_tested[indiv_index] = 0;
    countdowns[indiv_index] = countdown_timesteps;
}

int Sequelae::indiv_status(int indiv_index) {
    return positives[indiv_index];
}

TimestepProbSequelae::TimestepProbSequelae(
    SequelaeType sequelae_type, SequelaeModelType sequelae_model_type, int num_individuals,
    double base_probability, SequelaeProbTimeUnit prob_timescale,
    int min_age_test, int min_infection, bool retest_tested_indivs,
    int countdown_timesteps, bool status_end_countdown, bool use_raw_infection_for_test,
    double average_age
) : Sequelae(
        sequelae_type, sequelae_model_type, num_individuals, base_probability,
        prob_timescale,  min_age_test, min_infection, 
        retest_tested_indivs, countdown_timesteps,
        status_end_countdown, use_raw_infection_for_test
    ),
    average_age(average_age)
{}

double TimestepProbSequelae::get_probability(double timestep_years, int days_in_year, int infection_level) {
    double exponent = 1.0 / 365.0;//timestep_years;
    if (prob_timescale == SequelaeProbTimeUnit::Year)
        exponent *= 365.0;//days_in_year;
    exponent *= 1 / average_age;
    return 1 - pow(1 - base_probability, exponent);
}

ExponentialProbSequelae::ExponentialProbSequelae(
    SequelaeType sequelae_type, SequelaeModelType sequelae_model_type, int num_individuals,
    double base_probability, SequelaeProbTimeUnit prob_timescale,
    int min_age_test, int min_infection, bool retest_tested_indivs,
    int countdown_timesteps, bool status_end_countdown, bool use_raw_infection_for_test,
    double prob_intercept, double prob_slope
) : Sequelae(
        sequelae_type, sequelae_model_type, num_individuals, base_probability,
        prob_timescale,  min_age_test, min_infection, 
        retest_tested_indivs, countdown_timesteps,
        status_end_countdown, use_raw_infection_for_test
    ),
    prob_intercept(prob_intercept),
    prob_slope(prob_slope)
{}

double ExponentialProbSequelae::get_probability(double timestep_years, int days_in_year, int mf_count) {
    double exponent = 1.0 / 365.0;//timestep_years;
    if (prob_timescale == SequelaeProbTimeUnit::Year)
        exponent *= 365.0;//days_in_year;
    return (
        1 - 
        pow(
            1 - (prob_intercept * std::exp(prob_slope * mf_count)),
            exponent
        )
    );
}

PowerLawProbSequelae::PowerLawProbSequelae(
    SequelaeType sequelae_type, SequelaeModelType sequelae_model_type, int num_individuals,
    double base_probability, SequelaeProbTimeUnit prob_timescale,
    int min_age_test, int min_infection, bool retest_tested_indivs,
    int countdown_timesteps, bool status_end_countdown, bool use_raw_infection_for_test,
    double prob_intercept, double prob_slope
) : Sequelae(
        sequelae_type, sequelae_model_type, num_individuals, base_probability,
        prob_timescale,  min_age_test, min_infection, 
        retest_tested_indivs, countdown_timesteps,
        status_end_countdown, use_raw_infection_for_test
    ),
    prob_intercept(prob_intercept),
    prob_slope(prob_slope)
{}

double PowerLawProbSequelae::get_probability(double timestep_years, int days_in_year, int mf_count) {
    return std::exp(prob_intercept + prob_slope * std::log((double)mf_count + 1.0));
}

OAESequelae::OAESequelae(
    SequelaeType sequelae_type, SequelaeModelType sequelae_model_type, int num_individuals,
    double base_probability, SequelaeProbTimeUnit prob_timescale,
    int min_age_test, int min_infection, bool retest_tested_indivs,
    int countdown_timesteps, bool status_end_countdown, bool use_raw_infection_for_test, 
    double prob_intercept, double prob_slope,
    std::mt19937& gen, int max_age_to_test,
    std::vector<double> initial_ages
) : PowerLawProbSequelae(
        sequelae_type, sequelae_model_type, num_individuals,
        base_probability, prob_timescale,
        min_age_test, min_infection, retest_tested_indivs,
        countdown_timesteps, status_end_countdown, use_raw_infection_for_test, 
        prob_intercept, prob_slope
    ),
    gen(gen),
    ages_to_test(num_individuals)
{
    random_age_to_test = std::uniform_real_distribution<double>((double)min_age_test, (double)max_age_to_test);
    // Makes sure to test everyone at the proper time when initiating the model
    for (int i = 0; i < num_individuals; ++i) {
        bool outside_of_age_range = initial_ages[i] > max_age_to_test;
        ages_to_test[i] = outside_of_age_range ? round(initial_ages[i] + 1) : round(random_age_to_test(gen));
    }
}

bool OAESequelae::should_test_individual(
    int indiv_index, double age, int infection_level, bool mating_worm_pair, bool male_female_worm_pair
) {
    bool has_not_been_tested = has_been_tested[indiv_index] == 0;
    bool is_age_to_test = round(age) == ages_to_test[indiv_index];
    return has_not_been_tested && is_age_to_test && male_female_worm_pair;
}

void OAESequelae::process_death(int indiv_index) {
    ages_to_test[indiv_index] = round(random_age_to_test(gen));
    Sequelae::process_death(indiv_index);
}
