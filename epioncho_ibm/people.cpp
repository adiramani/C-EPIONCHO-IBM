#include "people.hpp"
#include <algorithm>
#include <numeric>
#include <cassert>
#include <cmath>
#include <random>


static double beta_distribution(std::mt19937& gen, double alpha, double beta) {
    std::gamma_distribution<double> dist_alpha(alpha, 1.0);
    std::gamma_distribution<double> dist_beta(beta, 1.0);
    double x = dist_alpha(gen);
    double y = dist_beta(gen);
    return x / (x + y);
}

People::People(int population_size_)
    : population_size(population_size_)
{}

void People::initialize_from_params(
    std::mt19937& gen,
    double timestep_years,
    double delta_time_days,
    double k_e,
    const HumanParams& human_params,
    const WormParams& worm_params,
    const MicrofilariaeParams& mf_params,
    const BlackflyParams& blackfly_params,
    const std::vector<std::unique_ptr<SequelaeParams>>& sequelae_params
)
{
    mean_age = human_params.mean_human_age;
    max_age = human_params.max_human_age;
    gender_ratio = human_params.gender_ratio;
    prop_serorevert_fast = human_params.prop_serorevert_fast;

    gender_dist = std::bernoulli_distribution(0.5);
    gamma_dist = std::gamma_distribution<double>(k_e, 1.0 / k_e);
    death_dist = std::bernoulli_distribution(1.0 / mean_age * timestep_years);
    serorevert_fast_dist = std::bernoulli_distribution(prop_serorevert_fast);
    uniform_dist = std::uniform_real_distribution(0.0, 1.0);

    ages.resize(population_size);
    sex.resize(population_size);
    exposure_heterogeneity.resize(population_size);
    ov16_serostatus.resize(population_size, false);
    number_of_treatments.resize(population_size, 0);
    time_of_last_treatment.resize(population_size, -1);
    compliance.resize(population_size, 0);
    diagnostic_random_vals_for_timestep.resize(population_size, 1);
    serorevert_fast.resize(population_size, false);

    // Age burn-in
    std::vector<double> curr_age(population_size, 0.0);
    for (int i = 0; i < 75000; ++i) {
        for (int j = 0; j < population_size; ++j) {
            curr_age[j] += timestep_years;
            if (death_dist(gen) || curr_age[j] > max_age)
                curr_age[j] = 0.0;
        }
    }
    ages = curr_age;

    // set sex, serorevert fast, and exposure heterogeneity
    for (int i = 0; i < population_size; ++i) {
        sex[i] = gender_dist(gen);
        exposure_heterogeneity[i] = gamma_dist(gen);
        serorevert_fast[i] = serorevert_fast_dist(gen);
    }

    // All four worm types share the same compartment structure, movement between them differ
    fertile_female_worms = WormPopulation(population_size, worm_params, worm_params.initial_worms);
    infertile_female_worms = WormPopulation(population_size, worm_params, worm_params.initial_worms);
    sterilized_female_worms = WormPopulation(population_size, worm_params, worm_params.initial_worms);
    male_worms = WormPopulation(population_size, worm_params, worm_params.initial_worms);
    fertile_female_worms.update_death_dists(timestep_years, nullptr);
    infertile_female_worms.update_death_dists(timestep_years, nullptr);
    sterilized_female_worms.update_death_dists(timestep_years, nullptr);
    male_worms.update_death_dists(timestep_years, nullptr);

    microfilariae = MFPopulation(population_size, mf_params, mf_params.initial_mf);

    l3 = L3Population(population_size, 
        (int)std::max(
            std::round(worm_params.l3_delay / delta_time_days),
            1.0
        )
    );

    blackflies = BlackflyPopulation(
        population_size, 
        blackfly_params,
        exposure_heterogeneity,
        delta_time_days
    );

    // Temporary storage buffers, used to reduce number of times space needs to be allocated
    const int worm_comps = worm_params.worm_age_stages;
    temp_to_fertile.assign(population_size * worm_comps, 0.0);
    temp_to_infertile.assign(population_size * worm_comps, 0.0);
    temp_zeros.assign(population_size, 0.0);
    temp_mf_loads.assign(population_size, 0.0);

    // Initialize Sequelae
    int num_sequelae = sequelae_params.size();
    sequelae.resize(num_sequelae);

    for (int s = 0; s < num_sequelae; ++s) {
        switch(sequelae_params[s]->sequelae_model_type) {
            case SequelaeModelType::TimestepProb: {
                const TimestepProbSequelaeParams& tsp = static_cast<TimestepProbSequelaeParams&>(*sequelae_params[s]);
                sequelae[s] = std::make_unique<TimestepProbSequelae>(
                    TimestepProbSequelae(
                        tsp.sequelae_type, tsp.sequelae_model_type, population_size,
                        tsp.base_probability, tsp.prob_timescale, tsp.min_age_test,
                        tsp.min_infection, tsp.retest_tested_indivs, tsp.countdown_timesteps,
                        tsp.status_end_countdown, tsp.use_raw_infection_for_test,
                        tsp.average_age
                    )
                );
                break;
            }
            case SequelaeModelType::ExponentialProb: {
                const ExponentialProbSequelaeParams& esp = static_cast<ExponentialProbSequelaeParams&>(*sequelae_params[s]);
                sequelae[s] = std::make_unique<ExponentialProbSequelae>(
                    ExponentialProbSequelae(
                        esp.sequelae_type, esp.sequelae_model_type, population_size,
                        esp.base_probability, esp.prob_timescale, esp.min_age_test,
                        esp.min_infection, esp.retest_tested_indivs, esp.countdown_timesteps,
                        esp.status_end_countdown, esp.use_raw_infection_for_test,
                        esp.prob_intercept, esp.prob_slope
                    )
                );
                break;
            }
            case SequelaeModelType::PowerLawProb: {
                const PowerLawProbSequelaeParams& psp = static_cast<PowerLawProbSequelaeParams&>(*sequelae_params[s]);
                sequelae[s] = std::make_unique<PowerLawProbSequelae>(
                    PowerLawProbSequelae(
                        psp.sequelae_type, psp.sequelae_model_type, population_size,
                        psp.base_probability, psp.prob_timescale, psp.min_age_test,
                        psp.min_infection, psp.retest_tested_indivs, psp.countdown_timesteps,
                        psp.status_end_countdown, psp.use_raw_infection_for_test, 
                        psp.prob_intercept, psp.prob_slope
                    )
                );
                break;
            }
            case SequelaeModelType::OAE: {
                const OAESequelaeParams& osp = static_cast<OAESequelaeParams&>(*sequelae_params[s]);
                sequelae[s] = std::make_unique<OAESequelae>(
                    OAESequelae(
                        osp.sequelae_type, osp.sequelae_model_type, population_size,
                        osp.base_probability, osp.prob_timescale, osp.min_age_test,
                        osp.min_infection, osp.retest_tested_indivs, osp.countdown_timesteps,
                        osp.status_end_countdown, osp.use_raw_infection_for_test,
                        osp.prob_intercept, osp.prob_slope,
                        gen, osp.max_age_test, ages
                    )
                );
                break;
            }
            default: {
                const SequelaeParams& sp = *sequelae_params[s];
                sequelae[s] = std::make_unique<Sequelae>(
                    Sequelae(
                        sp.sequelae_type, sp.sequelae_model_type, population_size,
                        sp.base_probability, sp.prob_timescale, sp.min_age_test,
                        sp.min_infection, sp.retest_tested_indivs, sp.countdown_timesteps,
                        sp.status_end_countdown, sp.use_raw_infection_for_test
                    )
                );
            }
        }
    }
}

void People::generate_random_vals_for_diagnostic(std::mt19937& gen) {
    for (int ind = 0; ind < population_size; ++ind) {
        diagnostic_random_vals_for_timestep[ind] = uniform_dist(gen);
    }
}

double People::calculate_individual_compliance(std::mt19937& gen) {
    double alpha = _current_cov * (1.0 - _current_rho) / _current_rho;
    double beta = (1.0 - _current_cov) * (1.0 - _current_rho) / _current_rho;
    return beta_distribution(gen, alpha, beta);
}

void People::calculate_population_compliance(std::mt19937& gen, double rho, double total_population_coverage, double proportion_never_treated) {
    _current_rho = rho;
    _current_cov = total_population_coverage;
    std::bernoulli_distribution dist(proportion_never_treated);
    for (auto& individual : compliance) {
        if (proportion_never_treated > 0 && dist(gen))
            individual = -1;
        else
            individual = calculate_individual_compliance(gen);
        if (individual == 0)
            printf("have a 0 compliance guy");
    }
}

void People::update_compliance(std::mt19937& gen, double rho, double total_population_coverage) {
    _current_rho = rho;
    _current_cov = total_population_coverage;

    // Individuals who are marked as "never treated" don't get a new compliance value
    std::vector<int> valid_indices;
    for (std::size_t i = 0; i < compliance.size(); ++i)
        if (compliance[i] != -1)
            valid_indices.push_back(i);

    int n = valid_indices.size();
    std::vector<double> new_compliance(n);
    for (auto& val : new_compliance)
        val = calculate_individual_compliance(gen);
    std::sort(new_compliance.begin(), new_compliance.end());

    // Individuals who have a higher old compliance should have the higher
    // new compliance values
    std::vector<int> rank_indices(n);
    std::iota(rank_indices.begin(), rank_indices.end(), 0);
    std::sort(rank_indices.begin(), rank_indices.end(),
        [&](int a, int b) { return compliance[valid_indices[a]] < compliance[valid_indices[b]]; }
    );
    for (int i = 0; i < n; ++i)
        compliance[valid_indices[rank_indices[i]]] = new_compliance[i];
}

void People::apply_treatment_round(std::mt19937& gen, int minimum_age_of_treatment, int current_time_step) {
    for (int i = 0; i < population_size; ++i) {
        bool treatable_age = ages[i] >= minimum_age_of_treatment;
        bool passed_compliance = uniform_dist(gen) < compliance[i];
        bool treated = treatable_age && passed_compliance;
        number_of_treatments[i] += treated;
        time_of_last_treatment[i] = treated ? current_time_step : time_of_last_treatment[i];
    }
}

std::vector<double> People::get_exposure(
    const ExposureParams& exp_param,
    double sex_ratio
) const {
    const double exp_females = 1.0 / (sex_ratio * (exp_param.Q - 1.0) + 1.0);
    const double exp_males = exp_param.Q * exp_females;

    const int n_females = (int)std::accumulate(sex.begin(), sex.end(), 0);
    const int n_males = population_size - n_females;

    assert(n_males > 0 && n_females > 0);

    std::vector<double> age_sex_exp(population_size);
    double tmp_alpha_m = 0.0, tmp_alpha_f = 0.0;

    for (int i = 0; i < population_size; ++i) {
        if (!sex[i]) {
            age_sex_exp[i] = std::exp(-exp_param.male_exposure_exponent * ages[i]);
            tmp_alpha_m += age_sex_exp[i];
        } else {
            age_sex_exp[i] = std::exp(-exp_param.female_exposure_exponent * ages[i]);
            tmp_alpha_f += age_sex_exp[i];
        }
    }

    const double male_scale = exp_males / (tmp_alpha_m / n_males);
    const double female_scale = exp_females / (tmp_alpha_f / n_females);

    for (int i = 0; i < population_size; ++i)
        age_sex_exp[i] *= sex[i] ? female_scale : male_scale;

    const double mean_het = (
        std::accumulate(exposure_heterogeneity.begin(), exposure_heterogeneity.end(), 0.0)
        / population_size
    );
    
    for (int i = 0; i < population_size; ++i)
        age_sex_exp[i] *= exposure_heterogeneity[i] / mean_het;

    const double mean_exp = (
        std::accumulate(age_sex_exp.begin(), age_sex_exp.end(), 0.0)
        / population_size
    );
    for (double& v : age_sex_exp)
        v /= mean_exp;

    return age_sex_exp;
}

std::vector<double> People::get_new_worms() {
    std::vector<double> new_worms(population_size);
    for (int i = 0; i < population_size; ++i)
        new_worms[i] = l3.get_new_worms(i);
    return new_worms;
}

double People::wplus1_rate(
    double mean_l3, double delta_h_zero, double delta_h_inf,
    double c_h, double ABR, double exposure,
    double timestep_years
) const {
    const double delta_h = (
        (delta_h_zero + delta_h_inf * c_h * ABR * mean_l3 * exposure)
        / (1.0 + c_h * ABR * mean_l3 * exposure)
    );
    return timestep_years * ABR * delta_h * exposure * mean_l3;
}

void People::new_established_l3(
    std::mt19937& gen,
    double delta_h_zero,
    double delta_h_inf,
    double c_h,
    double ABR,
    const std::vector<double>& exposure,
    double timestep_years
) {
    const double mean_l3 = blackflies.mean_l3();
    for (int i = 0; i < population_size; ++i) {
        const double rate = wplus1_rate(
            mean_l3, delta_h_zero, delta_h_inf,
            c_h, ABR, exposure[i], timestep_years
        );
        const int new_l3 = std::poisson_distribution<int>(rate)(gen);
        l3.add_new_established_l3(i, new_l3);
    }
}

void People::update_ov16_status_individual(int indiv_index) {
    bool has_male_worm = male_worms.get_raw_load(indiv_index) > 0;
    bool has_fertile_female_worm = fertile_female_worms.get_raw_load(indiv_index) > 0;
    bool has_infertile_female_worm = infertile_female_worms.get_raw_load(indiv_index) > 0;
    bool has_any_mf = microfilariae.get_raw_load(indiv_index) > 0;
    bool has_any_l3 = l3.get_raw_load(indiv_index) > 0;

    bool seropositive = (has_male_worm && has_fertile_female_worm) && has_any_mf;
    bool slow_seroreversion = !has_male_worm && !has_fertile_female_worm && !has_infertile_female_worm && !has_any_l3;
    bool fast_seroreversion = !has_male_worm || !has_fertile_female_worm || !has_any_mf;
    bool seroreversion_condition_to_use = serorevert_fast[indiv_index] ? fast_seroreversion : slow_seroreversion;

    ov16_serostatus[indiv_index] = !ov16_serostatus[indiv_index] ? seropositive : ov16_serostatus[indiv_index];

    ov16_serostatus[indiv_index] = ov16_serostatus[indiv_index] ? !seroreversion_condition_to_use : ov16_serostatus[indiv_index];
}

void People::update_all_status(
    std::mt19937& generator,
    double timestep_years, int timesteps_per_year,
    int skin_snip_weight, int num_skin_snips
) {
    std::vector<int> skin_snip_count(population_size);
    std::vector<int> raw_mf_count(population_size);
    bool mating_worm_pair = false;
    bool male_female_worm_pair = false;
    for (int i = 0; i < population_size; ++i) {
        update_ov16_status_individual(i);
        skin_snip_count[i] = round(microfilariae.get_skin_snip_load_person(generator, i, skin_snip_weight, num_skin_snips));
        raw_mf_count[i] = round(microfilariae.get_raw_load(i));
        bool has_male_worm = male_worms.get_raw_load(i) > 0;
        bool has_fertile_female_worm = fertile_female_worms.get_raw_load(i) > 0;
        bool has_infertile_female_worm = infertile_female_worms.get_raw_load(i) > 0;

        mating_worm_pair = has_male_worm && has_fertile_female_worm;
        male_female_worm_pair = has_male_worm && (has_fertile_female_worm || has_infertile_female_worm);
    }
    for (auto& s : sequelae) {
        s->update_all_individuals(
            generator, uniform_dist, timestep_years, 
            timesteps_per_year, population_size, ages,
            raw_mf_count, skin_snip_count,
            mating_worm_pair, male_female_worm_pair
        );
    }

}

bool People::sample_serostatus_individual(int indiv_index, double sens, double spec) {
    return (
        ov16_serostatus[indiv_index] ?
        diagnostic_random_vals_for_timestep[indiv_index] <= sens :
        diagnostic_random_vals_for_timestep[indiv_index] > spec
    );
}

static double get_elapsed_time(clock_t start_time) {
    return (double)(clock() - start_time) / CLOCKS_PER_SEC;
}

void People::update_worm_dists(double timestep_years) {
    DrugParams* drug_params = nullptr;
    if (treatment_params.has_value()) {
        drug_params = &treatment_params->drug_params;
    }
    fertile_female_worms.update_death_dists(timestep_years, drug_params);
    infertile_female_worms.update_death_dists(timestep_years, drug_params);
    sterilized_female_worms.update_death_dists(timestep_years, drug_params);
    male_worms.update_death_dists(timestep_years, drug_params);
}

void People::age(
    std::mt19937& gen, 
    int current_timestep, double timestep_years,
    const std::vector<double>& new_male_worms,
    const std::vector<double>& new_female_worms,
    std::vector<double>& age_timers
) {
    clock_t start = clock();
    DrugParams* drug_params = nullptr;
    if (treatment_params.has_value()) {
        drug_params = &treatment_params->drug_params;
    }
    age_timers[0] += get_elapsed_time(start);
    start = clock();
    
    // Age MF
    microfilariae.age(
        gen,
        false, // stochastic
        current_timestep, timestep_years, male_worms, fertile_female_worms,
        drug_params, number_of_treatments, time_of_last_treatment
    );
    age_timers[1] += get_elapsed_time(start);
    start = clock();

    std::fill(temp_to_fertile.begin(), temp_to_fertile.end(), 0.0);
    std::fill(temp_to_infertile.begin(), temp_to_infertile.end(), 0.0);
    for (int p = 0; p < population_size; ++p) {
        male_worms.age_single(
            p, gen, WormType::Male, current_timestep, timestep_years, new_male_worms,
            drug_params, number_of_treatments, time_of_last_treatment
        );
        infertile_female_worms.age_single_swap(
            p, gen, WormType::Infertile_Female, current_timestep, timestep_years, new_female_worms, temp_to_fertile,
            drug_params, number_of_treatments, time_of_last_treatment
        );
        fertile_female_worms.age_single_swap(
            p, gen, WormType::Fertile_Female, current_timestep, timestep_years, temp_zeros, temp_to_infertile,
            drug_params, number_of_treatments, time_of_last_treatment
        );
    }
    age_timers[4] += get_elapsed_time(start);
    start = clock();

    // Apply swapped worms to their destination populations
    fertile_female_worms.age_helper_swapped_worms(temp_to_fertile);
    infertile_female_worms.age_helper_swapped_worms(temp_to_infertile);
    age_timers[5] += get_elapsed_time(start);
    start = clock();

    // Advance L3 delay buffer
    l3.age_all();
    age_timers[6] += get_elapsed_time(start);
    start = clock();

    // Increment ages
    for (double& a : ages)
        a += timestep_years;
    age_timers[7] += get_elapsed_time(start);
    start = clock();
}

void People::process_deaths(std::mt19937& gen) {
    for (int i = 0; i < population_size; ++i) {
        bool is_dead = death_dist(gen) || (ages[i] > max_age);
        if (!is_dead) continue;

        // Reset demographics
        ages[i] = 0.0;
        sex[i] = gender_dist(gen);
        exposure_heterogeneity[i] = gamma_dist(gen);
        ov16_serostatus[i] = false;
        serorevert_fast[i] = serorevert_fast_dist(gen);

        // Reset any treatment_values
        if (_current_rho != -1 && _current_cov != -1) {
            compliance[i] = calculate_individual_compliance(gen);
            number_of_treatments[i] = 0;
            time_of_last_treatment[i] = -1;
        }

        // Reset all parasite populations for this person
        fertile_female_worms.process_death(i);
        infertile_female_worms.process_death(i);
        sterilized_female_worms.process_death(i);
        male_worms.process_death(i);
        microfilariae.process_death(i);
        l3.process_death(i);
        blackflies.process_death(i);
        for (auto& s : sequelae)
            s->process_death(i);
    }
}

double People::mean_l1_per_blackfly() const { return blackflies.mean_l1(); }
double People::mean_l2_per_blackfly() const { return blackflies.mean_l2(); }
double People::mean_l3_per_blackfly() const { return blackflies.mean_l3(); }
double People::mean_l3_prevalence_blackflies()  const { return blackflies.mean_l3_prevalence(); }

std::vector<int> People::get_sub_population(int age_start, int age_end) {
    std::vector<int> indeces_of_sub_pop;
    for (int i = 0; i < population_size; ++i) {
        if (age_start <= ages[i] && ages[i] < age_end) {
            indeces_of_sub_pop.push_back(i);
        }
    }
    return indeces_of_sub_pop;
}

// Cloning helper functions
People::People(const People& other)
    : population_size(other.population_size),
      mean_age(other.mean_age),
      max_age(other.max_age),
      gender_ratio(other.gender_ratio),
      prop_serorevert_fast(other.prop_serorevert_fast),
      _current_rho(other._current_rho),
      _current_cov(other._current_cov),
      gamma_dist(other.gamma_dist),
      gender_dist(other.gender_dist),
      death_dist(other.death_dist),
      serorevert_fast_dist(other.serorevert_fast_dist),
      uniform_dist(other.uniform_dist),
      ages(other.ages),
      sex(other.sex),
      exposure_heterogeneity(other.exposure_heterogeneity),
      ov16_serostatus(other.ov16_serostatus),
      number_of_treatments(other.number_of_treatments),
      time_of_last_treatment(other.time_of_last_treatment),
      compliance(other.compliance),
      serorevert_fast(other.serorevert_fast),
      fertile_female_worms(other.fertile_female_worms),
      infertile_female_worms(other.infertile_female_worms),
      sterilized_female_worms(other.sterilized_female_worms),
      male_worms(other.male_worms),
      microfilariae(other.microfilariae),
      l3(other.l3),
      blackflies(other.blackflies),
      treatment_params(other.treatment_params),
      diagnostic_random_vals_for_timestep(other.diagnostic_random_vals_for_timestep),
      temp_to_fertile(other.temp_to_fertile),
      temp_to_infertile(other.temp_to_infertile),
      temp_zeros(other.temp_zeros),
      temp_mf_loads(other.temp_mf_loads)
{
    sequelae.reserve(other.sequelae.size());
    for (const auto& s : other.sequelae)
        sequelae.push_back(s->clone());
}

People& People::operator=(const People& other) {
    if (this == &other) return *this;
    population_size = other.population_size;
    mean_age = other.mean_age;
    max_age = other.max_age;
    gender_ratio = other.gender_ratio;
    prop_serorevert_fast = other.prop_serorevert_fast;
    _current_rho = other._current_rho;
    _current_cov = other._current_cov;
    gamma_dist = other.gamma_dist;
    gender_dist = other.gender_dist;
    death_dist = other.death_dist;
    serorevert_fast_dist = other.serorevert_fast_dist;
    uniform_dist = other.uniform_dist;
    ages = other.ages;
    sex = other.sex;
    exposure_heterogeneity = other.exposure_heterogeneity;
    ov16_serostatus = other.ov16_serostatus;
    number_of_treatments = other.number_of_treatments;
    time_of_last_treatment = other.time_of_last_treatment;
    compliance = other.compliance;
    serorevert_fast = other.serorevert_fast;
    fertile_female_worms = other.fertile_female_worms;
    infertile_female_worms = other.infertile_female_worms;
    sterilized_female_worms = other.sterilized_female_worms;
    male_worms = other.male_worms;
    microfilariae = other.microfilariae;
    l3 = other.l3;
    blackflies = other.blackflies;
    treatment_params = other.treatment_params;
    diagnostic_random_vals_for_timestep = other.diagnostic_random_vals_for_timestep;
    temp_to_fertile = other.temp_to_fertile;
    temp_to_infertile = other.temp_to_infertile;
    temp_zeros = other.temp_zeros;
    temp_mf_loads = other.temp_mf_loads;
    sequelae.clear();
    sequelae.reserve(other.sequelae.size());
    for (const auto& s : other.sequelae)
        sequelae.push_back(s->clone());
    return *this;
}
