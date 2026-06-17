#include "parasite.hpp"
#include <numeric>
#include <cmath>
#include <algorithm>
#include <cassert>

ParasitePopulation::ParasitePopulation(
    int n_people_,
    int compartments_,
    float ypc,
    float y,
    float d,
    double initial_value
) : n_people(n_people_),
    compartments(compartments_),
    years_per_compartment(ypc),
    y_mortality_rate(y),
    d_mortality_rate(d),
    parasites(n_people_ * compartments_, initial_value),
    age_categories(compartments_)
{
    for (int c = 0; c < compartments_; ++c)
        age_categories[c] = c * (double)ypc;
}

double ParasitePopulation::weibull_mortality(int c, double timestep_years) const {
    return timestep_years
         * std::pow((double)y_mortality_rate, (double)d_mortality_rate)
         * (double)d_mortality_rate
         * std::pow(age_categories[c], (double)d_mortality_rate - 1.0);
}

void ParasitePopulation::process_death(int person_idx) {
    int base = person_idx * compartments;
    std::fill(
        parasites.begin() + base,
        parasites.begin() + base + compartments, 0.0
    );
}

double ParasitePopulation::get_raw_load(int person_idx) const {
    int base = person_idx * compartments;
    return std::accumulate(
        parasites.begin() + base,
        parasites.begin() + base + compartments, 0.0
    );
}

void ParasitePopulation::get_all_raw_loads(std::vector<double>& out) const {
    out.resize(n_people);
    for (int p = 0; p < n_people; ++p)
        out[p] = get_raw_load(p);
}


WormPopulation::WormPopulation(int n_people_, const WormParams& params, int initial_worms)
    : ParasitePopulation(
        n_people_,
        params.worm_age_stages,
        (float)params.q_w,
        (float)params.y_w,
        (float)params.d_w,
        (double)initial_worms
    ),
    mf_production_rate((float)params.epsilon),
    infertile_to_fertile_rate((float)params.omega),
    fertile_to_infertile_rate((float)params.lambda_zero),
    F((float)params.fecundity_worms_F),
    G((float)params.fecundity_worms_G)
{}

double WormPopulation::fecundity_rate(int c) const {
    return (
        mf_production_rate
        * F
        / (F + std::pow((double)G, -1.0 * age_categories[c]) - 1.0)
    );
}

int WormPopulation::fecundity_movement(
    double alive_worms,
    double sterility_swap_rate,
    std::mt19937& gen
) const {
    int n = std::max(0, (int)alive_worms);
    if (n == 0) return 0;

    std::binomial_distribution<int> dist(n, sterility_swap_rate);// infertile_to_fertile_rate * timestep_years);
    return dist(gen);
}

void WormPopulation::age(
    std::mt19937& gen,
    WormType type,
    double current_timestep,
    double timestep_years,
    const std::vector<double>& new_worms,
    std::vector<double>* swapped_out,
    DrugParams* drug_params,
    std::vector<int>& number_of_treatments,
    std::vector<int>& time_of_last_treatment
) {
    const double age_rate = timestep_years / (double)years_per_compartment;
    bool is_fertile_female = type == WormType::Fertile_Female;
    bool is_fertile_or_tmp_infertile_female = is_fertile_female || (type == WormType::Infertile_Female);

    // Pre-populate vectors for distributions. Using pre-computed bernoulli distributions vs
    // creating new binomial distributions for each timestep saves 50% run time for the whole model
    // TODO: since they only change if timestep changes, possibly move the distributions to the population
    // class.
    std::vector<std::bernoulli_distribution> death_dists;
    std::vector<std::bernoulli_distribution> death_dists_permanent_sterilization;
    death_dists.resize(compartments);
    death_dists_permanent_sterilization.resize(compartments);
    for (int c = 0; c < compartments; ++c) {
        double base_mortality = weibull_mortality(c, timestep_years);
        death_dists[c] = std::bernoulli_distribution(base_mortality);
        if (drug_params) {
            death_dists_permanent_sterilization[c] = std::bernoulli_distribution(
                std::min(
                    base_mortality + drug_params->permanent_infertility,
                    1.0
                )
            );
        }
    }
    std::bernoulli_distribution ageing_dist(age_rate);

    for (int person = 0; person < n_people; ++person) {
        const int base = person * compartments;
        double incoming = new_worms[person];
        
        // Effects of Treatment
        bool was_treated = number_of_treatments[person] > 0;
        double time_since_treatment_years = was_treated ? (current_timestep - time_of_last_treatment[person]) * timestep_years : -1.0;
        
        double treatment_induced_temp_sterility = (was_treated && is_fertile_female) ? drug_params->embryostatic_lambda_max * exp(-drug_params->embryostatic_phi * time_since_treatment_years) : 0;
        bool apply_permanent_infertility = ((number_of_treatments[person] >= 2) && time_since_treatment_years == 0 && is_fertile_or_tmp_infertile_female);

        for (int c = 0; c < compartments; ++c) {
            const int n_worms = std::max(0, (int)std::round(parasites[base + c]));
            if (n_worms <= 0) {
                parasites[base + c] += incoming;
                incoming = 0;
                if (swapped_out) {
                    (*swapped_out)[base + c] = 0;
                }
                continue;
            }
            double dead_worms = 0.0;
            for(int w = 0; w < n_worms; ++w) {
                dead_worms += apply_permanent_infertility ? death_dists_permanent_sterilization[c](gen) : death_dists[c](gen);
            }


            const int n_after_death = std::max(0, n_worms - (int)dead_worms);
            double aging_worms = 0.0;
            if (n_after_death > 0) {
                for (int nad = 0; nad < n_after_death; ++nad) {
                    aging_worms += ageing_dist(gen);
                }
            }

            double outgoing_swapped = 0.0;
            if (type != WormType::Male && type != WormType::Sterilized_Female) {
                const double current_alive = parasites[base + c] - dead_worms - aging_worms;
                double sterility_swap_rate = (
                    type == WormType::Fertile_Female ? 
                    (fertile_to_infertile_rate * timestep_years) + treatment_induced_temp_sterility :
                    (infertile_to_fertile_rate * timestep_years)
                );
                outgoing_swapped = (double)fecundity_movement(
                    current_alive, std::min(sterility_swap_rate, 1.0), gen
                );
            }

            parasites[base + c] += incoming - dead_worms - aging_worms - outgoing_swapped;
            incoming = aging_worms;

            if (swapped_out) {
                (*swapped_out)[base + c] = outgoing_swapped;
            }
        }
    }
}

void WormPopulation::age_helper_swapped_worms(const std::vector<double>& incoming) {
    // incoming has a flat layout: [person * compartments + c]
    const int total = n_people * compartments;
    for (int i = 0; i < total; ++i)
        parasites[i] += incoming[i];
}

static double derivmf_one(
    int male_present,
    const WormPopulation& ff_worms, int person,
    double mf_current,
    double mf_mort, double mf_move, double k_in
) {
    // Sum fecundity contribution across all worm compartments for this person
    double new_incoming = 0.0;
    const int base = person * ff_worms.compartments;
    for (int c = 0; c < ff_worms.compartments; ++c)
        new_incoming += ff_worms.parasites[base + c] * ff_worms.fecundity_rate(c);

    const double die = std::max(mf_mort * (mf_current + k_in), 0.0);
    const double move_to_next = std::max(mf_move * (mf_current + k_in), 0.0);
    return (double)male_present * new_incoming - die - move_to_next;
}

static double derivmf_rest(
    double mf_mort, double mf_move,
    double mf_current, double mf_previous, double k_in
) {
    const double prev_to_curr = std::max(mf_previous * mf_move, 0.0);
    const double die = std::max(mf_mort * (mf_current + k_in), 0.0);
    const double move_to_next = std::max(mf_move * (mf_current + k_in), 0.0);
    return prev_to_curr - die - move_to_next;
}

MFPopulation::MFPopulation(
    int n_people_, const MicrofilariaeParams& params, double initial_mf
) : ParasitePopulation(n_people_,
        params.mf_age_stages,
        (float)params.q_m,
        (float)params.y_m,
        (float)params.d_m,
        initial_mf
    ),
    mf_move_rate(params.mf_move_rate),
    kmf_const(params.kmf_const),
    use_kmf_const(params.use_kmf_const),
    initial_kmf(params.initial_kmf),
    slope_kmf(params.slope_kmf)
{}

void MFPopulation::calc_new_mf_for_person(
    int person,
    double timestep_years,
    double treatment_microfilaricidal_effect,
    double male_load,
    const WormPopulation& ff_worms
) {
    const int male_present = (male_load > 0.0) ? 1 : 0;
    const double mortality = weibull_mortality(0, 1.0) + treatment_microfilaricidal_effect;
    double& mf0 = parasites[person * compartments + 0];

    const double k1 = derivmf_one(
        male_present, ff_worms, person, mf0,
        mortality, mf_move_rate, 0.0
    );
    const double k2 = derivmf_one(
        male_present, ff_worms, person, mf0,
        mortality, mf_move_rate, timestep_years * k1 / 2.0
    );
    const double k3 = derivmf_one(
        male_present, ff_worms, person, mf0,
        mortality, mf_move_rate, timestep_years * k2 / 2.0
    );
    const double k4 = derivmf_one(
        male_present, ff_worms, person, mf0,
        mortality, mf_move_rate, timestep_years * k3
    );

    mf0 += (timestep_years / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
}

void MFPopulation::age_exiting_mf_for_person(
    int person, int c,
    double timestep_years, double prev_mf,
    double treatment_microfilaricidal_effect
) {
    const double mortality = weibull_mortality(c, 1.0) + treatment_microfilaricidal_effect;
    double& mfc = parasites[person * compartments + c];

    const double k1 = derivmf_rest(mortality, mf_move_rate, mfc, prev_mf, 0.0);
    const double k2 = derivmf_rest(mortality, mf_move_rate, mfc, prev_mf, timestep_years * k1 / 2.0);
    const double k3 = derivmf_rest(mortality, mf_move_rate, mfc, prev_mf, timestep_years * k2 / 2.0);
    const double k4 = derivmf_rest(mortality, mf_move_rate, mfc, prev_mf, timestep_years * k3);

    mfc += (timestep_years / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
}

void MFPopulation::age(
    int current_timestep,
    double timestep_years,
    const WormPopulation& male_worms,
    const WormPopulation& fertile_female_worms,
    DrugParams* drug_params,
    std::vector<int>& number_of_treatments,
    std::vector<int>& time_of_last_treatment
) {
    for (int person = 0; person < n_people; ++person) {
        // Effects of Treatment
        bool was_treated = number_of_treatments[person] > 0;
        int time_since_treatment_years = was_treated ? (current_timestep - time_of_last_treatment[person]) * timestep_years : -1;
        
        double treatment_microfilaricidal_effect = was_treated ? pow(time_since_treatment_years + drug_params->microfilaricidal_upsilon, -drug_params->microfilaricidal_kappa) : 0;
        
        const double male_load = male_worms.get_raw_load(person);
        const int base = person * compartments;

        // Save the old compartment-0 value before RK4 updates it in place,
        // then chain that saved value as 'prev_mf' into compartment 1, etc.
        const double prev_mf_0 = parasites[base + 0];
        calc_new_mf_for_person(person, timestep_years, treatment_microfilaricidal_effect, male_load, fertile_female_worms);

        double prev_mf = prev_mf_0;
        for (int c = 1; c < compartments; ++c) {
            const double curr_mf = parasites[base + c];  // save before update
            age_exiting_mf_for_person(person, c, timestep_years, prev_mf, treatment_microfilaricidal_effect);
            prev_mf = curr_mf;
        }
    }
}

double MFPopulation::get_skin_snip_load_person(
    std::mt19937& gen, int person_idx,
    int skin_snip_weight, int num_skin_snips
) {
    double mf_load = get_raw_load(person_idx);
    if (mf_load == 0.0)
        return 0.0;
    double kmf_to_use = kmf_const;
    // if (!use_kmf_const)...

    double skin_snip_load = 0.0;
    double mu = skin_snip_weight * mf_load;
    double p = kmf_to_use / (kmf_to_use + mu);
    std::negative_binomial_distribution<int> skin_snip_dist(kmf_to_use, p);
    for (int snip = 0; snip < num_skin_snips; ++snip) {
        skin_snip_load += skin_snip_dist(gen);
    }
    return (double)skin_snip_load / (skin_snip_weight * num_skin_snips);
}


L3Population::L3Population(int n_people_, int compartments_)
    : n_people(n_people_),
      compartments(compartments_),
      current_index(0),
      parasites(n_people_ * compartments_, 0.0)
{}

void L3Population::age_all() {
    current_index = (current_index + 1 >= compartments) ? 0 : current_index + 1;
}

double L3Population::get_new_worms(int person_idx) {
    double& slot = parasites[person_idx * compartments + current_index];
    const double val = slot;
    slot = 0.0;
    return val;
}

double L3Population::get_raw_load(int person_idx) {
    int base = person_idx * compartments;
    return std::accumulate(
        parasites.begin() + base,
        parasites.begin() + base + compartments, 0.0
    );
}

void L3Population::add_new_established_l3(int person_idx, int num_l3) {
    parasites[person_idx * compartments + current_index] = (double)num_l3;
}

void L3Population::process_death(int person_idx) {
    int base = person_idx * compartments;
    std::fill(
        parasites.begin() + base,
        parasites.begin() + base + compartments, 0.0
    );
}
