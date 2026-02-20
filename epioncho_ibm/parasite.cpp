#include "parasite.hpp"
#include <numeric>
#include <random>
#include <iostream>

double derivmf_one(int male_present, Worm fetile_female_worms, double mf_current, double mf_mort, double mf_move, double k_in) {
    double new_incoming = 0.0;
    for (long i = 1; i < fetile_female_worms.parasites.size(); ++i){
        new_incoming += fetile_female_worms.parasites[i] * fetile_female_worms.fecundity_rate(i);
    }

    double die = std::max(mf_mort * (mf_current + k_in), 0.0);
    double move_to_next = std::max(mf_move * (mf_current + k_in), 0.0);


    return male_present * new_incoming - die - move_to_next;
}

double derivmf_rest(double mf_mort, double mf_move, double mf_current, double mf_previous, double k_in) {
    double new_current = 0.0;
    double previous_to_current = std::max(mf_previous * mf_move, 0.0);

    double die = std::max(mf_mort * (mf_current + k_in), 0.0);
    double move_to_next = std::max(mf_move * (mf_current + k_in), 0.0);

    return previous_to_current - die - move_to_next;
}

// Parasite constructor
Parasite::Parasite(
    int compartments_,
    float years_per_compartment_,
    float y_mortality_rate_by_age_,
    float d_mortality_rate_by_age_,
    int initial_parasites_
)
    : compartments(compartments_),
      years_per_compartment(years_per_compartment_),
      y_mortality_rate_by_age(y_mortality_rate_by_age_),
      d_mortality_rate_by_age(d_mortality_rate_by_age_),
      parasites(compartments_, initial_parasites_),
      age_categories(compartments_)
{
    for (int i = 0; i < compartments_; ++i) {
        age_categories[i] = i * years_per_compartment;
    }
}

// Worm constructor
Worm::Worm(
    int compartments_,
    float years_per_compartment_,
    float y_mortality_rate_by_age_,
    float d_mortality_rate_by_age_,
    int initial_parasites_,
    float mf_production_rate_,
    float infertile_to_fertile_rate_,
    float fertile_to_infertile_rate_,
    float F_,
    float G_
)
    : Parasite(
          compartments_,
          years_per_compartment_,
          y_mortality_rate_by_age_,
          d_mortality_rate_by_age_,
          initial_parasites_
      ),
      mf_production_rate(mf_production_rate_),
      infertile_to_fertile_rate(infertile_to_fertile_rate_),
      fertile_to_infertile_rate(fertile_to_infertile_rate_),
      F(F_),
      G(G_)
{}


// MF constructor
MF::MF(
    int compartments_,
    float years_per_compartment_,
    float y_mortality_rate_by_age_,
    float d_mortality_rate_by_age_,
    double mf_move_rate_,
    int initial_parasites_
)
    : Parasite(
          compartments_,
          years_per_compartment_,
          y_mortality_rate_by_age_,
          d_mortality_rate_by_age_,
          initial_parasites_
      ),
      mf_move_rate(mf_move_rate_)
{}

L3::L3(
    int compartments_,
    float years_per_compartment_
)
    : Parasite(
        compartments_,
        years_per_compartment_,
        0,
        0,
        0
    ),
    current_index(0)
{}

double Parasite::weibull_mortality(int compartment, double timestep_in_years) {
    return(
        timestep_in_years * std::pow(y_mortality_rate_by_age, d_mortality_rate_by_age) * d_mortality_rate_by_age * std::pow(age_categories[compartment], d_mortality_rate_by_age - 1)
    );
}

void Parasite::process_death() {
    for (auto& compartment : parasites) {
        compartment = 0;
    }
}

double Worm::fecundity_rate(int compartment) {
   return mf_production_rate * F / (F + std::pow(G, -1 * age_categories[compartment]) - 1);
}

std::bernoulli_distribution Worm::fecundity_movement(std::string worm_type, double timestep_years) {
    // TODO: Effects of treatment
    if (worm_type == "fertile_female") {
        return std::bernoulli_distribution(fertile_to_infertile_rate * timestep_years);
    }
    return std::bernoulli_distribution(infertile_to_fertile_rate * timestep_years);
}


std::vector<double> Worm::age(
    std::bernoulli_distribution age_dist, std::mt19937& gen, 
    std::string worm_type, double timestep_years, double new_worms
) {

    // double dead_worms = 0;
    // double aging_worms = 0;
    std::vector<double> female_worms_swapped_out(compartments);

    // for females
    // determine MDA values as needed, update mortality parameters
    // calculate dead fertile, infertile, and sterile
    double incoming_worms = new_worms;
    std::bernoulli_distribution fecundity_movement_gen = fecundity_movement(worm_type, timestep_years);
    for (int i = 0; i < compartments; ++i) {
        double mortality_rate = weibull_mortality(i, timestep_years);
        double dead_worms = 0;
        double aging_worms = 0;
        for (int j = 0; j < parasites[i]; ++j) {
            int is_dead = std::bernoulli_distribution(mortality_rate)(gen);
            dead_worms += is_dead;
            if (is_dead == 0) {
                int did_age = age_dist(gen);
                aging_worms += did_age;
            }
        }
        
        // if (dead_worms > 0) {
        //     std::cout << "Dead: " << dead_worms << " Age " << aging_worms << "\n";
        // }

        double outgoing_worms_swapped = 0.0;
        if (worm_type == "fertile_female" || worm_type == "infertile_female") {
            double current_alive_worms = parasites[i] - dead_worms - aging_worms;
            for (int j = 0; j < current_alive_worms; ++j) {
                outgoing_worms_swapped += fecundity_movement_gen(gen);
            }
        }
        parasites[i] += incoming_worms - dead_worms - aging_worms - outgoing_worms_swapped;
        incoming_worms = aging_worms;
        female_worms_swapped_out[i] = outgoing_worms_swapped;
    }
    return female_worms_swapped_out;
}

void Worm::age_helper_female_swapped_worms(std::vector<double> incoming_swapped_worms) {
    for(int i = 0; i < compartments; ++i) {
        parasites[i] += incoming_swapped_worms[i];
    }
}

void MF::calc_new_mf(double timestep_years, double male_worms, Worm& fertile_female_worms) {
    int male_present = male_worms > 0;
    double mortality = weibull_mortality(0, 1);
    double k1 = derivmf_one(
        male_present, fertile_female_worms, parasites[0],
        mortality, mf_move_rate, 0.0
    );
    double k2 = derivmf_one(
        male_present, fertile_female_worms, parasites[0],
        mortality, mf_move_rate, timestep_years * k1 / 2
    );
    double k3 = derivmf_one(
        male_present, fertile_female_worms, parasites[0],
        mortality, mf_move_rate, timestep_years * k2 / 2
    );
    double k4 = derivmf_one(
        male_present, fertile_female_worms, parasites[0],
        mortality, mf_move_rate, timestep_years * k3
    );

    parasites[0] += (timestep_years / 6) * (k1 + 2 * k2 + 2 * k3 + k4);
}

void MF::age_exiting_mf(int age_category, double timestep_years) {
    double mortality = weibull_mortality(age_category, 1);
    double k1 = derivmf_rest(
        mortality, mf_move_rate,
        parasites[age_category], parasites[age_category-1], 0
    );
    double k2 = derivmf_rest(
        mortality, mf_move_rate,
        parasites[age_category], parasites[age_category-1], timestep_years * k1 / 2
    );
    double k3 = derivmf_rest(
        mortality, mf_move_rate,
        parasites[age_category], parasites[age_category-1], timestep_years * k2 / 2
    );
    double k4 = derivmf_rest(
        mortality, mf_move_rate,
        parasites[age_category], parasites[age_category-1], timestep_years * k3
    );

    parasites[age_category] += (timestep_years / 6) * (k1 + 2 * k2 + 2 * k3 + k4);
}

void MF::age(double timestep_years, double male_worms, Worm& fertile_female_worms) {
    // if there is a male and fertile feamle
    // TODO: impact of MDA

    // RDK4 method for calculating initial mf
    calc_new_mf(timestep_years, male_worms, fertile_female_worms);
    
    // RDK4 method for calculating aging of mf > 1
    for (int i = 1; i < parasites.size(); ++i) {
        age_exiting_mf(i, timestep_years);
    }
}

double Parasite::get_parasite_load() {
    double total = std::accumulate(parasites.begin(), parasites.end(), 0.0);
    return total;
}

void L3::age() {
    if (current_index + 1 >= parasites.size()) {
        current_index = 0;
    } else {
        current_index += 1;
    }
}

void L3::add_new_established_l3(int num_l3) {
    parasites[current_index] = num_l3;
}

double L3::get_new_worms() {
    double number_new_worms = parasites[current_index];
    parasites[current_index] = 0;
    return(
        number_new_worms
    );
}