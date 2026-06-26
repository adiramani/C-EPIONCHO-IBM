#ifndef PARASITE_HPP
#define PARASITE_HPP

#include <vector>
#include <random>
#include "params.hpp"

class ParasitePopulation {
public:
    int n_people = 0;
    int compartments = 0;
    float years_per_compartment = 0.0f;
    float y_mortality_rate = 0.0f;
    float d_mortality_rate = 0.0f;

    // Stores all n_people × compartments values in a single flat vector:
    //   parasites[person * compartments + compartment]
    std::vector<double> parasites;

    // Shared age-category lookup (same for every person)
    std::vector<double> age_categories;  // [compartments]

    ParasitePopulation() = default;
    ParasitePopulation(
        int n_people_,
        int compartments_,
        float years_per_compartment_,
        float y_mortality_rate_,
        float d_mortality_rate_,
        double initial_value
    );

    // Convenience indexed accessors
    inline double& at(int person, int c) {
        return parasites[person * compartments + c];
    }
    inline double at(int person, int c) const {
        return parasites[person * compartments + c];
    }

    double weibull_mortality(int c, double timestep_years) const;

    // Zero all compartments for one person (on birth/death)
    virtual void process_death(int person_idx);

    double get_raw_load(int person_idx) const;

    void get_all_raw_loads(std::vector<double>& out) const;
};


class WormPopulation : public ParasitePopulation {
public:
    float mf_production_rate = 0.0f;  // epsilon
    float infertile_to_fertile_rate = 0.0f;  // omega
    float fertile_to_infertile_rate = 0.0f;  // lambda_zero
    float F = 0.0f;
    float G = 0.0f;
    std::vector<double> death_dists;
    std::vector<double> death_dists_permanent_sterilization;
    double ageing_dist;

    WormPopulation() = default;
    WormPopulation(int n_people_, const WormParams& params, int initial_worms);

    double fecundity_rate(int c) const;

    // Returns how many worms transition between fertile/infertile states
    int fecundity_movement(
        double alive_worms,
        double treatment_induced_temp_sterility,
        std::mt19937& gen
    ) const;

    void update_death_dists(double timestep_years, DrugParams* drug_params);

    // Ages every person in one pass.
    // new_worms[i] = newly established L3 entering person i this timestep
    // swapped_out = if non-null, size must be n_people * compartments;
    // filled with per-person per-compartment counts that
    // switch state (fertile↔infertile)
    void age(
        std::mt19937& gen,
        WormType type,
        int current_timestep,
        double timestep_years,
        const std::vector<double>& new_worms,
        std::vector<double>* swapped_out,
        DrugParams* drug_params,
        std::vector<int>& number_of_treatments,
        std::vector<int>& time_of_last_treatment
    );

    void age_single(
        int person_index,
        std::mt19937& gen,
        WormType type,
        int current_timestep,
        double timestep_years,
        const std::vector<double>& new_worms,
        DrugParams* drug_params,
        std::vector<int>& number_of_treatments,
        std::vector<int>& time_of_last_treatment
    );

    void age_single_swap(
        int person_index,
        std::mt19937& gen,
        WormType type,
        int current_timestep,
        double timestep_years,
        const std::vector<double>& new_worms,
        std::vector<double>& swapped_out,
        DrugParams* drug_params,
        std::vector<int>& number_of_treatments,
        std::vector<int>& time_of_last_treatment
    );

    // Add incoming worms that switched from the complementary population
    void age_helper_swapped_worms(const std::vector<double>& incoming);

    void process_death(int person_idx) override;
};

class MFPopulation : public ParasitePopulation {
public:
    double mf_move_rate = 0.0;
    double kmf_const;
    bool use_kmf_const;
    double initial_kmf;
    double slope_kmf;

    MFPopulation() = default;
    MFPopulation(int n_people_, const MicrofilariaeParams& params, double initial_mf);

    void age(
        std::mt19937& gen,
        bool stochastic,
        int current_timestep,
        double timestep_years,
        const WormPopulation& male_worms,
        const WormPopulation& fertile_female_worms,
        DrugParams* drug_params,
        std::vector<int>& number_of_treatments,
        std::vector<int>& time_of_last_treatment
    );

    double get_skin_snip_load_person(std::mt19937& gen, int person_idx, int skin_snip_weight, int num_skin_snips);

private:
    int calc_new_mf_for_person(
        std::mt19937& gen,
        bool stochastic,
        int person,
        double timestep_years,
        double treatment_microfilaricidal_effect,
        double male_load,
        const WormPopulation& ff_worms
    );

    int age_exiting_mf_for_person(
        std::mt19937& gen,
        bool stochastic,
        int person, int c, double timestep_years, 
        double prev_mf, double treatment_microfilaricidal_effect
    );
};

class L3Population {
public:
    int n_people = 0;
    int compartments = 0;
    int current_index = 0;

    // index = person * compartments + slot
    std::vector<double> parasites;

    L3Population() = default;
    L3Population(int n_people_, int compartments_);

    // Advance the shared circular index (call once per timestep)
    void age_all();

    double get_raw_load(int person_idx);

    double get_new_worms(int person_idx);

    // Write newly established L3 into the current delay slot for one person
    void add_new_established_l3(int person_idx, int num_l3);

    void process_death(int person_idx);
};

#endif
