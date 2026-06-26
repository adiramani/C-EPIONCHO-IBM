#ifndef PEOPLE_HPP
#define PEOPLE_HPP

#include "parasite.hpp"
#include "disease_vector.hpp"
#include "params.hpp"
#include "sequelae.hpp"
#include <vector>
#include <random>
#include <memory>

class People {
public:
    int population_size = 0;
    double mean_age = 0.0;
    double max_age = 0.0;
    double gender_ratio = 0.5;
    double prop_serorevert_fast = 0.5;
    double _current_rho = -1;
    double _current_cov = -1;

    std::gamma_distribution<double> gamma_dist;
    std::bernoulli_distribution gender_dist;
    std::bernoulli_distribution death_dist;
    std::bernoulli_distribution serorevert_fast_dist;
    std::uniform_real_distribution<double> uniform_dist;

    // Tracking status in population
    std::vector<double> ages;
    std::vector<bool> sex; // false=male, true=female
    std::vector<double> exposure_heterogeneity;
    std::vector<bool> ov16_serostatus;
    std::vector<int> number_of_treatments;
    std::vector<int> time_of_last_treatment;
    std::vector<double> compliance;
    std::vector<bool> serorevert_fast;

    // Population objects — one per parasite type
    WormPopulation fertile_female_worms;
    WormPopulation infertile_female_worms;
    WormPopulation sterilized_female_worms;
    WormPopulation male_worms;
    MFPopulation microfilariae;
    L3Population l3;
    BlackflyPopulation blackflies;
    std::vector<std::unique_ptr<Sequelae>> sequelae;

    // Other useful structs to be tracked
    std::optional<TreatmentParams> treatment_params = std::nullopt;
    std::vector<double> diagnostic_random_vals_for_timestep;

    // Buffers — allocated once, reused every call to age().
    // Size = population_size * worm_compartments.
    std::vector<double> temp_to_fertile;
    std::vector<double> temp_to_infertile;
    // Vector of zeros for the fertile-female new_worms argument (no direct L3 input)
    std::vector<double> temp_zeros;
    // Temporary storage buffer for MF loads (populated before blackfly update)
    std::vector<double> temp_mf_loads;

    explicit People(int population_size_);
    People(const People& other);
    People& operator=(const People& other);
    
    void initialize_from_params(
        std::mt19937& gen,
        double timestep_years,
        double delta_time_days,
        double k_e,
        const HumanParams& human_params,
        const WormParams& worm_params,
        const MicrofilariaeParams& mf_params,
        const BlackflyParams& blackfly_params,
        const std::vector<std::unique_ptr<SequelaeParams>>& sequelae_params
    );

    void generate_random_vals_for_diagnostic(std::mt19937& gen);

    // Calculates (or updates) the compliance values for each individual based on the
    // total population and the correlation coefficient, rho.
    double calculate_individual_compliance(std::mt19937& gen);
    void calculate_population_compliance(std::mt19937& gen, double rho, double total_population_coverage, double proportion_never_treated);
    void update_compliance(std::mt19937& gen, double rho, double total_population_coverage);
    void apply_treatment_round(std::mt19937& gen, int minimum_age_of_treatment, int current_time_step);

    // Returns exposure values for every person (age- and sex-weighted,
    // heterogeneity-normalised).
    std::vector<double> get_exposure(
        const ExposureParams& exp_param,
        double sex_ratio
    ) const;

    // Drains the L3 delay buffer for every person and returns newly
    // established worms ready to be sex-split this timestep.
    std::vector<double> get_new_worms();

    // Convenience wrappers for mean blackfly larval loads
    double mean_l1_per_blackfly() const;
    double mean_l2_per_blackfly() const;
    double mean_l3_per_blackfly() const;
    double mean_l3_prevalence_blackflies() const;

    // Returns the mean-L3 value used by wplus1_rate
    double mean_l3_in_blackfly() const;

    double wplus1_rate(
        double mean_l3, double delta_h_zero, double delta_h_inf,
        double c_h, double ABR, double exposure,
        double timestep_years
    ) const;

    // Samples new L3 arrivals from a Poisson distribution and adds them
    // to the L3 delay buffer for each person.
    void new_established_l3(
        std::mt19937& gen,
        double delta_h_zero,
        double delta_h_inf,
        double c_h,
        double ABR,
        const std::vector<double>& exposure,
        double timestep_years
    );

    void update_ov16_status_individual(int indiv_index);

    void update_all_sequelae_individual(
        std::mt19937& generator, int indiv_index,
        double timestep_years, int timesteps_per_year,
        int skin_snip_weight, int num_skin_snips
    );

    void update_all_status(
        std::mt19937& generator,
        double timestep_years, int timesteps_per_year,
        int skin_snip_weight, int num_skin_snips
    );

    bool sample_serostatus_individual(int indiv_index, double sens, double spec);

    void update_worm_dists(double timestep_years);

    void age(
        std::mt19937& gen,
        int current_timestep, double timestep_years,
        const std::vector<double>& new_male_worms,
        const std::vector<double>& new_female_worms,
        std::vector<double>& age_timers
    );

    // Processes deaths: Resets all parasite
    // counts and resamples demographic attributes for dead individuals.
    void process_deaths(std::mt19937& gen);

    std::vector<int> get_sub_population(int age_start, int age_end);
};

#endif
