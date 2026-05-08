#ifndef SEQUELAE_HPP
#define SEQUELAE_HPP

#include "params.hpp"
#include <random>

class Sequelae {
    private:
        // size = num_sequelae
        std::vector<SequelaeType> sequelae_types;
        std::vector<double> probabilities;
        std::vector<SequelaeProbTimeUnit> timescales;
        std::vector<int> countdown_timesteps;
        // do individuals turn positive or negative at the end of the countdown
        std::vector<bool> end_countdown;
        // used for conversion to timestep probabilitiy for certain Sequelae
        std::vector<double> average_ages;
        
        // size = individuals * num_sequelae
        std::vector<bool> positives;
        std::vector<bool> has_been_tested;
        std::vector<int> countdowns;

        int num_sequelae;
        int num_individuals;
        

    public:
        Sequelae() = default;

        Sequelae(
            std::vector<SequelaeType> sequelae_types, int num_individuals,
            std::vector<double> probabilities, std::vector<SequelaeProbTimeUnit> timescales,
            std::vector<int> countdown_timesteps, std::vector<double> average_ages
        );

        double calculate_timestep_probability(int sequelae_index, double timestep_years, int days_in_year);

        void decrease_countdown_individual(int indiv_index);

        bool should_test_individual(
            SequelaeType st, int person_sequelae_index, double age, 
            double raw_mf_count, double observed_mf_count
        );

        void update_sequelae(SequelaeType st, int person_sequelae_index, bool positive);

        void update_all_sequelae_indiv(
            std::mt19937& generator, std::uniform_real_distribution<double>& uniform_dist,
            int indiv_index, double timestep_years, int days_in_year,
            double age, double raw_mf_count, double observed_mf_count
        );

        void process_death(int indiv_index);

        int indiv_status(SequelaeType sequelae_type, int indiv_index);
};

#endif
