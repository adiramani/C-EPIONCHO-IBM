#ifndef SEQUELAE_HPP
#define SEQUELAE_HPP

#include "params.hpp"
#include <random>

class Sequelae {
    protected:
        SequelaeType sequelae_type;
        SequelaeModelType smt;

        int num_individuals;
        double base_probability;
        SequelaeProbTimeUnit prob_timescale;
        int min_age_test;
        int min_infection = 0;
        bool retest_tested_indivs;
        bool use_raw_infection_for_test = false;

        // number of timesteps to countdown for
        int countdown_timesteps;
        // do individuals turn positive or negative at the end of the countdown
        bool status_end_countdown;

        // size = individuals
        std::vector<bool> positives;
        // 0 = not tested or can be retested, 1 = can't be retested, 2 = can't be retested, tested positive
        std::vector<double> has_been_tested;
        std::vector<int> countdowns;        

    public:
        virtual ~Sequelae() = default;
        Sequelae() = default;

        Sequelae(
            SequelaeType sequelae_type, SequelaeModelType sequelae_model_type, int num_individuals,
            double base_probability, SequelaeProbTimeUnit prob_timescale,
            int min_age_test, int min_infection, bool retest_tested_indivs,
            int countdown_timesteps, bool status_end_countdown, bool use_raw_infection_for_test
        );

        SequelaeType get_type();
        SequelaeModelType get_model_type();
        bool use_raw_infection();

        virtual double get_probability(double timestep_years, int timesteps_per_year, int infection_level);

        void decrease_countdown_individual(int indiv_index);

        virtual bool should_test_individual(
            int indiv_index, double age,
            int infection_level,
            bool mating_worm_pair, bool male_female_worm_pair
        );

        virtual void update_all_individuals(
            std::mt19937& generator, std::uniform_real_distribution<double>& uniform_dist,
            double timestep_years, int timesteps_per_year, int num_individuals,
            const std::vector<double>& ages, const std::vector<int>& infection_levels_raw,
            const std::vector<int>& infection_levels_ss, bool mating_worm_pair, bool male_female_worm_pair
        );

        virtual void process_death(int indiv_index);

        int indiv_status(int indiv_index);
};

class TimestepProbSequelae : public Sequelae {
    protected:
        // used for conversion to timestep probability for certain Sequelae
        double average_age;
    public:
        virtual ~TimestepProbSequelae() = default;

        TimestepProbSequelae(
            SequelaeType sequelae_type, SequelaeModelType sequelae_model_type, int num_individuals,
            double base_probability, SequelaeProbTimeUnit prob_timescale,
            int min_age_test, int min_infection, bool retest_tested_indivs,
            int countdown_timesteps, bool status_end_countdown, bool use_raw_infection_for_test,
            double average_age
        );

        double get_probability(double timestep_years, int timesteps_per_year, int infection_level) override;
};

class ExponentialProbSequelae : public Sequelae {
    protected:
        // used for sequelae that have calculated probability based on mf count
        double prob_intercept;
        double prob_slope;
    public:
        virtual ~ExponentialProbSequelae() = default;

        ExponentialProbSequelae(
            SequelaeType sequelae_type, SequelaeModelType sequelae_model_type, int num_individuals,
            double base_probability, SequelaeProbTimeUnit prob_timescale,
            int min_age_test, int min_infection, bool retest_tested_indivs,
            int countdown_timesteps, bool status_end_countdown, bool use_raw_infection_for_test, 
            double prob_intercept, double prob_slope
        );

        double get_probability(double timestep_years, int timesteps_per_year, int infection_level) override;
};

class PowerLawProbSequelae : public Sequelae {
    protected:
        // used for sequelae that have calculated probability based on mf count
        double prob_intercept;
        double prob_slope;
    public:
        virtual ~PowerLawProbSequelae() = default;

        PowerLawProbSequelae(
            SequelaeType sequelae_type, SequelaeModelType sequelae_model_type, int num_individuals,
            double base_probability, SequelaeProbTimeUnit prob_timescale,
            int min_age_test, int min_infection, bool retest_tested_indivs,
            int countdown_timesteps, bool status_end_countdown, bool use_raw_infection_for_test, 
            double prob_intercept, double prob_slope
        );

        double get_probability(double timestep_years, int timesteps_per_year, int infection_level) override;
};

class OAESequelae : public PowerLawProbSequelae {
    protected:
        std::mt19937& gen;
        std::uniform_real_distribution<double> random_age_to_test;
        std::vector<int> ages_to_test;
    public:
        virtual ~OAESequelae() = default;

        OAESequelae(
            SequelaeType sequelae_type, SequelaeModelType sequelae_model_type, int num_individuals,
            double base_probability, SequelaeProbTimeUnit prob_timescale,
            int min_age_test, int min_infection, bool retest_tested_indivs,
            int countdown_timesteps, bool status_end_countdown, bool use_raw_infection_for_test, 
            double prob_intercept, double prob_slope,
            std::mt19937& gen, int max_age_test,
            std::vector<double> initial_ages
        );

        bool should_test_individual(
            int indiv_index, double age,
            int infection_level,
            bool mating_worm_pair, bool male_female_worm_pair
        ) override;

        void process_death(int indiv_index) override;
};

#endif
