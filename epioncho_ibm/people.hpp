#ifndef PEOPLE_HPP
#define PEOPLE_HPP

#include "parasite.hpp"
#include "vector.hpp"
#include "params.hpp"
#include <cstddef>
#include <map>
#include <vector>
#include <string>
#include <random>

class People {
    private:
        const std::vector<std::string> worm_types = {"fertile_female", "infertile_female", "sterilized_female", "male"};
    public:
        int population_size;
        double mean_age;
        std::gamma_distribution<double> gamma_dist;
        std::bernoulli_distribution gender_dist;
        std::bernoulli_distribution death_dist;

        std::vector<int> ages;
        std::vector<bool> sex; // 0 = male, 1 = female
        std::vector<double> exposure_heterogeneity;
        std::vector<bool> ov16_serostatus;
        std::vector<bool> ov16_serostatus_finite_seroreversion;

        std::vector<L3> l3;
        std::map<std::string, std::vector<Worm>> worms;
        std::vector<MF> microfilariae;
        std::vector<Blackfly> blackflies;
        
        People(int population_size_);

        void initialize_from_params(
            std::mt19937& gen,
            double timestep_years,
            int mean_human_age_,
            int max_human_age_,
            double k_e,
            WormParams worm_params,
            MicrofilariaeParams mf_params,
            BlackflyParams blackfly_params
        );

        std::vector<double> get_exposure(ExposureParams exp_param, double sex_ratio);

        std::vector<double> get_new_worms();

        double mean_l3_per_blackfly();

        double wplus1_rate(double mean_l3, double delta_h_zero, double delta_h_inf, double c_h, double ABR, double exposure, double timestep_in_years);

        void new_established_l3(std::mt19937& gen, double delta_h_zero, double delta_h_inf, double c_h, double ABR, std::vector<double>& exposure, double timestep_in_years);

        void age(std::mt19937 gen, double timestep_years, std::vector<double>& new_male_worms, std::vector<double>& new_female_worms);

        void process_deaths(std::mt19937 gen);
};

#endif
