#ifndef WORM_HPP
#define WORM_HPP

#include <vector>
#include <string>
#include <random>
#include "enums.hpp"

class Parasite {
    public:
        int compartments; // c_max
        float years_per_compartment; // q_W | q_M

        float y_mortality_rate_by_age; // y_W | w_M
        float d_mortality_rate_by_age; // d_W | d_M

        std::vector<double> parasites;
        std::vector<double> age_categories;

        Parasite() = default;

        Parasite(
            int compartments_,
            float years_per_compartment_,
            float y_mortality_rate_by_age_,
            float d_mortality_rate_by_age_,
            int initial_parasites_
        );

        double get_parasite_load();

        void process_death();

        double weibull_mortality(int compartment, double timestep_years);
};

class Worm: public Parasite {
    public:
        float mf_production_rate; // ε*
        float infertile_to_fertile_rate; // ω
        float fertile_to_infertile_rate; // λ0
        float F;
        float G;
    
        Worm() = default;

        Worm(
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
        );

        double get_worm_load() {
            return get_parasite_load();
        }

        double fecundity_rate(int compartment);

        int fecundity_movement(int num_worms, WormType worm_type, double timestep_years, std::mt19937& gen);

        void age(
            std::mt19937& gen, WormType worm_type, double timestep_years, double new_worms,
            std::vector<double>* swapped_out
        );

        void age_helper_female_swapped_worms(
            std::vector<double>& incoming_swapped_worms
        );
};

class MF: public Parasite {
    public:
        MF() = default;
        double mf_move_rate;

        MF(
            int compartments_,
            float years_per_compartment_,
            float y_mortality_rate_by_age_,
            float d_mortality_rate_by_age_,
            double mf_move_rate_,
            int initial_parasites_
        );

        void calc_new_mf(double timestep_in_years, double male_worms, Worm& fertile_female_worms);

        void age_exiting_mf(int age_category, double timestep_years, double prev_mf);

        void age(double timestep_years, double male_worms, Worm& fertile_female_worms);

        double get_mf_load() {
            return get_parasite_load();
        };
};

class L3: public Parasite {
    public:
        int current_index;

        L3() = default;

        L3(
            int compartments,
            float years_per_compartment
        );

        void age();

        double get_new_worms();

        void add_new_established_l3(int num_l3);
};

#endif