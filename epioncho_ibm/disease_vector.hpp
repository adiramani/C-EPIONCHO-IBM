#ifndef DISEASE_VECTOR_HPP
#define DISEASE_VECTOR_HPP

#include <vector>
#include "params.hpp"

class VectorPopulation {
public:
    int n_people = 0;

    double v_1 = 0.0;  // L1→L2 per-capita development rate
    double v_2 = 0.0;  // L2→L3 per-capita development rate
    double mu_v = 0.0;  // blackfly per-capita mortality
    double k0 = 0.0; // intercept for l3 intensity/prevalence relationship
    double k1 = 0.0; // slope for l3 intensity/prevalence relationship

    // length = population size
    std::vector<double> l1;
    std::vector<double> l2;
    std::vector<double> l3;

    VectorPopulation() = default;
    
    VectorPopulation(
        int n_people_,
        double v_1_,
        double v_2_,
        double mu_v_,
        double k0_,
        double k1_,
        double init_L1,
        double init_L2,
        double init_L3
    );

    double mean_l1() const;
    double mean_l2() const;
    double mean_l3() const;
    double mean_l3_prevalence() const;
};

class BlackflyPopulation : public VectorPopulation {
public:
    int delay_size = 0;  // l1_delay from BlackflyParams

    // Per-person delay buffer position
    std::vector<int> delay_index;

    // Flat delay buffers: [person * delay_size + slot]
    std::vector<double> delay_l1;
    std::vector<double> delay_mf;
    std::vector<double> delay_exposure;

    // Shared blackfly parameters
    double delta_v0 = 0.0;  // density dependence baseline
    double c_v = 0.0;  // density dependence severity
    double alpha_v = 0.0;  // mf-induced blackfly mortality
    double tou_v = 0.0;  // L1 delay duration

    BlackflyPopulation() = default;

    BlackflyPopulation(
        int n_people_,
        const BlackflyParams& params,
        const std::vector<double>& exposure_heterogeneity,
        const double delta_time_days
    );

    void process_death(int person_idx);

    void update_all(
        double days_in_year,
        double beta,
        const std::vector<double>& exposure_vals,
        const std::vector<double>& mf_loads,
        double a_H,
        double gonotrophic_cycle_length,
        double mu_l3
    );
};

#endif
