#include "disease_vector.hpp"
#include <cmath>
#include <numeric>
#include <algorithm>

VectorPopulation::VectorPopulation(
    int n_people_,
    double v_1_,
    double v_2_,
    double mu_v_,
    double k0_,
    double k1_,
    double init_L1,
    double init_L2,
    double init_L3
) :     
    n_people(n_people_),
    v_1(v_1_),
    v_2(v_2_),
    mu_v(mu_v_),
    k0(k0_),
    k1(k1_),
    l1(n_people_, init_L1),
    l2(n_people_, init_L2),
    l3(n_people_, init_L3)
{}


BlackflyPopulation::BlackflyPopulation(
    int n_people_,
    const BlackflyParams& params,
    const std::vector<double>& exposure_heterogeneity
) : 
    VectorPopulation(
        n_people_,
        params.l1_l2_per_larva_per_year,
        params.l2_l3_per_larva_per_year,
        params.blackfly_mort_per_fly_per_year,
        params.k0,
        params.k1,
        params.initial_L1,
        params.initial_L2,
        params.initial_L3
    ),
    delay_size((int)params.l1_delay),
    delay_index(n_people_, 0),
    delay_l1(n_people_ * (int)params.l1_delay, params.initial_L1),
    delay_mf(n_people_ * (int)params.l1_delay, 0.0),
    delay_exposure(n_people_ * (int)params.l1_delay, 0.0),
    delta_v0(params.delta_v0),
    c_v(params.c_v),
    alpha_v(params.blackfly_mort_from_mf_per_fly_per_year),
    tou_v(params.l1_delay)
{
    // Set the exposure heterogeneity to the initial value
    for (int p = 0; p < n_people_; ++p) {
        for (int s = 0; s < delay_size; ++s)
            delay_exposure[p * delay_size + s] = exposure_heterogeneity[p];
    }
}

void BlackflyPopulation::process_death(int person_idx) {
    l1[person_idx] = 0.0;
    // l2 and l3 are not reset on death in the original model
    const int base = person_idx * delay_size;
    for (int s = 0; s < delay_size; ++s) {
        delay_l1[base + s] = 0.0;
        delay_mf[base + s] = 0.0;
        // exposure delay is preserved, as it is just individual heterogeneity
    }
}

void BlackflyPopulation::update_all(
    double timestep_years,
    double beta,
    const std::vector<double>& exposure_vals,
    const std::vector<double>& mf_loads,
    double a_H,
    double gonotrophic_cycle_length,
    double mu_l3
) {
    for (int p = 0; p < n_people; ++p) {
        const int person_loc_flat = p * delay_size;
        const int delay_loc = delay_index[p];
        const double mf = mf_loads[p];
        const double exp = exposure_vals[p];

        // ---------- calc_L1 ----------
        const double pim = mu_v + alpha_v * mf * exp;
        const double density_dep = delta_v0 / (1.0 + c_v * mf * exp);
        const double delay_pim = mu_v + alpha_v * delay_mf[person_loc_flat + delay_loc] * delay_exposure[person_loc_flat + delay_loc];
        const double new_l1 = (density_dep * beta * exp * mf)
                            / (pim + v_1 * std::exp(-(tou_v * timestep_years) * delay_pim));

        // ---------- calc_L2 (uses delay slot values before overwrite) ----------
        const double new_l2 = (
            (delay_l1[person_loc_flat + delay_loc]
            * v_1 * std::exp(-(tou_v * timestep_years) * delay_pim))
        ) / (v_2 + mu_v);

        // ---------- calc_L3 (uses old l2 before this timestep's update) ----------
        const double new_l3 = (v_2 * l2[p])
                            / (mu_v + mu_l3 + (a_H / gonotrophic_cycle_length));


        l1[p] = new_l1;
        l2[p] = new_l2;
        l3[p] = new_l3;

        for (int d = 0; d < delay_size; ++d)
            delay_l1[person_loc_flat + d] = new_l1;
        delay_mf[person_loc_flat + delay_loc] = mf;
        delay_exposure[person_loc_flat + delay_loc] = exp;

        delay_index[p] = (delay_loc + 1 >= delay_size) ? 0 : delay_loc + 1;
    }
}

double VectorPopulation::mean_l1() const {
    return std::accumulate(l1.begin(), l1.end(), 0.0) / n_people;
}
double VectorPopulation::mean_l2() const {
    return std::accumulate(l2.begin(), l2.end(), 0.0) / n_people;
}
double VectorPopulation::mean_l3() const {
    return std::accumulate(l3.begin(), l3.end(), 0.0) / n_people;
}
double VectorPopulation::mean_l3_prevalence() const {
    double mean_l3_per_blackfly = mean_l3();
    double k = k0 + k1 * mean_l3_per_blackfly;
    return 1 - pow((1 + mean_l3_per_blackfly / k), -k);
}
