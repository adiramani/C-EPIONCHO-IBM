#include "vector.hpp"
#include <vector>
#include <numeric>
#include <iostream>

Vector::Vector(
  double delta_v0_,
  double c_v_,
  double v_1_,
  double v_2_,
  double mu_v_,
  double alpha_v_,
  double tou_v_,
  double init_l1,
  double init_l2,
  double init_l3
)
  : delta_v0(delta_v0_),
    c_v(c_v_),
    v_1(v_1_),
    v_2(v_2_),
    mu_v(mu_v_),
    alpha_v(alpha_v_),
    tou_v(tou_v_),
    l1(init_l1),
    l2(init_l2),
    l3(init_l3)
{}

Blackfly::Blackfly(
  double delta_v0_,
  double c_v_,
  double v_1_,
  double v_2_,
  double mu_v_,
  double alpha_v_,
  double tou_v_,
  double init_l1,
  double init_l2,
  double init_l3,
  double init_exposure_heterogeneity
)
  : Vector(
      delta_v0_,
      c_v_,
      v_1_,
      v_2_,
      mu_v_,
      alpha_v_,
      tou_v_,
      init_l1,
      init_l2,
      init_l3
    ),
    delays(tou_v_, {init_l1, 0.0, init_exposure_heterogeneity}),
    delay_index(0)
{}

// Only l1 goes to 0 upon death
void Vector::process_death() {
  l1 = 0;
}

void Blackfly::update_delay_index(double l1, double mf, double exposure) {
  delays[delay_index] = {l1, mf, exposure};
  if (delay_index + 1 >= delays.size()) {
    delay_index = 0;
  } else {
    delay_index += 1;
  }
}

double Blackfly::calc_density_dependence_i(double mf, double exposure) {
  return(
    delta_v0 /
    (
      1 + (c_v * mf * exposure)
    )
  );
}

void Blackfly::calc_L1(double timestep_years, double beta, double exposure, float mf) {
  // for (int i = 0; i < l1.size(); ++i) {
  double parasite_induced_mortality = mu_v + alpha_v * mf * exposure;
  l1 = (
    beta * 
    calc_density_dependence_i(
      mf, exposure
    ) * 
    exposure * mf
  ) /
  (
    (parasite_induced_mortality) +
    v_1*std::exp(-(tou_v * timestep_years) * (mu_v + alpha_v * delays[delay_index].microfilariae * delays[delay_index].exposure))
  );
  // }
}

void Blackfly::calc_L2(double timestep_years) {
  double parasite_induced_mortality = mu_v + alpha_v * delays[delay_index].microfilariae * delays[delay_index].exposure;
  l2 = (delays[delay_index].l1 * (v_1 * std::exp(-(tou_v * timestep_years) * (parasite_induced_mortality)))) / (v_2 + mu_v);
}

void Blackfly::calc_L3(double a_H, double g, double mu_l3) {
  l3 = (v_2 * l2) / (mu_v + mu_l3 + (a_H / g));
}

void Blackfly::process_death() {
  Vector::process_death();
  for (auto& triplet : delays) {
    triplet = {0, 0, 0};
  }
}
