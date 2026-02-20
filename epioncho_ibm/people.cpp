#include "people.hpp"
#include <random>
#include <algorithm>
#include <numeric>
#include <assert.h>
#include <iostream>

Worm create_new_worm(WormParams worm_params) {
  return(
    Worm(
      worm_params.worm_age_stages,
      worm_params.q_w,
      worm_params.y_w,
      worm_params.d_w,
      worm_params.initial_worms,
      worm_params.epsilon,
      worm_params.omega,
      worm_params.lambda_zero,
      worm_params.fecundity_worms_F,
      worm_params.fecundity_worms_G
    )
  );
}

People::People(int population_size_)
  : population_size(population_size_)
{}

void People::initialize_from_params(
  std::mt19937& gen,
  double timestep_years,
  int mean_human_age_,
  int max_human_age_,
  double k_e,
  WormParams worm_params,
  MicrofilariaeParams mf_params,
  BlackflyParams blackfly_params
)
{
  mean_age = mean_human_age_;
  
  ages.resize(population_size);
  std::exponential_distribution<int> age_exp_dist(1 / mean_human_age_);
  gender_dist = std::bernoulli_distribution(0.5);
  gamma_dist = std::gamma_distribution<double>(k_e, 1 / k_e);
  death_dist = std::bernoulli_distribution((1/mean_age) * timestep_years);


  for (auto worm_type : worm_types) {
    worms[worm_type].resize(population_size);
  }

  microfilariae.resize(population_size);
  l3.resize(population_size);
  blackflies.resize(population_size);
  sex.resize(population_size);
  exposure_heterogeneity.resize(population_size);

  for (int i = 0; i < population_size; ++i) {
    ages[i] = std::min(age_exp_dist(gen), max_human_age_);
    sex[i] = gender_dist(gen);
    exposure_heterogeneity[i] = gamma_dist(gen);

    for (auto worm_type : worm_types) {
      worms[worm_type][i] = create_new_worm(worm_params);
    }
    
    microfilariae[i] = MF(
      mf_params.mf_age_stages,
      mf_params.q_m,
      mf_params.y_m,
      mf_params.d_m,
      mf_params.mf_move_rate,
      0
    );
    
    blackflies[i] = Blackfly(
      blackfly_params.delta_v0,
      blackfly_params.c_v,
      blackfly_params.l1_l2_per_larva_per_year,
      blackfly_params.l2_l3_per_larva_per_year,
      blackfly_params.blackfly_mort_per_fly_per_year,
      blackfly_params.blackfly_mort_from_mf_per_fly_per_year,
      blackfly_params.l1_delay,
      blackfly_params.initial_L1,
      blackfly_params.initial_L2,
      blackfly_params.initial_L3,
      exposure_heterogeneity[i]
    );
    
    l3[i] = L3(
      worm_params.l3_delay,
      worm_params.l3_delay / 365
    );
  }
}

std::vector<double> People::get_exposure(ExposureParams exp_param, double sex_ratio) {
  double exp_females = 1 / (sex_ratio * (exp_param.Q - 1) + 1);
  double exp_males = (exp_param.Q) * exp_females;
  int n_females = std::count(sex.begin(), sex.end(), true);
  int n_males = population_size - n_females;

  assert(n_males > 0 && n_females > 0);

  std::vector<double> age_sex_exp(population_size);
  double tmp_alpha_m = 0.0, tmp_alpha_f = 0.0;
  
  for (int i = 0; i < population_size; ++i) {
    if (sex[i] == 0) {
      // Males
      age_sex_exp[i] = std::exp(-exp_param.male_exposure_exponent * ages[i]);
      tmp_alpha_m += age_sex_exp[i];
    } else {
      // Females
      age_sex_exp[i] = std::exp(-exp_param.female_exposure_exponent * ages[i]);
      tmp_alpha_f += age_sex_exp[i];
    }
  }

  double male_scale = exp_males / (tmp_alpha_m / n_males);
  double female_scale = exp_females / (tmp_alpha_f / n_females);

  for (int i = 0; i < population_size; ++i) {
    age_sex_exp[i] *= (sex[i] == 0) ? male_scale : female_scale;
  }

  double mean_exp_het = std::accumulate(exposure_heterogeneity.begin(), exposure_heterogeneity.end(), 0.0) / population_size;
  for (int i = 0; i < population_size; ++i) {
    age_sex_exp[i] *= (exposure_heterogeneity[i] / mean_exp_het);
  }

  double mean_age_sex_exp = std::accumulate(age_sex_exp.begin(), age_sex_exp.end(), 0.0) / population_size;
  for (double& val : age_sex_exp) {
    val /= mean_age_sex_exp;
  }

  return age_sex_exp;
}

double People::wplus1_rate(
  double mean_l3, double delta_h_zero, double delta_h_inf, double c_h,
  double ABR, double exposure, double timestep_in_years
)
{
  double delta_h = (
    delta_h_zero + delta_h_inf * c_h * ABR * mean_l3 * exposure
  ) /
  (
    1 + c_h * ABR * mean_l3 * exposure
  );
  return timestep_in_years * ABR * delta_h * exposure * mean_l3;
}

void People::new_established_l3(
  std::mt19937& gen, double delta_h_zero, double delta_h_inf, double c_h,
  double ABR, std::vector<double>& exposure, double timestep_in_years
) {
  double mean_l3 = mean_l3_per_blackfly();
  for (int i = 0; i < population_size; ++i) {
    double l3_in_rate = wplus1_rate(mean_l3, delta_h_zero, delta_h_inf, c_h, ABR, exposure[i], timestep_in_years);
    int new_l3_in = std::poisson_distribution(l3_in_rate)(gen);
    // if (new_l3_in > 0) {
    //   std::cout << l3[i].current_index << " New Worms " << new_l3_in << "\n";
    // }
    l3[i].add_new_established_l3(new_l3_in);
  }
}

void People::age(std::mt19937 gen, double timestep_years, std::vector<double>& new_male_worms, std::vector<double>& new_female_worms) {
  for (int i = 0; i < population_size; ++i) {
    std::vector<double> fertile_female_swapped_worms;
    std::vector<double> infertile_female_swapped_worms;
    
    for (auto worm_type : worm_types) {
      double new_worms = 0.0;
      if (worm_type == "infertile_female") {
        new_worms = new_female_worms[i];
      } else if (worm_type == "male") {
        new_worms = new_male_worms[i];
      }

      double aging_rate = timestep_years / worms[worm_type][i].years_per_compartment;
      std::vector<double> swapped_worms = worms[worm_type][i].age(
        std::bernoulli_distribution(aging_rate),
        gen,
        worm_type,
        timestep_years,
        new_worms
      );

      if (worm_type == "infertile_female") {
        fertile_female_swapped_worms = swapped_worms;
      } else if (worm_type == "fertile_female") {
        infertile_female_swapped_worms = swapped_worms;
      }
    }
    worms["fertile_female"][i].age_helper_female_swapped_worms(fertile_female_swapped_worms);
    worms["infertile_female"][i].age_helper_female_swapped_worms(infertile_female_swapped_worms);
    microfilariae[i].age(timestep_years, worms["male"][i].get_worm_load(), worms["fertile_female"][i]);
    l3[i].age();
    ages[i] += 1;
  }
}

void People::process_deaths(std::mt19937 gen) {
  std::vector<int> people_to_die(population_size);

  int num_dead = 0;
  for (int i = 0; i < population_size; ++i) {
    bool is_dead = death_dist(gen);
    people_to_die[i] = is_dead;
    num_dead += is_dead;
    if (is_dead) {
      ages[i] = 0;
      sex[i] = gender_dist(gen);
      exposure_heterogeneity[i] = gamma_dist(gen);

      for (auto worm_type : worm_types) {
        worms[worm_type][i].process_death();
      }

      microfilariae[i].process_death();
      l3[i].process_death();
      blackflies[i].process_death();
    }
  }
  // if (num_dead > 0) {
  //   std::cout << "Dead: " << num_dead << "\n";
  // }
}

std::vector<double> People::get_new_worms() {
  std::vector<double> new_worms(population_size, 0.0);
  for(int i = 0; i < population_size; ++i) {
    new_worms[i] = l3[i].get_new_worms();
  }
  return new_worms;
}

double People::mean_l3_per_blackfly() {
  double total_l3 = 0.0;
  for(int i = 0; i < population_size; ++i) {
    total_l3 += blackflies[i].l3;
  }
  return total_l3 / population_size;
}