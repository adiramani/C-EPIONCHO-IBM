#ifndef PARAMS_HPP
#define PARAMS_HPP

#include <string>
#include <vector>
#include <optional>

// -------------------- Treatment Parameters --------------------

struct InterventionParams {
    int start_time = 0;
    int end_time = 0;
    int interval_years = 1;
    std::string intervention_name = "";
    std::string intervention_type = "";
    std::vector<int> application_timesteps;

    InterventionParams(int start, int end, int interval,
                       const std::string& name = "",
                       const std::string& type = "")
        : start_time(start), end_time(end), interval_years(interval),
          intervention_name(name), intervention_type(type)
    {
        calculate_application_timesteps();
    }

    void calculate_application_timesteps() {
        application_timesteps.clear();
        if (interval_years <= 0 || end_time < start_time) return;

        for (int t = start_time; t <= end_time; t += interval_years) {
            application_timesteps.push_back(t);
        }
    }
};

struct TreatmentParams : public InterventionParams {
    double microfilaricidal_nu = 0.0096;
    double microfilaricidal_omega = 1.25;
    double embryostatic_lambda_max = 32.4;
    double embryostatic_phi = 19.6;
    double permanent_infertility = 0.345;
    int min_age_of_treatment = 5;
    double correlation = 0.3;
    double total_population_coverage = 0.65;

    TreatmentParams();

        TreatmentParams(
        int start, int end, int interval,
        std::string intervention_name,
        std::string intervention_type,
        double nu = 0.0096,
        double omega = 1.25,
        double lambda_max = 32.4,
        double phi = 19.6,
        double perm_infertility = 0.345,
        int min_age = 5,
        double corr = 0.3,
        double coverage = 0.65
    )
    : InterventionParams(start, end, interval, "", ""),
      microfilaricidal_nu(nu),
      microfilaricidal_omega(omega),
      embryostatic_lambda_max(lambda_max),
      embryostatic_phi(phi),
      permanent_infertility(perm_infertility),
      min_age_of_treatment(min_age),
      correlation(corr),
      total_population_coverage(coverage)
    {}
};

struct VectorControlParams: public InterventionParams {
    double efficacy = 0;
    int bounce_back_timesteps = 1;

    VectorControlParams();
    VectorControlParams(
        int start, int end, int interval,
        std::string intervention_name,
        std::string intervention_type,
        double efficacy, int bounce_back_timesteps
    )
    : InterventionParams(start, end, interval, "", ""),
      efficacy(efficacy),
      bounce_back_timesteps(bounce_back_timesteps)
    {}
};

// -------------------- Worm Parameters --------------------
struct WormParams {
    double y_w = 0.09953;
    double d_w = 6.00569;
    int initial_worms = 1;
    int worm_age_stages = 21;
    double q_w = 1.0; // year(s) spend in each worm compartmenmt
    int max_worm_age = 21;
    double fecundity_worms_F = 70.0;
    double fecundity_worms_G = 0.72;
    double omega = 0.59;
    double lambda_zero = 0.33;
    double sex_ratio = 0.5;
    double epsilon = 1.158305; // mf production per worm
    double l3_delay = 10.0 * 30;  // months, assuming 30 days per month
};

// -------------------- Blackfly Parameters --------------------
struct BlackflyParams {
    double delta_h_zero = 0.1864987;
    double delta_h_inf = 0.002772749;
    double blackfly_mort_per_fly_per_year = 26.0;
    double blackfly_mort_from_mf_per_fly_per_year = 0.39;
    double mu_L3 = 52.0;
    double a_H = 0.8;
    double l1_l2_per_larva_per_year = 201.6189;
    double l2_l3_per_larva_per_year = 207.7384;
    double delta_v0 = 0.0166;
    double c_v = 0.0205;
    double initial_L3 = 0.03;
    double initial_L2 = 0.03;
    double initial_L1 = 0.03;
    double human_blood_index = 0.63;
    double gonotrophic_cycle_length = 1.0 / 104.0;
    double c_h = 0.004900419;
    double bite_rate_per_person_per_year = 1000.0;
    double l1_delay = 4.0;   // days
};

// -------------------- Microfilaria Parameters --------------------
struct MicrofilariaeParams {
    double mf_move_rate = 8.13333;
    int mf_age_stages = 21;
    double q_m = 0.125; // year(s) spend in each mf compartmenmt
    double max_mf_age = 2.5;
    double initial_mf = 0.0;
    double y_m = 1.089;
    double d_m = 1.428;
    double kmf_const = 15;
    // double slope_kmf = 0.0478;
    // double initial_kmf = 0.313;
};

// -------------------- Exposure Parameters --------------------
struct ExposureParams {
    double Q = 1.2;
    double male_exposure_exponent = 0.007;
    double female_exposure_exponent = -0.023;
};

// -------------------- Human Parameters --------------------
struct HumanParams {
    int min_skinsnip_age = 5;
    int max_human_age = 80;
    int mean_human_age = 50;
    int skin_snip_weight = 2;
    int skin_snip_number = 2;
    double gender_ratio = 0.5;
};

// -------------------- Base Parameters --------------------
struct BaseParams {
    int seed = 0;
    int n_people = 1000;
    double k_E = 0.3;
    double n_treatments_bin_size = 1.0;
    double delta_time_days = 1.0;
    double year_length_days = 365.0;
    double month_length_days = 28.0;
    std::vector<int> sequela_active;  // placeholder for SequelaType
};

// -------------------- Top-level Params Struct --------------------
struct Params {
    BaseParams base;
    WormParams worms;
    BlackflyParams blackfly;
    MicrofilariaeParams mf;
    ExposureParams exposure;
    HumanParams human;

    Params() = default;
};

struct InputParams {
    // std::vector<TreatmentParams> treatments;
    // std::vector<VectorControlParams> vector_control;
    Params params;

    InputParams() = default;
};

#endif
