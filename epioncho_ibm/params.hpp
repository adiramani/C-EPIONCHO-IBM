#ifndef PARAMS_HPP
#define PARAMS_HPP

#include <string>
#include <vector>
#include <optional>
#include <stdexcept>
#include <format>
#include "enums.hpp"

// -------------------- Treatment Parameters --------------------

struct InterventionParams {
    // Inclusive of start time, but not end time
    int start_time = 0;
    int end_time = 1;
    int interval_years = 1;
    bool pre_converted_timesteps = false;
    std::string intervention_name = "";
    InterventionType intervention_type;
    std::vector<int> application_times;

    InterventionParams(
        int start, int end, int interval_years,
        const std::string& name,
        InterventionType type
    ) 
    : start_time(start), end_time(end), interval_years(interval_years),
      intervention_name(name), intervention_type(type)
    {
        calculate_application_times();
    }

    InterventionParams(
        std::vector<int> application_times,
        bool pre_converted_timesteps,
        const std::string& name,
        InterventionType type
    ) 
    : pre_converted_timesteps(pre_converted_timesteps),
      intervention_name(name), intervention_type(type),
      application_times(application_times)
    {}

    void calculate_application_times() {
        application_times.clear();
        if (interval_years <= 0 || end_time < start_time) return;

        for (int t = start_time; t < end_time; t += interval_years) {
            application_times.push_back(t);
        }
    }
};

struct DrugParams {
    double microfilaricidal_upsilon = 0.0096;
    double microfilaricidal_kappa = 1.25;
    double embryostatic_lambda_max = 32.4;
    double embryostatic_phi = 19.6;
    double permanent_infertility = 0.345;
};

struct TreatmentParams : public InterventionParams {
    DrugParams drug_params;
    int min_age_of_treatment = 5;
    double rho = 0.3;
    double total_population_coverage = 0.65;
    double proportion_never_treated = 0.0;

    TreatmentParams(
        int start, int end, int interval_years,
        std::string intervention_name = "IVM",
        DrugParams drug_params = {},
        int min_age = 5,
        double rho = 0.3,
        double proportion_never_treated = 0.0,
        double coverage = 0.65
    )
    : InterventionParams(start, end, interval_years, intervention_name, InterventionType::MDA),
      drug_params(drug_params),
      min_age_of_treatment(min_age),
      rho(rho),
      total_population_coverage(coverage),
      proportion_never_treated(proportion_never_treated)
    {}

    TreatmentParams(
        std::vector<int> application_times,
        bool pre_converted_timesteps,
        std::string intervention_name = "IVM",
        DrugParams drug_params = {},
        int min_age = 5,
        double rho = 0.3,
        double proportion_never_treated = 0.0,
        double coverage = 0.65
    )
    : InterventionParams(application_times, pre_converted_timesteps, intervention_name, InterventionType::MDA),
      drug_params(drug_params),
      min_age_of_treatment(min_age),
      rho(rho),
      total_population_coverage(coverage),
      proportion_never_treated(proportion_never_treated)
    {}
};

struct VectorControlParams: public InterventionParams {
    std::vector<double> efficacies;

    VectorControlParams(
        int start, int end, int interval_years,
        double efficacy, int bounce_back_intervals = 0,
        std::string intervention_name = ""
    )
    : InterventionParams(
        start, end + 1 + (bounce_back_intervals * interval_years), interval_years,
        intervention_name, InterventionType::VectorControl
      )
    {
        for (std::size_t t = 0; t < application_times.size(); ++t)
            efficacies.push_back(efficacy);
    }

    VectorControlParams(
        int start, int end, int interval_years,
        std::vector<double> efficacy, int bounce_back_intervals = 0,
        std::string intervention_name = ""
    )
    : InterventionParams(
        start, end + (bounce_back_intervals * interval_years), interval_years,
        intervention_name, InterventionType::VectorControl
      ),
      efficacies(efficacy)
    {
        if (efficacies.size() != application_times.size()) {
            throw std::invalid_argument(
                "Efficacy size must match application_times size."
            );
        }
    }

    VectorControlParams(
        std::vector<int> application_times,
        bool pre_converted_timesteps,
        std::vector<double> efficacy,
        std::string intervention_name = ""
    )
    : InterventionParams(
        application_times,
        pre_converted_timesteps,
        intervention_name, InterventionType::VectorControl
      ),
      efficacies(efficacy)
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
    double l3_delay = 10.0 * 28;  // days, assuming 28 days per month
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
    bool use_kmf_const = true;
    double slope_kmf = 0.0478;
    double initial_kmf = 0.313;
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
    Params params;
    std::vector<TreatmentParams> treatments;
    std::vector<VectorControlParams> vector_control;

    InputParams(
        Params params, 
        std::vector<TreatmentParams> treatments = {},
        std::vector<VectorControlParams> vector_control = {}
    )
    : params(params),
      treatments(treatments),
      vector_control(vector_control)
    {}
};

// --------------- Output Struct ----------------------------
struct OutputInfo {
    std::vector<ModelOutputTypes> outputs_to_track;
    std::vector<double> output_time_years;
    int start_age = 0;
    int end_age = 80;

    OutputInfo(
        double end_time_years, double interval_years = 1,
        int start_age = 0, int end_age = 80,
        std::vector<ModelOutputTypes> outputs_to_track = {ModelOutputTypes::mf_prevalence, ModelOutputTypes::ov16_seroprevalence}
    )
    : outputs_to_track(outputs_to_track),
      start_age(start_age), end_age(end_age)
    {
        for (int t = 0; t < end_time_years; t += interval_years) {
            output_time_years.push_back(t);
        }
    }

    OutputInfo(
        std::vector<double> output_time_years,
        int start_age = 0, int end_age = 80,
        std::vector<ModelOutputTypes> outputs_to_track = {ModelOutputTypes::mf_prevalence, ModelOutputTypes::ov16_seroprevalence}
    )
    : outputs_to_track(outputs_to_track),
      output_time_years(output_time_years),
      start_age(start_age), end_age(end_age)
    {}
};

#endif
