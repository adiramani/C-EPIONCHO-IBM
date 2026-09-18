#ifndef CONFIG_PARSER_HPP
#define CONFIG_PARSER_HPP

#include <nlohmann/json.hpp>
#include "params.hpp"
#include "oncho_params.hpp"
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <iostream>

using json = nlohmann::json;

class DrugRegistry {
private:
    static std::map<std::string, DrugParams> existing_drugs;
    static bool initialized;
    
    static void initialize() {
        if (initialized) return;
        existing_drugs["IVM"] = DrugParamsIVM();
        existing_drugs["MOX"] = DrugParamsMOX();
        initialized = true;
    }

public:
    static bool drug_exists(const std::string& drug_name) {
        initialize();
        return existing_drugs.find(drug_name) != existing_drugs.end();
    }

    static DrugParams get_drug(const std::string& drug_name) {
        initialize();
        auto it = existing_drugs.find(drug_name);
        if (it != existing_drugs.end()) {
            return it->second;
        }
        throw std::invalid_argument("Unknown drug: " + drug_name + ". Known drugs: IVM, MOX");
    }

    static std::vector<std::string> list_drugs() {
        initialize();
        std::vector<std::string> names;
        for (const auto& pair : existing_drugs) {
            names.push_back(pair.first);
        }
        return names;
    }
};

std::map<std::string, DrugParams> DrugRegistry::existing_drugs;
bool DrugRegistry::initialized = false;

class Validator {
public:
    template <typename T>
    static void validate_in_array(T value, std::vector<T> valid_values, const std::string& param_name) {
        if (std::find(valid_values.begin(), valid_values.end(), value) == valid_values.end()) {
            throw std::invalid_argument(
                param_name + " input value of " + std::to_string(value) +
                "is not valid."
            );
        }
    }

    static void validate_range(double value, double min_val, double max_val, const std::string& param_name, bool inclusive = true) {
        bool invalid = inclusive ? (value < min_val || value > max_val) : (value <= min_val || value >= max_val);
        if (invalid) {
            throw std::invalid_argument(
                param_name + " must be in range [" + std::to_string(min_val) + ", " + 
                std::to_string(max_val) + "], got " + std::to_string(value)
            );
        }
    }

    static void validate_positive(double value, const std::string& param_name) {
        if (value <= 0) {
            throw std::invalid_argument(param_name + " must be positive, got " + std::to_string(value));
        }
    }

    static void validate_non_negative(double value, const std::string& param_name) {
        if (value < 0) {
            throw std::invalid_argument(param_name + " must be non-negative, got " + std::to_string(value));
        }
    }

    static void validate_positive_int(int value, const std::string& param_name) {
        if (value <= 0) {
            throw std::invalid_argument(param_name + " must be positive, got " + std::to_string(value));
        }
    }

    static void validate_non_negative_int(int value, const std::string& param_name) {
        if (value < 0) {
            throw std::invalid_argument(param_name + " must be non-negative, got " + std::to_string(value));
        }
    }

    static void validate_probability(double value, const std::string& param_name = "probability") {
        validate_range(value, 0.0, 1.0, param_name);
    }
};

DrugParams parse_drug_from_json(const json& drug_json) {
    DrugParams drug;

    if (drug_json.contains("name")) {
        std::string name = drug_json["name"].get<std::string>();
        if (DrugRegistry::drug_exists(name)) {
            drug = DrugRegistry::get_drug(name);
        } else {
            throw std::invalid_argument("Drug name not recognized: " + name);
        }
    } else {
        throw std::invalid_argument("Drug must have a 'name' field");
    }

    // Override builtin parameters if custom values provided
    if (drug_json.contains("microfilaricidal_upsilon")) {
        drug.microfilaricidal_upsilon = drug_json["microfilaricidal_upsilon"].get<double>();
        Validator::validate_non_negative(drug.microfilaricidal_upsilon, "microfilaricidal_upsilon");
    }
    if (drug_json.contains("microfilaricidal_kappa")) {
        drug.microfilaricidal_kappa = drug_json["microfilaricidal_kappa"].get<double>();
        Validator::validate_positive(drug.microfilaricidal_kappa, "microfilaricidal_kappa");
    }
    if (drug_json.contains("embryostatic_lambda_max")) {
        drug.embryostatic_lambda_max = drug_json["embryostatic_lambda_max"].get<double>();
        Validator::validate_positive(drug.embryostatic_lambda_max, "embryostatic_lambda_max");
    }
    if (drug_json.contains("embryostatic_phi")) {
        drug.embryostatic_phi = drug_json["embryostatic_phi"].get<double>();
        Validator::validate_positive(drug.embryostatic_phi, "embryostatic_phi");
    }
    if (drug_json.contains("permanent_infertility")) {
        drug.permanent_infertility = drug_json["permanent_infertility"].get<double>();
        Validator::validate_probability(drug.permanent_infertility, "permanent_infertility");
    }
    if (drug_json.contains("allocation_proportion")) {
        drug.allocation_proportion = drug_json["allocation_proportion"].get<double>();
        Validator::validate_probability(drug.allocation_proportion, "allocation_proportion");
    }

    return drug;
}

std::vector<TreatmentParams> parse_treatments_from_json(const json& treatments_json) {
    std::vector<TreatmentParams> treatments;

    if (!treatments_json.is_array()) {
        throw std::invalid_argument("treatments must be an array");
    }

    for (const auto& treatment_json : treatments_json) {
        if (!treatment_json.contains("name")) {
            throw std::invalid_argument("Treatment must have a 'name' field");
        }

        std::string name = treatment_json["name"].get<std::string>();
        std::vector<double> application_times;
        bool has_interval_spec = treatment_json.contains("start_year") && treatment_json.contains("end_year");
        bool has_explicit_spec = treatment_json.contains("application_years");

        if (has_interval_spec && has_explicit_spec) {
            throw std::invalid_argument(
                "Treatment '" + name + "' cannot have both interval specification (start_year/end_year) "
                "and explicit application_years"
            );
        }

        if (!has_interval_spec && !has_explicit_spec) {
            throw std::invalid_argument(
                "Treatment '" + name + "' must have either (start_year + end_year + interval_years) "
                "or application_years"
            );
        }

        // Parse and validate drug(s) for a given treatment
        std::vector<DrugParams> drug_params;
        if (!treatment_json.contains("drugs")) {
            throw std::invalid_argument("Treatment '" + name + "' must have a 'drugs' field");
        }

        const auto& drugs_json = treatment_json["drugs"];
        if (!drugs_json.is_array()) {
            throw std::invalid_argument("drugs must be an array for treatment '" + name + "'");
        }

        double total_allocation = 0.0;
        for (const auto& drug_json : drugs_json) {
            DrugParams drug = parse_drug_from_json(drug_json);
            drug_params.push_back(drug);
            total_allocation += drug.allocation_proportion;
        }

        if (std::abs(total_allocation - 1.0) > 0) {
            throw std::invalid_argument(
                "Drug allocation_proportion must sum to 1.0 for treatment '" + name + 
                "', got " + std::to_string(total_allocation)
            );
        }

        // Parse and validate treatment parameters
        int min_age = treatment_json.value("min_age", 5);
        Validator::validate_non_negative_int(min_age, "min_age");

        double rho = treatment_json.value("rho", 0.3);
        Validator::validate_probability(rho, "rho");

        double coverage = treatment_json.value("coverage", 0.65);
        Validator::validate_probability(coverage, "coverage");

        double proportion_never_treated = treatment_json.value("proportion_never_treated", 0.0);
        Validator::validate_probability(proportion_never_treated, "proportion_never_treated");

        bool use_infection = treatment_json.value("use_infection", false);
        double infection_threshold = treatment_json.value("infection_threshold", 0.0);
        Validator::validate_non_negative(infection_threshold, "infection_threshold");

        TreatmentParams tp([&]() {
            if (has_interval_spec) {
                int start = treatment_json["start_year"].get<int>();
                int end = treatment_json["end_year"].get<int>();
                double interval = treatment_json.value("interval_years", 1.0);
                
                Validator::validate_non_negative_int(start, "start_year");
                Validator::validate_non_negative_int(end, "end_year");
                Validator::validate_positive(interval, "interval_years");

                if (start >= end) {
                    throw std::invalid_argument("start_year must be less than end_year");
                }

                return TreatmentParams(
                    start, end, interval, name, drug_params,
                    min_age, rho, proportion_never_treated, coverage,
                    use_infection, infection_threshold
                );
            } else {
                application_times = treatment_json["application_years"].get<std::vector<double>>();
                if (application_times.empty()) {
                    throw std::invalid_argument("application_years cannot be empty for treatment '" + name + "'");
                }
                for (double t : application_times) {
                    Validator::validate_non_negative(t, "application_year");
                }
                return TreatmentParams(
                    application_times, false, name, drug_params,
                    min_age, rho, proportion_never_treated, coverage,
                    use_infection, infection_threshold
                );
            }
        }());

        treatments.push_back(tp);
    }

    return treatments;
}

std::vector<VectorControlParams> parse_vector_control_from_json(const json& vc_json) {
    std::vector<VectorControlParams> vector_controls;

    if (!vc_json.is_array()) {
        throw std::invalid_argument("vector_control must be an array");
    }

    for (const auto& vc_item : vc_json) {
        if (!vc_item.contains("name")) {
            throw std::invalid_argument("Vector control intervention must have a 'name' field");
        }

        std::string name = vc_item["name"].get<std::string>();
        bool has_interval_spec = vc_item.contains("start_year") && vc_item.contains("end_year");
        bool has_explicit_spec = vc_item.contains("application_years");

        if (has_interval_spec && has_explicit_spec) {
            throw std::invalid_argument(
                "Vector control '" + name + "' cannot have both interval and explicit application_years"
            );
        }

        if (!has_interval_spec && !has_explicit_spec) {
            throw std::invalid_argument(
                "Vector control '" + name + "' must have either (start_year/end_year/interval_years) "
                "or application_years"
            );
        }

        std::vector<double> efficacies;
        if (!vc_item.contains("efficacies")) {
            throw std::invalid_argument("Vector control '" + name + "' must have an 'efficacies' field");
        }

        efficacies = vc_item["efficacies"].get<std::vector<double>>();
        if (efficacies.empty()) {
            throw std::invalid_argument("efficacies cannot be empty for vector control '" + name + "'");
        }

        for (double eff : efficacies) {
            Validator::validate_probability(eff, "efficacy");
        }

        int bounce_back_interval = vc_item.value("bounce_back_interval", 0);
        Validator::validate_non_negative_int(bounce_back_interval, "bounce_back_interval");

        VectorControlParams vcp([&]() {
            if (has_interval_spec) {
                int start = vc_item["start_year"].get<int>();
                int end = vc_item["end_year"].get<int>();
                double interval = vc_item.value("interval_years", 1.0);

                Validator::validate_non_negative_int(start, "start_year");
                Validator::validate_non_negative_int(end, "end_year");
                Validator::validate_positive(interval, "interval_years");

                if (start >= end) {
                    throw std::invalid_argument("start_year must be less than end_year in vector control '" + name + "'");
                }

                if (efficacies.size() != 1 && efficacies.size() != static_cast<size_t>(std::ceil((end - start) / interval))) {
                    throw std::invalid_argument(
                        "Number of efficacies must be 1 or match number of application times for vector control '" + name + "'"
                    );
                }

                return VectorControlParams(
                    start, end, interval, efficacies, bounce_back_interval, name
                );
            } else {
                std::vector<double> application_years = vc_item["application_years"].get<std::vector<double>>();
                if (application_years.empty()) {
                    throw std::invalid_argument("application_years cannot be empty for vector control '" + name + "'");
                }

                if (application_years.size() != efficacies.size()) {
                    throw std::invalid_argument(
                        "application_years and efficacies must have the same length for vector control '" + name + "'"
                    );
                }
                bool pre_converted_timesteps = vc_item.value("pre_converted_timesteps", false);

                for (double t : application_years) {
                    Validator::validate_non_negative(t, "application_year");
                }

                return VectorControlParams(
                    application_years, pre_converted_timesteps, efficacies, name
                );
            }
        }());

        vector_controls.push_back(vcp);
    }

    return vector_controls;
}

class ConfigParser {
public:
    static InputParams parse_config_file(const std::string& filename) {
        std::ifstream config_file(filename);
        if (!config_file.is_open()) {
            throw std::runtime_error("Cannot open config file: " + filename);
        }

        json config;
        try {
            config_file >> config;
        } catch (const json::parse_error& e) {
            throw std::runtime_error("JSON parse error in " + filename + ": " + e.what());
        }

        return parse_config_json(config);
    }

    static InputParams parse_config_json(const json& config) {
        Params params;

        if (config.contains("base")) {
            const auto& base_json = config["base"];
            
            if (base_json.contains("seed")) {
                params.base.seed = base_json["seed"].get<int>();
                Validator::validate_non_negative_int(params.base.seed, "seed");
            }
            if (base_json.contains("n_people")) {
                params.base.n_people = base_json["n_people"].get<int>();
                Validator::validate_positive_int(params.base.n_people, "n_people");
            }
            if (base_json.contains("k_E")) {
                params.base.k_E = base_json["k_E"].get<double>();
                std::vector<double> valid_ke = {0.2, 0.3, 0.4};
                Validator::validate_in_array(params.base.k_E, valid_ke, "kE");
            }
            if (base_json.contains("delta_time_days")) {
                params.base.delta_time_days = base_json["delta_time_days"].get<double>();
                Validator::validate_positive(params.base.delta_time_days, "delta_time_days");
            }
            if (base_json.contains("year_length_days")) {
                params.base.year_length_days = base_json["year_length_days"].get<double>();
                Validator::validate_positive(params.base.year_length_days, "year_length_days");
            }
            if (base_json.contains("month_length_days")) {
                params.base.month_length_days = base_json["month_length_days"].get<double>();
                Validator::validate_positive(params.base.month_length_days, "month_length_days");
            }
        }

        if (config.contains("blackfly")) {
            const auto& bf_json = config["blackfly"];
            
            if (bf_json.contains("bite_rate_per_person_per_year")) {
                params.blackfly.bite_rate_per_person_per_year = bf_json["bite_rate_per_person_per_year"].get<double>();
                Validator::validate_positive(params.blackfly.bite_rate_per_person_per_year, "bite_rate_per_person_per_year");
            }
            if (bf_json.contains("human_blood_index")) {
                params.blackfly.human_blood_index = bf_json["human_blood_index"].get<double>();
                Validator::validate_probability(params.blackfly.human_blood_index, "human_blood_index");
            }
            if (bf_json.contains("use_density_dependence")) {
                params.blackfly.use_density_dependence = bf_json["use_density_dependence"].get<bool>();
            }
            if (bf_json.contains("manual_density_dependence_params")) {
                params.blackfly.manual_density_dependence_params = bf_json["manual_density_dependence_params"].get<bool>();
                if (params.blackfly.manual_density_dependence_params) {
                    params.blackfly.c_h = bf_json["c_h"].get<double>();
                    Validator::validate_probability(params.blackfly.c_h, "c_h");
                    params.blackfly.delta_h_inf = bf_json["delta_h_inf"].get<double>();
                    Validator::validate_probability(params.blackfly.delta_h_inf, "delta_h_inf");
                    params.blackfly.delta_h_zero = bf_json["delta_h_zero"].get<double>();
                    Validator::validate_probability(params.blackfly.delta_h_zero, "delta_h_zero");
                }
            }
            if (bf_json.contains("k0")) {
                params.blackfly.k0 = bf_json["k0"].get<double>();
                Validator::validate_positive(params.blackfly.k0, "k0");
            }
            if (bf_json.contains("k1")) {
                params.blackfly.k1 = bf_json["k1"].get<double>();
                Validator::validate_positive(params.blackfly.k1, "k1");
            }
            if (bf_json.contains("x1")) {
                params.blackfly.x1 = bf_json["x1"].get<double>();
                Validator::validate_positive(params.blackfly.x1, "x1");
            }
            if (bf_json.contains("hbi_lb")) {
                params.blackfly.hbi_lb = bf_json["hbi_lb"].get<double>();
                Validator::validate_probability(params.blackfly.hbi_lb, "hbi_lb");
            }
            // Add other blackfly parameters as needed
        }

        if (config.contains("worms")) {
            const auto& worm_json = config["worms"];

            if (worm_json.contains("worm_age_stages")) {
                params.worms.worm_age_stages = worm_json["worm_age_stages"].get<double>();
                Validator::validate_positive(params.worms.worm_age_stages, "worm_age_stages");
            }
            if (worm_json.contains("q_w")) {
                params.worms.q_w = worm_json["q_w"].get<int>();
                Validator::validate_positive_int(params.worms.q_w, "q_w");
            }
            if (worm_json.contains("max_worm_age")) {
                params.worms.max_worm_age = worm_json["max_worm_age"].get<int>();
                Validator::validate_positive_int(params.worms.max_worm_age, "max_worm_age");
            }
            // Add other worm parameters as needed
        }

        if (config.contains("microfilaria")) {
            const auto& mf_json = config["microfilaria"];

            if (mf_json.contains("slope_kmf")) {
                params.mf.slope_kmf = mf_json["slope_kmf"].get<double>();
            }
            if (mf_json.contains("initial_kmf")) {
                params.mf.initial_kmf = mf_json["initial_kmf"].get<double>();
                Validator::validate_positive(params.mf.initial_kmf, "initial_kmf");
            }
            if (mf_json.contains("mf_age_stages")) {
                params.mf.mf_age_stages = mf_json["mf_age_stages"].get<int>();
                Validator::validate_positive_int(params.mf.mf_age_stages, "mf_age_stages");
            }
            if (mf_json.contains("q_m")) {
                params.mf.q_m = mf_json["q_m"].get<double>();
                Validator::validate_positive(params.mf.q_m, "q_m");
            }
            if (mf_json.contains("max_mf_age")) {
                params.mf.max_mf_age = mf_json["max_mf_age"].get<double>();
                Validator::validate_positive(params.mf.max_mf_age, "max_mf_age");
            }
            if (mf_json.contains("mf_move_rate")) {
                params.mf.mf_move_rate = mf_json["mf_move_rate"].get<double>();
                Validator::validate_positive(params.mf.mf_move_rate, "mf_move_rate");
            }
            if (mf_json.contains("use_kmf_const")) {
                params.mf.use_kmf_const = mf_json["use_kmf_const"].get<bool>();
            }
            if (mf_json.contains("kmf_const")) {
                params.mf.kmf_const = mf_json["kmf_const"].get<double>();
                Validator::validate_positive(params.mf.kmf_const, "kmf_const");
            }
            // Add other mf parameters as needed
        }

        if (config.contains("exposure")) {
            const auto& exp_json = config["exposure"];

            if (exp_json.contains("Q")) {
                params.exposure.Q = exp_json["Q"].get<double>();
                Validator::validate_positive(params.exposure.Q, "Q");
            }
            if (exp_json.contains("use_onchosim_exposure")) {
                params.exposure.use_onchosim_exposure = exp_json["use_onchosim_exposure"].get<bool>();
            }
            // Add other exposure parameters as needed
        }

        if (config.contains("human")) {
            const auto& human_json = config["human"];
            if (human_json.contains("min_skinsnip_age")) {
                params.human.min_skinsnip_age = human_json["min_skinsnip_age"].get<int>();
                Validator::validate_positive_int(params.human.min_skinsnip_age, "min_skinsnip_age");
            }
            if (human_json.contains("max_human_age")) {
                params.human.max_human_age = human_json["max_human_age"].get<int>();
                Validator::validate_positive_int(params.human.max_human_age, "max_human_age");
            }
            if (human_json.contains("mean_human_age")) {
                params.human.mean_human_age = human_json["mean_human_age"].get<int>();
                Validator::validate_positive_int(params.human.mean_human_age, "mean_human_age");
            }
            if (human_json.contains("skin_snip_weight")) {
                params.human.skin_snip_weight = human_json["skin_snip_weight"].get<int>();
                Validator::validate_positive_int(params.human.skin_snip_weight, "skin_snip_weight");
            }
            if (human_json.contains("skin_snip_number")) {
                params.human.skin_snip_number = human_json["skin_snip_number"].get<int>();
                Validator::validate_positive_int(params.human.skin_snip_number, "skin_snip_number");
            }
            if (human_json.contains("gender_ratio")) {
                params.human.gender_ratio = human_json["gender_ratio"].get<double>();
                Validator::validate_probability(params.human.gender_ratio, "gender_ratio");
            }
            if (human_json.contains("prop_serorevert_fast")) {
                params.human.prop_serorevert_fast = human_json["prop_serorevert_fast"].get<double>();
                Validator::validate_probability(params.human.prop_serorevert_fast, "prop_serorevert_fast");
            }
        }

        // TODO: Parse sequelae; 
        // if (!config.contains("sequelae") || config["sequelae"].is_null()) {
        //     params.sequelae_params = get_all_oncho_sequelae_params();
        // }

        // Setup InputParams object used by model
        std::vector<TreatmentParams> treatments;
        std::vector<VectorControlParams> vector_controls;

        if (config.contains("treatments")) {
            treatments = parse_treatments_from_json(config["treatments"]);
        }

        if (config.contains("vector_control")) {
            vector_controls = parse_vector_control_from_json(config["vector_control"]);
        }

        return InputParams(std::move(params), treatments, vector_controls);
    }
};

#endif
