#ifndef TOOLS_HPP
#define TOOLS_HPP

#include <iostream>

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

static std::string output_types_to_string(ModelOutputOption metric) {
    switch (metric) {
        case ModelOutputOption::population_size: return "population_size";
        case ModelOutputOption::mf_prevalence: return "mf_prevalence";
        case ModelOutputOption::true_ov16_seroprevalence: return "true_ov16_seroprevalence";
        case ModelOutputOption::adjusted_ov16_seroprevalence: return "adjusted_ov16_seroprevalence";
        case ModelOutputOption::mf_intensity: return "mf_intensity";
        case ModelOutputOption::worm_load: return "mean_worm_load";
        case ModelOutputOption::male_worm_load: return "mean_male_worm_load";
        case ModelOutputOption::female_worm_load: return "mean_female_worm_load";
        case ModelOutputOption::fertile_female_worm_load: return "mean_fertile_female_worm_load";
        case ModelOutputOption::infertile_female_worm_load: return "mean_infertile_female_worm_load";
        case ModelOutputOption::perm_sterile_female_worm_load: return "mean_perm_sterile_female_worm_load";
        case ModelOutputOption::compliance_percent: return "never_treated_percent";
        case ModelOutputOption::severe_itch_prevalence: return "severe_itch_prevalence";
        case ModelOutputOption::rsd_prevalence: return "rsd_prevalence";
        case ModelOutputOption::atrophy_prevalence: return "atrophy_prevalence";
        case ModelOutputOption::hanging_groin_prevalence: return "hanging_groin_prevalence";
        case ModelOutputOption::depigmentation_prevalence: return "depigmentation_prevalence";
        case ModelOutputOption::blindness_prevalence: return "blindness_prevalence";
        case ModelOutputOption::visual_impairment_prevalence: return "visual_impairment_prevalence";
        case ModelOutputOption::oae_prevalence: return "oae_prevalence";
        case ModelOutputOption::l3_per_blackfly: return "l3_per_blackfly";
        case ModelOutputOption::l3_prevalence_blackflies: return "l3_prevalence_blackflies";
        default: return "unknown";
    }
}

static ModelOutputOption string_to_output_types(const std::string& metric) {
    static const std::unordered_map<std::string, ModelOutputOption> lookup = {
        {"population_size", ModelOutputOption::population_size},
        {"mf_prevalence", ModelOutputOption::mf_prevalence},
        {"true_ov16_seroprevalence", ModelOutputOption::true_ov16_seroprevalence},
        {"adjusted_ov16_seroprevalence", ModelOutputOption::adjusted_ov16_seroprevalence},
        {"mf_intensity", ModelOutputOption::mf_intensity},
        {"mean_worm_load", ModelOutputOption::worm_load},
        {"mean_male_worm_load", ModelOutputOption::male_worm_load},
        {"mean_female_worm_load", ModelOutputOption::female_worm_load},
        {"mean_fertile_female_worm_load", ModelOutputOption::fertile_female_worm_load},
        {"mean_infertile_female_worm_load", ModelOutputOption::infertile_female_worm_load},
        {"mean_perm_sterile_female_worm_load", ModelOutputOption::perm_sterile_female_worm_load},
        {"never_treated_percent", ModelOutputOption::compliance_percent},
        {"severe_itch_prevalence", ModelOutputOption::severe_itch_prevalence},
        {"rsd_prevalence", ModelOutputOption::rsd_prevalence},
        {"atrophy_prevalence", ModelOutputOption::atrophy_prevalence},
        {"hanging_groin_prevalence", ModelOutputOption::hanging_groin_prevalence},
        {"depigmentation_prevalence", ModelOutputOption::depigmentation_prevalence},
        {"blindness_prevalence", ModelOutputOption::blindness_prevalence},
        {"visual_impairment_prevalence", ModelOutputOption::visual_impairment_prevalence},
        {"oae_prevalence", ModelOutputOption::oae_prevalence},
        {"l3_per_blackfly", ModelOutputOption::l3_per_blackfly},
        {"l3_prevalence_blackflies", ModelOutputOption::l3_prevalence_blackflies}
    };
    
    auto item = lookup.find(metric);
    if (item != lookup.end()) {
        return item->second;
    }
    throw std::invalid_argument(
        "Metric " + metric + " is not a valid option."
    );
}

#endif