#include "model_outputs.hpp"
#include <fstream>
#include <string>


ModelOutputs::ModelOutputs(const OutputInfo output_info) 
    : output_info(output_info) 
{
    n_storage_times = output_info.output_time_years.size();
    outputs.resize(output_info.output_time_years.size() * output_info.outputs_to_track.size());
}

bool ModelOutputs::should_update(int current_timestep, double timestep_years) {
    return round(output_info.output_time_years[output_index] / timestep_years) == current_timestep;
}

static double return_metric(StateSummary& s, ModelOutputOption metric) {
    switch(metric) {
        case ModelOutputOption::population_size: return s.total_individuals;
        case ModelOutputOption::mf_prevalence: return s.mean_mf_prevalence();
        case ModelOutputOption::true_ov16_seroprevalence: return s.mean_raw_ov16_seroprevalence();
        case ModelOutputOption::adjusted_ov16_seroprevalence: return s.mean_adjusted_ov16_seroprevalence();
        case ModelOutputOption::mf_intensity: return s.mean_skin_snip_mf_count();
        case ModelOutputOption::worm_load: return s.mean_worm_load();
        case ModelOutputOption::male_worm_load: return s.mean_male_worm_load();
        case ModelOutputOption::female_worm_load: return s.mean_female_worm_load();
        case ModelOutputOption::fertile_female_worm_load: return s.mean_fertile_female_worm_load();
        case ModelOutputOption::infertile_female_worm_load: return s.mean_infertile_female_worm_load();
        case ModelOutputOption::perm_sterile_female_worm_load: return s.mean_perm_sterile_female_worm_load();
        case ModelOutputOption::compliance_percent: return s.percent_never_treated();
        case ModelOutputOption::severe_itch_prevalence: return s.sequelae_prevalence(SequelaeType::SevereItch);
        case ModelOutputOption::rsd_prevalence: return s.sequelae_prevalence(SequelaeType::ReactiveSkinDisease);
        case ModelOutputOption::atrophy_prevalence: return s.sequelae_prevalence(SequelaeType::Atrophy);
        case ModelOutputOption::hanging_groin_prevalence: return s.sequelae_prevalence(SequelaeType::HangingGroin);
        case ModelOutputOption::depigmentation_prevalence: return s.sequelae_prevalence(SequelaeType::Depigmentation);
        case ModelOutputOption::blindness_prevalence: return s.sequelae_prevalence(SequelaeType::Blindness);
        case ModelOutputOption::visual_impairment_prevalence: return std::min(s.sequelae_prevalence(SequelaeType::Blindness) * 1.78, 1.0);
        case ModelOutputOption::oae_prevalence: return s.sequelae_prevalence(SequelaeType::OAE);
        default: return 0.0;
    }
}

void ModelOutputs::update(State& state) {
    std::vector<int> inds = state.get_sub_population(output_info.start_age, output_info.end_age);
    StateSummary summary = StateSummary(state.params.base.sequela_active);
    state.update_state_summary(inds, summary);
    for(std::size_t i = 0; i < output_info.outputs_to_track.size(); ++i) {
        ModelOutputOption curr_type = output_info.outputs_to_track[i];
        outputs[(i * n_storage_times) + output_index] = return_metric(summary, curr_type);
    }
    output_index++;
}

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
        default: return "unknown";
    }
}

void ModelOutputs::write(const std::string& path) const {
    std::ofstream f(path);

    // Header
    f << "output_year";
    for (const auto& output_type : output_info.outputs_to_track) {
        f << "," << output_types_to_string(output_type);
    }
    f << "\n";

    // Rows — one per output timepoint
    for (std::size_t t = 0; t < output_info.output_time_years.size(); ++t) {
        f << output_info.pretty_time_years[t];
        for (std::size_t i = 0; i < output_info.outputs_to_track.size(); ++i) {
            f << "," << outputs[(i * output_info.output_time_years.size()) + t];
        }
        f << "\n";
    }
}