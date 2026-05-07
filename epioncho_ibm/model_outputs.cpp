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

static const std::unordered_map<ModelOutputTypes, std::function<double(State&, std::vector<int>&)>> metrics_map = {
    { ModelOutputTypes::population_size, [](State& s, std::vector<int>& inds) { return (double)inds.size(); }},
    { ModelOutputTypes::mf_intensity, [](State& s, std::vector<int>& inds) { return s.mf_intensity(inds); }},
    { ModelOutputTypes::mf_prevalence, [](State& s, std::vector<int>& inds) { return s.mf_prevalence(inds); }},
};

void ModelOutputs::update(State& state) {
    std::vector<int> inds = state.get_sub_population(output_info.start_age, output_info.end_age);
    for(std::size_t i = 0; i < output_info.outputs_to_track.size(); ++i) {
        ModelOutputTypes curr_type = output_info.outputs_to_track[i];
        auto metric_iter = metrics_map.find(output_info.outputs_to_track[i]);
        if (metric_iter != metrics_map.end()) {
            outputs[(i * n_storage_times) + output_index] = metric_iter->second(state, inds);
        }
    }
    output_index++;
}

static std::string output_types_to_string(ModelOutputTypes t) {
    switch (t) {
        case ModelOutputTypes::population_size: return "population_size";
        case ModelOutputTypes::mf_intensity:   return "mf_intensity";
        case ModelOutputTypes::mf_prevalence:  return "mf_prevalence";
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