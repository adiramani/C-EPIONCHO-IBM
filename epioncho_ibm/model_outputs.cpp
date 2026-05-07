#include "model_outputs.hpp"
#include <fstream>
#include <string>


ModelOutputs::ModelOutputs(const OutputInfo output_info) 
    : output_info(output_info) 
{
    outputs.resize(output_info.output_time_years.size() * output_info.outputs_to_track.size());
}

void ModelOutputs::update(State& state) {
    std::vector<int> inds = state.get_sub_population(output_info.start_age, output_info.end_age);
    for(std::size_t i = 0; i < output_info.outputs_to_track.size(); ++i) {
        ModelOutputTypes curr_type = output_info.outputs_to_track[i];
        if (curr_type == ModelOutputTypes::mf_intensity) {
            // Flat array structure, so i * array index is the location of output index "output_index" for output type i
            outputs[(i * output_info.output_time_years.size()) + output_index] = state.mf_intensity(inds);
        }
        if (curr_type == ModelOutputTypes::mf_prevalence) {
            // Flat array structure, so i * array index is the location of output index "output_index" for output type i
            outputs[(i * output_info.output_time_years.size()) + output_index] = state.mf_prevalence(inds);
        }
    }
    output_index++;
}

static std::string output_types_to_string(ModelOutputTypes t) {
    switch (t) {
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
        f << output_info.output_time_years[t];
        for (std::size_t i = 0; i < output_info.outputs_to_track.size(); ++i) {
            f << "," << outputs[(i * output_info.output_time_years.size()) + t];
        }
        f << "\n";
    }
}