#ifndef MODEL_OUTPUTS_HPP
#define MODEL_OUTPUTS_HPP

#include "params.hpp"
#include "state.hpp"

class ModelOutputs {
    private:
        OutputInfo output_info;
        int n_storage_times;
        std::vector<double> outputs; // flat container for outputs
        int output_index = 0;
    public:
        explicit ModelOutputs(const OutputInfo output_info);

        bool should_update(int current_timestep, double timestep_years);

        void update(State& state);

        void write(const std::string& path) const;
};

#endif