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
        // To differentiate multiple iterations of the same model in the output
        int model_iteration;
    public:
        explicit ModelOutputs(const OutputInfo output_info, int model_iteration = 1);

        bool should_update(int current_timestep, double timestep_years);

        void update(State& state);

        void write(const std::string& path, bool update = false) const;
};

#endif