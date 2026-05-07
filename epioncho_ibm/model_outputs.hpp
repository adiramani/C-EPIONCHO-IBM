#ifndef MODEL_OUTPUTS_HPP
#define MODEL_OUTPUTS_HPP

#include "params.hpp"
#include "state.hpp"

class ModelOutputs {
    public:
        OutputInfo output_info;
        std::vector<double> outputs; // flat container for outputs
        int output_index = 0;

        explicit ModelOutputs(const OutputInfo output_info);

        void update(State& state);

        void write(const std::string& path) const;
};

#endif