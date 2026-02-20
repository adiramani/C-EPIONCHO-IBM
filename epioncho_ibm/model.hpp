#ifndef MODEL_HPP
#define MODEL_HPP

#include "state.hpp"
#include <memory>

class Model {
public:
    // std::unique_ptr<State> state;
    State state;

    Model(const InputParams& input_params);

    // Advance timestep function
    void advance_timestep();

};

#endif
