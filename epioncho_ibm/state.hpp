#ifndef STATE_HPP
#define STATE_HPP

#include "params.hpp"
#include "people.hpp"
#include "interventions.hpp"
#include <memory>
#include <vector>

class State {
public:
    std::mt19937 generator;

    int current_timestep = 0;
    double timestep_years;

    People people;
    std::vector<std::unique_ptr<Treatment>> treatments;
    std::vector<std::unique_ptr<VectorControl>> vector_control;
    Params params;

    State(const InputParams& input_params);

    // // Deleted copy constructor and assignment to avoid accidental copies
    // State(const State&) = delete;
    // State& operator=(const State&) = delete;

    // // Move constructor and assignment
    // State(State&&) = default;
    // State& operator=(State&&) = default;
};

#endif
