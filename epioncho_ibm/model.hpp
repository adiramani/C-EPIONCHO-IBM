#ifndef MODEL_HPP
#define MODEL_HPP

#include "state.hpp"
#include <memory>

class Model {
public:
    // std::unique_ptr<State> state;
    State state;
    double overall_time = 0.0;
    double get_new_worms_time = 0.0;
    double worm_sex_time = 0.0;
    double new_l3_time = 0.0;
    double age_time = 0.0;
    double update_blackfly_time = 0.0;
    double process_death_time = 0.0;
    double calc_outputs_time = 0.0;
    double print_time = 0.0;

    Model(const InputParams& input_params);

    // Advance timestep function
    void advance_timestep();

};

#endif
