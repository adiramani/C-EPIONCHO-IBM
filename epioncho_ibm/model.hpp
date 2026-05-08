#ifndef MODEL_HPP
#define MODEL_HPP

#include "state.hpp"
#include <json/json.h>

class Model {
private:
    void _print_yearly_stats(int current_timestep, double current_year, int num_vc_year, int num_mda_year);
public:
    State state;
    std::vector<std::unique_ptr<Treatment>> treatments;
    std::vector<std::unique_ptr<VectorControl>> vector_control;

    // Per-phase timing accumulators (wall-clock seconds)
    double overall_time = 0.0;
    double get_new_worms_time = 0.0;
    double worm_sex_time = 0.0;
    double new_l3_time = 0.0;
    double age_time = 0.0;
    double update_blackfly_time = 0.0;
    double process_death_time = 0.0;
    double calc_outputs_time = 0.0;
    double print_time = 0.0;

    explicit Model(const InputParams& input_params);

    void advance_timestep(bool verbose);

    InputParams params_from_json(const std::string& path);
};

#endif
