#include "state.hpp"

State::State(const InputParams& input_params)
    : params(input_params.params),
      people(input_params.params.base.n_people),
      timestep_years(input_params.params.base.delta_time_days / input_params.params.base.year_length_days)
{
  generator = std::mt19937(params.base.seed);


  people.initialize_from_params(
    generator,
    timestep_years,
    input_params.params.human.mean_human_age,
    input_params.params.human.max_human_age,
    input_params.params.base.k_E,
    input_params.params.worms,
    input_params.params.mf,
    input_params.params.blackfly
  );
  
  // for (int i = 0; i < input_params.treatments.size(); i++) {
  //   treatments.push_back(std::make_unique<Treatment>(input_params.treatments[i]));
  // } 
  // for (int i = 0; i < input_params.vector_control.size(); i++) {
  //   vector_control.push_back(std::make_unique<VectorControl>(input_params.vector_control[i]));
  // }
}
