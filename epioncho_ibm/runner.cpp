#include "model.hpp"

int main() {
    Params parameters;
    parameters.base.seed = 1;
    parameters.base.n_people = 500;
    parameters.blackfly.bite_rate_per_person_per_year = 1000;
    InputParams input_params;
    input_params.params = parameters;
    Model model(
        input_params
    );

    int timesteps = 365*100;
    for (int i = 0; i < timesteps; ++i) {
        model.advance_timestep();
    }

    printf("Total runtime: %f", model.overall_time);

    return 0;
}