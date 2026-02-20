#include "model.hpp"

int main() {
    Params parameters;
    parameters.base.seed = 1;
    parameters.base.n_people = 400;
    InputParams input_params;
    input_params.params = parameters;
    Model model(
        input_params
    );

    int timesteps = 365*50;
    for (int i = 0; i < timesteps; ++i) {
        model.advance_timestep();
    }

    return 0;
}