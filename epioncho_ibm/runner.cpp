#include "model.hpp"
#include "model_outputs.hpp"
#include <cstdio>
#include <ctime>

static double get_elapsed_time(clock_t start_time) {
    return (double)(clock() - start_time) / CLOCKS_PER_SEC;
}

// TODO: possibly use cxxopts for parsing params
int main(int argc, char* argv[]) {
    bool verbose = false;
    for (int i = 1; i < argc; ++i) {
        if (std::string(argv[i]) == "--verbose") {
            verbose = true;
        }
    }

    int total_years = 100;
    clock_t start = clock();
    Params parameters;
    parameters.base.seed = 1;
    parameters.base.n_people = 400;
    parameters.blackfly.bite_rate_per_person_per_year = 1000;

    VectorControlParams vcp = VectorControlParams(
        total_years - 40, total_years - 29, 1,
        1, 1,
        "slash and clear"
    );
    std::vector<VectorControlParams> vcps = {vcp};


    TreatmentParams tp = TreatmentParams(
        total_years - 30, total_years - 20, 1,
        "IVM",
        DrugParams(),
        5,
        0.3,
        0.0,
        0.65
    );
    std::vector<TreatmentParams> treatments = {tp};

    InputParams input_params(parameters, treatments, vcps);
    Model model(input_params);

    OutputInfo oi = OutputInfo(
        total_years, 0, 1,
        0, 80,
        1900,
        {ModelOutputTypes::mf_intensity, ModelOutputTypes::mf_prevalence, ModelOutputTypes::population_size}
    );

    ModelOutputs mo = ModelOutputs(oi);

    const int timesteps = 365 * total_years;
    for (int i = 0; i < timesteps; ++i) {
        if (mo.should_update(model.state.current_timestep, model.state.timestep_years)) {
            mo.update(model.state);
        }
        model.advance_timestep(verbose);
    }

    mo.write("test_output.csv");

    printf("Model runtime: %f\n", model.overall_time);
    printf("Total runtime: %f\n", get_elapsed_time(start));
    return 0;
}
