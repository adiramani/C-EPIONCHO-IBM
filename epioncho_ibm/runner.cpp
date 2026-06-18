#include "model.hpp"
#include "model_outputs.hpp"
#include "oncho_params.hpp"
#include <cstdio>
#include <ctime>
#include <sstream>
#include <omp.h>

static double get_elapsed_time(clock_t start_time) {
    return (double)(clock() - start_time) / CLOCKS_PER_SEC;
}

static double get_elapsed_time_omp(double start_time) {
    return omp_get_wtime() - start_time;
}

// TODO: possibly use cxxopts for parsing params
int main(int argc, char* argv[]) {
    bool verbose = false;
    double k_E = 0.3;
    // double seed = 1;
    double abr = 1000;
    int repeats = 10;
    int total_years = 100;
    std::string output_folder = "";
    bool enable_timing = false;
    for (int i = 1; i < argc; ++i) {
        if (std::string(argv[i]) == "--verbose") {
            verbose = true;
        } else if (std::string(argv[i]) == "--kE") {
            if (i + 1 >= argc)
                return 1;
            k_E = atof(argv[++i]);
        } else if (std::string(argv[i]) == "--abr") {
            if (i + 1 >= argc)
                return 1;
            abr = atoi(argv[++i]);
        } else if (std::string(argv[i]) == "--repeats") {
            if (i + 1 >= argc)
                return 1;
            repeats = atoi(argv[++i]);
        } else if (std::string(argv[i]) == "--output-folder") {
            if (i + 1 >= argc)
                return 1;
            output_folder = std::string(argv[++i]);
        } else if (std::string(argv[i]) == "--total-years") {
            if (i + 1 >= argc)
                return 1;
            total_years = atoi(argv[++i]);
        } else if (std::string(argv[i]) == "--enable-timing") {
            enable_timing = true;
        }
    }
    std::cout << "Starting Simulations for kE: " << k_E << " ABR: " << abr << "\n";
    std::vector<ModelOutputs> final_model_outputs;
    double overall_start = omp_get_wtime();
#pragma omp parallel for
    for (int seed = 1; seed <= repeats; ++seed) {
        clock_t start = clock();
        Params parameters;
        parameters.base.seed = seed;
        parameters.base.n_people = 400;
        parameters.base.k_E = k_E;
        parameters.blackfly.bite_rate_per_person_per_year = abr;
        parameters.sequelae_params = get_all_oncho_sequelae_params();

        VectorControlParams vcp = VectorControlParams(
            total_years - 40, total_years - 29, 1,
            1, 1,
            "slash and clear"
        );
        std::vector<VectorControlParams> vcps = {};


        TreatmentParams tp_aIVM = TreatmentParams(
            total_years - 40, total_years - 20, 1,
            "IVM",
            DrugParamsIVM(),
            5,
            0.3,
            0.0,
            0.65
        );

        TreatmentParams tp_bIVM = TreatmentParams(
            total_years - 40, total_years - 20, 0.5,
            "IVM",
            DrugParamsIVM(),
            5,
            0.3,
            0.0,
            0.65
        );

        TreatmentParams tp_aMOX = TreatmentParams(
            total_years - 40, total_years - 20, 1,
            "MOX",
            DrugParamsMOX(),
            4,
            0.3,
            0.0,
            0.65
        );

        TreatmentParams tp_bMOX = TreatmentParams(
            total_years - 40, total_years - 20, 0.5,
            "MOX",
            DrugParamsMOX(),
            4,
            0.3,
            0.0,
            0.65
        );

        std::vector<TreatmentParams> treatments = {tp_aMOX};
        parameters.base.delta_time_days = 0.5;

        const int total_timesteps = (parameters.base.year_length_days / parameters.base.delta_time_days) * total_years;


        InputParams input_params(
            std::move(parameters), 
            treatments, 
            vcps
        );
        Model model(std::move(input_params), enable_timing);


        std::vector<ModelOutputOption> all_outputs = {
            ModelOutputOption::mf_intensity, 
            ModelOutputOption::mf_prevalence, 
            ModelOutputOption::population_size, 
            ModelOutputOption::true_ov16_seroprevalence,
            ModelOutputOption::adjusted_ov16_seroprevalence,
            ModelOutputOption::worm_load,
            ModelOutputOption::female_worm_load,
            ModelOutputOption::male_worm_load,
            ModelOutputOption::fertile_female_worm_load,
            ModelOutputOption::infertile_female_worm_load,
            ModelOutputOption::perm_sterile_female_worm_load,
            ModelOutputOption::compliance_percent,
            ModelOutputOption::severe_itch_prevalence,
            ModelOutputOption::rsd_prevalence,
            ModelOutputOption::atrophy_prevalence,
            ModelOutputOption::hanging_groin_prevalence,
            ModelOutputOption::depigmentation_prevalence,
            ModelOutputOption::blindness_prevalence,
            ModelOutputOption::visual_impairment_prevalence,
            ModelOutputOption::oae_prevalence,
            ModelOutputOption::l3_per_blackfly,
            ModelOutputOption::l3_prevalence_blackflies
        };

        std::vector<int> age_starts = {
            0, 0, 5, 10, 15, 20, 30, 40, 50, 60, 70
        };
        std::vector<int> age_ends = {
            81, 5, 10, 15, 20, 30, 40, 50, 60, 70, 81
        };
        std::vector<ModelOutputs> all_model_outputs;

        for (size_t a = 0; a < age_starts.size(); ++a) {
            all_model_outputs.push_back(
                ModelOutputs(
                    OutputInfo(
                        total_years, 0, 1,
                        age_starts[a], age_ends[a],
                        0,
                        all_outputs
                    ),
                    seed
                )
            );
        }

        for (int i = 0; i < total_timesteps; ++i) {
            for (auto& mo : all_model_outputs) {
                if (mo.should_update(model.state.current_timestep, model.state.timestep_years))
                    mo.update(model.state);
            }
            model.advance_timestep(verbose);
        }
#pragma omp critical
{
        for (auto& mo : all_model_outputs)
            final_model_outputs.push_back(mo);
}
        if (enable_timing) {
            printf("Model runtime: %f\n", model.overall_time);
            printf("Total runtime: %f\n", get_elapsed_time(start));
        }
    }

    std::ostringstream oss;
    oss << output_folder << "test_output_abr_" << abr << "_kE_" << k_E << ".csv";

    int iter = 0;
    for (auto& mo : final_model_outputs) {
        mo.write(oss.str(), iter > 0);
        iter++;
    }

    printf("Overall Model Runtime: %f\n", get_elapsed_time_omp(overall_start));

    return 0;
}
