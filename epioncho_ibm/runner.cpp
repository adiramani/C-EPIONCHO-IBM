#include "model.hpp"
#include "model_outputs.hpp"
#include "oncho_params.hpp"
#include <cstdio>
#include <ctime>
#include <sstream>
#include <omp.h>
#include <unistd.h> 
#include <sys/wait.h>
#include <signal.h>
#include <filesystem>
#include <fstream>
#include <algorithm>

namespace fs = std::filesystem;

static double get_elapsed_time(clock_t start_time) {
    return (double)(clock() - start_time) / CLOCKS_PER_SEC;
}

static double get_elapsed_time_omp(double start_time) {
    return omp_get_wtime() - start_time;
}
static std::vector<pid_t> global_pids;

void signal_handler(int sig) {
    std::cout << "\nInterrupt received. Killing child processes...\n";
    
    for (pid_t pid : global_pids) {
        kill(pid, SIGTERM);
    }
    
    for (pid_t pid : global_pids) {
        waitpid(pid, nullptr, 0);
    }
    
    std::cout << "Child processes terminated.\n";
    exit(1);
}

void merge_output_csvs(const std::string& tmp_output_folder, const std::string& final_output_path) {
    std::vector<std::string> csv_files;
    
    try {
        for (const auto& entry : fs::directory_iterator(tmp_output_folder)) {
            if (entry.is_regular_file()) {
                std::string filename = entry.path().filename().string();
                
                if (
                    filename.find("tmp_output_") == 0 && 
                    filename.find(".csv") == filename.length() - 4
                ) {
                    csv_files.push_back(entry.path().string());
                }
            }
        }
    } catch (const std::exception& e) {
        std::cerr << "Error reading directory: " << e.what() << "\n";
        return;
    }
    
    if (csv_files.empty()) {
        std::cerr << "No CSV files found in " << tmp_output_folder << "\n";
        return;
    }
    
    std::sort(csv_files.begin(), csv_files.end());
    
    std::cout << "Found " << csv_files.size() << " CSV files to merge\n";
    
    std::string merged_file = final_output_path;
    std::ofstream merged(merged_file);
    
    if (!merged.is_open()) {
        std::cerr << "Failed to open " << merged_file << " for writing\n";
        return;
    }
    
    bool first_file = true;
    
    for (const auto& csv_file : csv_files) {
        std::ifstream infile(csv_file);
        if (!infile.is_open()) {
            std::cerr << "Warning: Could not open " << csv_file << "\n";
            continue;
        }
        
        std::string line;
        bool skip_header = !first_file;
        
        while (std::getline(infile, line)) {
            if (skip_header) {
                skip_header = false;
                continue;
            }
            merged << line << "\n";
        }
        
        infile.close();
        first_file = false;
    }
    
    merged.close();
    std::cout << "Successfully merged to: " << merged_file << "\n";
}

// TODO: possibly use cxxopts for parsing params
int main(int argc, char* argv[]) {
    signal(SIGINT, signal_handler);
    signal(SIGTERM, signal_handler);

    bool verbose = false;
    double k_E = 0.3;
    double abr = 1000;
    int repeats = 10;
    int n_cores = 1;
    int total_years = 100;
    std::string output_folder = "";
    std::string temporary_output_folder = "";
    bool enable_timing = false;
    bool use_benin_int_history = false;
    int additional_treatment_years = 5;
    std::string additional_treatment_name = "bIVM";
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
        } else if (std::string(argv[i]) == "--tmp-output-folder") {
            if (i + 1 >= argc)
                    return 1;
            temporary_output_folder = std::string(argv[++i]);
        } else if (std::string(argv[i]) == "--total-years") {
            if (i + 1 >= argc)
                return 1;
            total_years = atoi(argv[++i]);
        } else if (std::string(argv[i]) == "--enable-timing") {
            enable_timing = true;
        } else if (std::string(argv[i]) == "--n-cores") {
            if (i + 1 >= argc)
                return 1;
            n_cores = atoi(argv[++i]);
        } else if (std::string(argv[i]) == "--benin") {
            use_benin_int_history = true;
        } else if (std::string(argv[i]) == "--additional-trt-yrs") {
            if (i + 1 >= argc)
                return 1;
            additional_treatment_years = atoi(argv[++i]);
        } else if (std::string(argv[i]) == "--additional-trt-name") {
            if (i + 1 >= argc)
                return 1;
            additional_treatment_name = std::string(argv[++i]);
        } 
    }
    if (temporary_output_folder.length() == 0) {
        temporary_output_folder = output_folder;
    }
    std::string country = "ghana";
    if (use_benin_int_history)
        country = "benin";
    std::cout << "Starting Simulations for kE: " << k_E << " ABR: " << abr << "\n";
    double overall_start = omp_get_wtime();
    int repeats_per_process = repeats / n_cores;

    for (int i = 0; i < n_cores; ++i) {
        pid_t pid = fork();
        if (pid == -1) {
            std::cerr << "fork() failed\n";
            exit(1);
        } else if (pid == 0) {
            for (int seed = 1; seed <= repeats_per_process; ++seed) {
                std::vector<ModelOutputs> local_outputs;
                
                int true_seed = (seed-1) * n_cores + i;
                clock_t start = clock();
                Params parameters;
                parameters.base.seed = true_seed;
                parameters.base.n_people = 500;
                parameters.base.k_E = k_E;
                parameters.blackfly.bite_rate_per_person_per_year = abr;
                parameters.blackfly.use_density_dependence = true;
                parameters.blackfly.x1 = 0.00008627075;
                parameters.blackfly.hbi_lb = 0.1319683;
                parameters.human.prop_serorevert_fast = 0.8026316;
                parameters.sequelae_params = get_all_oncho_sequelae_params();


                std::vector<double> vce_vals = {
                    0.80, 0.0, 0.0, 0.0, 0.80, 0.80, 0.80, 0.80, 0.80, 0.80,
                    0.64, 0.48, 0.32, 0.16, 0
                };
                std::vector<double> vc_years = {
                    86, 87, 88, 89, 90, 91, 92, 93, 94, 95,
                    96, 97, 98, 99, 100
                };
                VectorControlParams vcp_ghana = VectorControlParams(
                    vc_years, 
                    false,
                    vce_vals, 
                    "larviciding"
                );
                VectorControlParams vcp_benin = VectorControlParams(
                    88, 107, 
                    1, 0.80,
                    0, "larviciding"
                );


                VectorControlParams vcp_after = VectorControlParams(
                    127, 127 + additional_treatment_years + 1,
                    1, 0.80,
                    0, "slashnclear"
                );

                std::vector<VectorControlParams> vcps = {vcp_ghana, vcp_after};
                if (use_benin_int_history) {
                    vcps = {vcp_benin, vcp_after};
                }


                TreatmentParams tp_ghana_1 = TreatmentParams(
                    87, 106, 1,
                    "IVM",
                    DrugParamsIVM(),
                    5,
                    0.2,
                    0.0,
                    0.68
                );
                TreatmentParams tp_ghana_2 = TreatmentParams(
                    106, 107, 1,
                    "IVM",
                    DrugParamsIVM(),
                    5,
                    0.2,
                    0.0,
                    0.702
                );
                TreatmentParams tp_ghana_3 = TreatmentParams(
                    107, 108, 1,
                    "IVM",
                    DrugParamsIVM(),
                    5,
                    0.2,
                    0.0,
                    0.724
                );
                TreatmentParams tp_ghana_4 = TreatmentParams(
                    108, 109, 1,
                    "IVM",
                    DrugParamsIVM(),
                    5,
                    0.2,
                    0.0,
                    0.746
                );
                TreatmentParams tp_ghana_5 = TreatmentParams(
                    109, 110, 1,
                    "IVM",
                    DrugParamsIVM(),
                    5,
                    0.2,
                    0.0,
                    0.768
                );
                TreatmentParams tp_ghana_6 = TreatmentParams(
                    110, 120, 0.5,
                    "IVM",
                    DrugParamsIVM(),
                    5,
                    0.2,
                    0.0,
                    0.79
                );
                TreatmentParams tp_ghana_7 = TreatmentParams(
                    121, 127, 0.5,
                    "IVM",
                    DrugParamsIVM(),
                    5,
                    0.2,
                    0.0,
                    0.79
                );
                std::vector<TreatmentParams> treatments = {
                    tp_ghana_1, tp_ghana_2, tp_ghana_3, tp_ghana_4, tp_ghana_5,
                    tp_ghana_6, tp_ghana_7
                };

                TreatmentParams tp_benin_1 = TreatmentParams(
                    96, 103, 1,
                    "IVM",
                    DrugParamsIVM(),
                    5,
                    0.2,
                    0.0,
                    0.80
                );
                TreatmentParams tp_benin_2 = TreatmentParams(
                    103, 113, 0.5,
                    "IVM",
                    DrugParamsIVM(),
                    5,
                    0.2,
                    0.0,
                    0.80
                );
                TreatmentParams tp_benin_3 = TreatmentParams(
                    113, 120, 1,
                    "IVM",
                    DrugParamsIVM(),
                    5,
                    0.2,
                    0.0,
                    0.80
                );
                TreatmentParams tp_benin_4 = TreatmentParams(
                    121, 127, 1,
                    "IVM",
                    DrugParamsIVM(),
                    5,
                    0.2,
                    0.0,
                    0.80
                );
                if (use_benin_int_history) {
                    treatments = {
                        tp_benin_1, tp_benin_2, tp_benin_3, tp_benin_4
                    };
                }

                TreatmentParams tp_scenario_bIVM = TreatmentParams(
                    127, 127 + additional_treatment_years, 0.5,
                    "IVM",
                    DrugParamsIVM(),
                    5,
                    0.2,
                    0.0,
                    0.80
                );
                TreatmentParams tp_scenario_aMOX = TreatmentParams(
                    127, 127 + additional_treatment_years, 1.0,
                    "IVM",
                    DrugParamsMOX(),
                    4,
                    0.2,
                    0.0,
                    0.80
                );
                TreatmentParams tp_scenario_bMOX = TreatmentParams(
                    127, 127 + additional_treatment_years, 0.5,
                    "MOX",
                    DrugParamsMOX(),
                    4,
                    0.2,
                    0.0,
                    0.80
                );
                
                if (additional_treatment_name == "bIVM") {
                    treatments.push_back(tp_scenario_bIVM);
                } else if (additional_treatment_name == "aMOX") {
                    treatments.push_back(tp_scenario_aMOX);
                    parameters.base.delta_time_days = 0.5;
                } else if (additional_treatment_name == "bMOX") {
                    treatments.push_back(tp_scenario_bMOX);
                    parameters.base.delta_time_days = 0.5;
                }

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
                    0, 0, 0, 3, 10, 15, 20, 30, 40, 50, 60, 70
                };
                std::vector<int> age_ends = {
                    81, 66, 3, 10, 15, 20, 30, 40, 50, 60, 70, 81
                };
                std::vector<ModelOutputs> all_model_outputs;

                for (size_t a = 0; a < age_starts.size(); ++a) {
                    double interval = 1.0;
                    if (age_starts[a] == 0 && age_ends[a] == 81)
                        interval = 1.0 / 8.0;
                    all_model_outputs.push_back(
                        ModelOutputs(
                            OutputInfo(
                                total_years, 65.0, interval,
                                age_starts[a], age_ends[a],
                                1900, 0.80, 0.99,
                                all_outputs
                            ),
                            true_seed
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
                for (auto& mo : all_model_outputs) {
                    local_outputs.push_back(mo);
                }
                if (enable_timing) {
                    printf("Model runtime: %f\n", model.overall_time);
                    printf("Total runtime: %f\n", get_elapsed_time(start));
                }

                std::ostringstream oss;

                oss << temporary_output_folder << "tmp_output_abr_" << abr << "_kE_" << k_E << "_trt_" << additional_treatment_name << "_yrs_" << additional_treatment_years << "_" << country << "_vc_" << true_seed << ".csv";

                printf("Writing output to %s\n", oss.str().c_str());


                int iter = 0;
                for (auto& mo : local_outputs) {
                    mo.write(oss.str(), iter > 0);
                    iter++;
                }
                printf("Total runtime seed %d: %f\n", true_seed, get_elapsed_time(start));
            }
            exit(0);
        }
        else {
            global_pids.push_back(pid);
        }
    }

    int completed = 0;
    for (pid_t pid : global_pids) {
        int status;
        waitpid(pid, &status, 0);
        
        if (WIFEXITED(status)) {
            completed++;
            std::cout << "Child " << completed << "/" << n_cores << " completed\n";
        } else {
            std::cerr << "Child process failed\n";
        }
    }

    std::cout << "All children done.\n";
    printf("Overall Model Runtime: %f\n", get_elapsed_time_omp(overall_start));
    std::cout << "Merging output files...\n";

    std::string final_output_path = (
        output_folder + "final_output_abr_" + std::to_string((int)abr) +
        "_kE_" + std::to_string(k_E) + "_trt_" + additional_treatment_name + "_yrs_" + 
        std::to_string(additional_treatment_years) + "_" + country + "_vc.csv"
    );

    merge_output_csvs(temporary_output_folder, final_output_path);

    return 0;
}
