#include "model.hpp"
#include "model_outputs.hpp"
#include "config_parser.hpp"
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
#include <string>
#include <iostream>

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

void check_output_paths(const std::string& tmp_output_folder, const std::string& final_output_folder) {
    if (!fs::exists(tmp_output_folder)) {
        std::cout << "Directory " << tmp_output_folder << " does not exist. Creating it now.\n";
        fs::create_directories(tmp_output_folder);
    }
    if (!fs::exists(final_output_folder)) {
        std::cout << "Directory " << final_output_folder << " does not exist. Creating it now.\n";
        fs::create_directories(final_output_folder);
    }
}

void merge_output_csvs(const std::string& tmp_output_folder, const std::string& final_output_path, bool delete_temp_after_processing) {
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
    if (delete_temp_after_processing) {
        std::cout << "Deleting temporary folder: " << tmp_output_folder << " ... ";
        try {
            fs::remove_all(tmp_output_folder);
        } catch (const std::exception& e) {
            std::cerr << "Error reading directory: " << e.what() << "\n";
            return;
        }
        std::cout << "Finished\n";
    }
}

// TODO: possibly use cxxopts for parsing params
int main(int argc, char* argv[]) {
    signal(SIGINT, signal_handler);
    signal(SIGTERM, signal_handler);

    std::string config_path = "";
    bool verbose = false;
    double k_E = 0.3;
    double abr = 1000.0;
    int repeats = 10;
    int total_years = 100;
    std::string output_folder = "test_final_output/";
    std::string temporary_output_folder = "temp_output/";
    bool delete_temp_after_processing = true;
    bool enable_timing = false;
    bool onchosim_exposure = false;
    int n_cores = 1;
    for (int i = 1; i < argc; ++i) {
        if (std::string(argv[i]) == "--config") {
            if (i + 1 >= argc)
                return 1;
            config_path = std::string(argv[++i]);
        } else if (std::string(argv[i]) == "--verbose") {
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
        } else if (std::string(argv[i]) == "--keep-tmp-folder") {
            delete_temp_after_processing = false;
        } else if (std::string(argv[i]) == "--total-years") {
            if (i + 1 >= argc)
                return 1;
            total_years = atoi(argv[++i]);
        } else if (std::string(argv[i]) == "--enable-timing") {
            enable_timing = true;
        } else if (std::string(argv[i]) == "--exposure") {
            if (i + 1 >= argc)
                return 1;
            onchosim_exposure = std::string(argv[++i]) == "onchosim";
        } else if (std::string(argv[i]) == "--n-cores") {
            if (i + 1 >= argc)
                return 1;
            n_cores = atoi(argv[++i]);
        }
    }
    if (temporary_output_folder.length() == 0) {
        temporary_output_folder = output_folder;
    }

    check_output_paths(temporary_output_folder, output_folder);
    double overall_start = omp_get_wtime();

    FullModelConfig model_config;
    bool loaded_params = false;
    std::string simulation_name = "Command line args";
    std::string simulation_name_for_output = simulation_name;
    if (config_path != "") {
        std::cout << "Loading configuration from: " << config_path << " ... ";
        try {
            model_config = ConfigParser::parse_config_file(config_path);
            loaded_params = true;
        } catch (const std::exception& e) {
            std::cerr << "Error loading config: " << e.what() << "\n";
            return 1;
        }
        std::cout << "Finished.\n";
        if (model_config.runtime_info.num_cores > 0) {
            n_cores = model_config.runtime_info.num_cores;
            std::cout << "Using num_cores set in config: " << n_cores << ".\n";
        }
        if (model_config.runtime_info.num_repeats > 0) {
            repeats = model_config.runtime_info.num_repeats;
        }
        simulation_name = model_config.simulation_name;
        simulation_name_for_output = simulation_name;
    }
    std::replace(simulation_name_for_output.begin(), simulation_name_for_output.end(), ' ', '_');
    std::transform(simulation_name_for_output.begin(), simulation_name_for_output.end(), simulation_name_for_output.begin(), ::tolower);

    int repeats_per_process = repeats / n_cores;
    std::cout << "Starting " << simulation_name << " Simulations with " << repeats << " repeats and " << n_cores << " cores.\n";


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

                InputParams input_params_seed;
                std::vector<ModelOutputs> all_model_outputs;
                if (!loaded_params) {
                    Params parameters;
                    parameters.base.seed = true_seed;
                    parameters.base.k_E = k_E;
                    parameters.blackfly.bite_rate_per_person_per_year = abr;
                    parameters.exposure.use_onchosim_exposure = onchosim_exposure;

                    input_params_seed = InputParams(
                        parameters, 
                        {}, 
                        {}
                    );
                    std::vector<ModelOutputOption> all_outputs = {
                        ModelOutputOption::mf_intensity, 
                        ModelOutputOption::mf_prevalence, 
                        ModelOutputOption::population_size, 
                        ModelOutputOption::true_ov16_seroprevalence,
                        ModelOutputOption::adjusted_ov16_seroprevalence,
                        ModelOutputOption::l3_per_blackfly,
                        ModelOutputOption::l3_prevalence_blackflies
                    };

                    std::vector<int> age_starts = {
                        5, 0, 5, 10, 15, 20, 30, 40, 50, 60, 70
                    };
                    std::vector<int> age_ends = {
                        81, 5, 10, 15, 20, 30, 40, 50, 60, 70, 81
                    };

                    for (size_t a = 0; a < age_starts.size(); ++a) {
                        double interval = 1.0;
                        all_model_outputs.push_back(
                            ModelOutputs(
                                OutputInfo(
                                    total_years, 50.0, interval,
                                    age_starts[a], age_ends[a],
                                    1900, 0.80, 0.99,
                                    all_outputs
                                ),
                                true_seed
                            )
                        );
                    }
                } else {
                    Params loaded_parameters = model_config.input_params.params;
                    loaded_parameters.base.seed = true_seed;
                    input_params_seed = InputParams(
                        loaded_parameters,
                        model_config.input_params.treatments,
                        model_config.input_params.vector_control
                    );
                    for (auto& output_info : model_config.output_infos) {
                        all_model_outputs.push_back(
                            ModelOutputs(
                                output_info,
                                true_seed
                            )
                        );
                    }
                    total_years = model_config.runtime_info.total_years;
                }

                const int total_timesteps = (input_params_seed.params.base.year_length_days / input_params_seed.params.base.delta_time_days) * total_years;
                
                Model model(std::move(input_params_seed), enable_timing);

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

                oss << temporary_output_folder << "tmp_output_" << simulation_name_for_output << "_" << true_seed << ".csv";
                printf("Writing output to %s\n", oss.str().c_str());
                int iter = 0;
                for (auto& mo : local_outputs) {
                    mo.write(oss.str(), iter > 0);
                    iter++;
                }
                printf("Total runtime seed %d: %f\n", true_seed, get_elapsed_time(start));
            }
            exit(0);
        } else {
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
        output_folder + "final_output_" + simulation_name_for_output + ".csv"
    );

    merge_output_csvs(temporary_output_folder, final_output_path, delete_temp_after_processing);

    return 0;
}
