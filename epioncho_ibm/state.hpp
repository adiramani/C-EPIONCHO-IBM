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
    double current_year = 0;
    Params params;
    double timestep_years = 0.0;
    People people;
    std::vector<Intervention*> current_interventions;
    

    State(const State&) = delete;
    State& operator=(const State&) = delete;
    State(State&&) = default;
    State& operator = (State&&) = default;

    explicit State(const Params& params);

    std::vector<int> get_sub_population(int age_start, int age_end);

    double mf_intensity(std::vector<int> filtered_individuals, int incubation_time_hours = 24);

    double mf_prevalence(std::vector<int> filtered_individuals, int incubation_time_hours = 24);
};

#endif
