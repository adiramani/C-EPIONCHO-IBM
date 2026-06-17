#ifndef STATE_HPP
#define STATE_HPP

#include "params.hpp"
#include "people.hpp"
#include "interventions.hpp"
#include <memory>
#include <vector>

struct StateSummary {
    int total_individuals = 0;
    double pop_raw_microfilarial_load = 0.0;
    double pop_skin_snip_microfilarial_load = 0.0;
    double pop_male_worm_load = 0.0;
    double pop_fertile_female_worm_load = 0.0;
    double pop_infertile_female_worm_load = 0.0;
    double pop_perm_sterile_female_worm_load = 0.0;
    double l3_per_blackfly = 0.0;
    double l3_prevalence_blackflies = 0.0;

    int skin_snip_positives = 0;
    int raw_ov16_seropositives = 0;
    int adjusted_ov16_seropositives = 0;
    int num_never_treated = 0;
    std::vector<SequelaeType> active_sequelae = {};
    std::vector<int> sequelae_positives = {};

    StateSummary() = default;
    StateSummary(std::vector<SequelaeType> active_sequelae)
    : active_sequelae(active_sequelae),
      sequelae_positives(active_sequelae.size(), 0)
    {}

    double mean_raw_mf_count() const { return total_individuals > 0 ? pop_raw_microfilarial_load / total_individuals : 0.0; }
    double mean_skin_snip_mf_count() const { return total_individuals > 0 ? pop_skin_snip_microfilarial_load / total_individuals : 0.0; }
    double mean_worm_load() const { return total_individuals > 0 ? (pop_male_worm_load + pop_fertile_female_worm_load + pop_infertile_female_worm_load + pop_perm_sterile_female_worm_load) / total_individuals : 0.0; }
    double mean_male_worm_load() const { return total_individuals > 0 ? pop_male_worm_load / total_individuals : 0.0; }
    double mean_female_worm_load() const { return total_individuals > 0 ? (pop_fertile_female_worm_load + pop_infertile_female_worm_load + pop_perm_sterile_female_worm_load) / total_individuals : 0.0; }
    double mean_fertile_female_worm_load() const { return total_individuals > 0 ? pop_fertile_female_worm_load / total_individuals : 0.0; }
    double mean_infertile_female_worm_load() const { return total_individuals > 0 ? pop_infertile_female_worm_load / total_individuals : 0.0; }
    double mean_perm_sterile_female_worm_load() const { return total_individuals > 0 ? pop_perm_sterile_female_worm_load / total_individuals : 0.0; }
    double mean_mf_prevalence() const { return total_individuals > 0 ? (double)skin_snip_positives / total_individuals : 0.0; }
    double mean_raw_ov16_seroprevalence() const { return total_individuals > 0 ? (double)raw_ov16_seropositives / total_individuals : 0.0; }
    double mean_adjusted_ov16_seroprevalence() const { return total_individuals > 0 ? (double)adjusted_ov16_seropositives / total_individuals : 0.0; }
    double percent_never_treated() const { return total_individuals > 0 ? (double)num_never_treated / total_individuals : 0.0; }
    double mean_l3_per_blackfly() const { return l3_per_blackfly; }
    double mean_l3_prevalence_blackflies() const { return l3_prevalence_blackflies; }

    double sequelae_prevalence(SequelaeType st) const { 
        for (size_t s = 0; s < active_sequelae.size(); ++s) {
            if (active_sequelae[s] == st)
                return total_individuals > 0 ? (double)sequelae_positives[s] / total_individuals : 0.0; 
        }
        return 0.0;
    };
};

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

        explicit State(Params params);

        std::vector<int> get_sub_population(int age_start, int age_end);

        double mf_intensity(const std::vector<int>& filtered_individuals, int incubation_time_hours = 24);

        double mf_prevalence(const std::vector<int>& filtered_individuals, int incubation_time_hours = 24);

        void update_state_summary(const std::vector<int>& filtered_individuals, StateSummary& summary, double sens = 0.80, double spec = 0.99);
};

#endif
