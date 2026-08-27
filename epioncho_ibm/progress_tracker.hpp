#ifndef PROGRESS_TRACKER_HPP
#define PROGRESS_TRACKER_HPP

#include <atomic>
#include <ctime>
#include <cstdio>
#include <iomanip>
#include <sstream>
#include <omp.h>

class ProgressTracker {
private:
    std::atomic<int> total_seeds_completed{0};
    int total_repeats;
    int timesteps_per_year;
    std::time_t start_time;
    double total_years;

public:

    ProgressTracker(int repeats, double total_years_)
        : total_repeats(repeats),
          total_years(total_years_),
          start_time(std::time(nullptr))
    {}

    void report_year(int seed, int current_year) {
        int years_completed = total_seeds_completed * (int)total_years + current_year;
        int total_years_across_all = total_repeats * (int)total_years;

        double elapsed = std::difftime(std::time(nullptr), start_time);
        double rate = years_completed > 0 ? elapsed / years_completed : 0.0;
        double eta_seconds = rate * (total_years_across_all - years_completed);

        double progress_percent = total_years_across_all > 0 
            ? (double)years_completed / total_years_across_all * 100.0 
            : 0.0;

        auto format_time = [](double seconds) -> std::string {
            int hours = (int)seconds / 3600;
            int minutes = ((int)seconds % 3600) / 60;
            int secs = (int)seconds % 60;
            std::ostringstream oss;
            oss << std::setfill('0')
                << std::setw(2) << hours << ":"
                << std::setw(2) << minutes << ":"
                << std::setw(2) << secs;
            return oss.str();
        };

#pragma omp critical
        {
            printf(
                "\r[%3d/%3d seeds] %4d/%4d years [%6.2f%%] | "
                "Elapsed: %s | ETA: %s",
                total_seeds_completed.load(), total_repeats,
                years_completed, total_years_across_all,
                progress_percent,
                format_time(elapsed).c_str(),
                format_time(eta_seconds).c_str()
            );
            fflush(stdout);
        }
    }

    void report_complete(int seed) {
        total_seeds_completed++;
        int completed = total_seeds_completed.load();

#pragma omp critical
        {
            double elapsed = std::difftime(std::time(nullptr), start_time);
            printf("\n");
            printf(
                "  [%d/%d] Seed %d complete (elapsed: %.0f s)\n",
                completed, total_repeats, seed, elapsed
            );
            fflush(stdout);
        }
    }

    void report_final() {
        double total_elapsed = std::difftime(std::time(nullptr), start_time);
        double avg_time_per_seed = total_elapsed / total_repeats;

#pragma omp critical
        {
            printf("\n");
            printf("╔════════════════════════════════════════════════════════╗\n");
            printf("║                   SIMULATION COMPLETE                  ║\n");
            printf("║  Seeds:           %3d                                  ║\n", total_repeats);
            printf("║  Total years:     %3.0f per seed                      ║\n", total_years);
            printf("║  Total elapsed:   %.0f seconds (%.2f minutes)        ║\n", total_elapsed, total_elapsed / 60.0);
            printf("║  Avg per seed:    %.0f seconds                       ║\n", avg_time_per_seed);
            printf("╚════════════════════════════════════════════════════════╝\n");
            fflush(stdout);
        }
    }
};

#endif