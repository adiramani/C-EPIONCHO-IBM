#ifndef VECTOR_HPP
#define VECTOR_HPP

#include <vector>

class Vector {
    public:
        double l1;
        double l2;
        double l3;

        double delta_v0;
        double c_v; // Severity of density-dependent 
        double v_1; // per capita development rate from L1 to L2
        double v_2; // per capita development rate from L2 to L3
        double mu_v; // per-capita mortality rate of blackflies
        double alpha_v; // per-capita mf induced mortality of blackflies
        double tou_v; // delay before L1 transition to L2

        Vector() = default;

        Vector(
            double delta_v0_,
            double c_v_,
            double v_1_,
            double v_2_,
            double mu_v_,
            double alpha_v_,
            double tou_v_,
            double init_l1,
            double init_l2,
            double init_l3
        );

        virtual void process_death();
};

class Blackfly: public Vector {
    public:

        struct Delay {
            double l1;
            double microfilariae;
            double exposure;
        };

        int delay_index;
        // [(l1 delay, mf delay, exposure delay), (...), ...]
        std::vector<Delay> delays;

        void update_delay_index(double l1, double mf, double exposure);

        double calc_density_dependence_i(double mf, double exposure);

        void calc_L1(double timestep_years, double beta, double exposure, float mf);

        void calc_L2(double timestep_years);

        void calc_L3(double a_H, double g, double mu_l3);

        Blackfly() = default;

        Blackfly(
            double delta_v0_,
            double c_v_,
            double v_1_,
            double v_2_,
            double mu_v_,
            double alpha_v_,
            double tou_v_,
            double init_l1,
            double init_l2,
            double init_l3,
            double init_exposure_heterogeneity
        );

        void process_death() override;
};

#endif