enum class WormType {
    Fertile_Female,
    Infertile_Female,
    Sterilized_Female,
    Male
};

enum class InterventionType {
    VectorControl,
    MDA
};

enum class ModelOutputTypes {
    population_size,
    mf_prevalence,
    ov16_seroprevalence,
    mf_intensity,
    worm_load,
    female_worm_load,
    male_worm_load,
    fertile_female_worm_load,
    infertile_female_worm_load,
    perm_sterile_female_worm_load,
    compliance_percent
};