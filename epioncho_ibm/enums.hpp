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

enum class SequelaeType {
    SevereItch,
    ReactiveSkinDisease,
    Atrophy,
    HangingGroin,
    Depigmentation,
    Blindness,
    // VisualImpairment, Just calculated as 1.2x Blindness prevalence
    OAE
};

enum class SequelaeProbTimeUnit {
    Year,
    Day
};

enum class ModelOutputOption {
    population_size,
    mf_prevalence,
    true_ov16_seroprevalence,
    adjusted_ov16_seroprevalence,
    mf_intensity,
    worm_load,
    female_worm_load,
    male_worm_load,
    fertile_female_worm_load,
    infertile_female_worm_load,
    perm_sterile_female_worm_load,
    compliance_percent,
    severe_itch_prevalence,
    rsd_prevalence,
    atrophy_prevalence,
    hanging_groin_prevalence,
    depigmentation_prevalence,
    blindness_prevalence,
    visual_impairment_prevalence,
    oae_prevalence
};