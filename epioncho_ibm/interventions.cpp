# include "interventions.hpp"

Treatment::Treatment(const TreatmentParams& params_)
    : Intervention(params_), treatment_params(params_) 
{}

VectorControl::VectorControl(const VectorControlParams& params_)
    : Intervention(params_), vc_params(params_) 
{}