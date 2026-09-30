import React, { useState } from 'react';
import './ConfigForm.css';

const ConfigForm = ({ onConfigGenerated, onModelResultsGenerated }) => {
  const [showAdvancedSettings, setShowAdvancedSettings] = useState({
    show_advanced_settings: false,
  });

  const getDefaultConfig = () => ({
    name: "Config",
    runtime_params: {
      total_years: 100,
      num_cores: 1,
      num_repeats: 1
    },
    model_params: {
      base: {
        seed: 1,
        n_people: 1000,
        k_E: 0.3,
        delta_time_days: 1.0,
        year_length_days: 366,
        month_length_days: 28,
      },
      blackfly: {
        bite_rate_per_person_per_year: 1000,
        human_blood_index: 0.63,
        use_density_dependence: true,
        manual_density_dependence_params: false,
        k0: 0.0054,
        k1: 0.1459,
        x1: -1,
        hbi_lb: 0.63,
        delta_h_zero: 0.1864987,
        delta_h_inf: 0.002772749,
        c_h: 0.004900419,
      },
      worms: {
        worm_age_stages: 21,
        q_w: 1.0,
        max_worm_age: 21,
      },
      microfilaria: {
        slope_kmf: 0.0478,
        initial_kmf: 0.313,
        mf_age_stages: 21,
        q_m: 0.125,
        max_mf_age: 2.5,
        mf_move_rate: 8.13333,
        use_kmf_const: false,
        kmf_const: 15,
      },
      exposure: {
        Q: 1.2,
        use_onchosim_exposure: false,
      },
      human: {
        min_skinsnip_age: 5,
        max_human_age: 80,
        mean_human_age: 50,
        skin_snip_weight: 2,
        skin_snip_number: 2,
        gender_ratio: 0.5,
        prop_serorevert_fast: 0.5,
      },
      treatments: [],
      vector_control: [],
    },
    output_params: []
  });

  const [config, setConfig] = useState(getDefaultConfig());
  const [lastSavedConfig, setLastSavedConfig] = useState(JSON.stringify(config, null, 2));
  const [isModelRunning, setIsModelRunning] = useState(false);

  const getDefaultOutput = () => ({
    end_time_years: config.runtime_params.total_years,
    start_time_years: 0,
    interval_years: 0.125,
    start_age: 5,
    end_age: 81,
    year_label_start: 1900,
    anti_ov16_test_sens: 0.80,
    anti_ov16_test_spec: 0.99,
    outputs_to_track: [
        "mf_intensity", "mf_prevalence", "population_size", "true_ov16_seroprevalence", 
        "adjusted_ov16_seroprevalence", "l3_per_blackfly", "l3_prevalence_blackflies"
    ]
  });

  const [newOutput, setNewOutput] = useState(getDefaultOutput());

  const getDefaultTreatment = () => ({
    name: '',
    start_year: 50,
    end_year: 100,
    interval_years: 1.0,
    min_age: 5,
    coverage: 0.65,
    rho: 0.3,
    proportion_never_treated: 0,
    use_infection: false,
    infection_threshold: 0,
    drugs: [{ name: 'IVM', allocation_proportion: 1.0 }],
  });

  const [newTreatment, setNewTreatment] = useState(getDefaultTreatment());

  const getDefaultVectorControl = () => ({
    name: '',
    start_year: 50,
    end_year: 60,
    interval_years: 1,
    efficacies: [0.5],
    bounce_back_interval: 1.0,
  });

  const [newVectorControl, setNewVectorControl] = useState(getDefaultVectorControl());

  const handleToggleAdvanceSettings = (field, value) => {
    setShowAdvancedSettings(prev => ({
      ...prev,
      [field]: value
    }));
    document.querySelectorAll(".advanced-setting").forEach((setting, _) => {
      setting.querySelector('input').disabled = !value;
      if (value) {
        setting.classList.remove("advanced-setting-hide");
      } else {
        setting.classList.add("advanced-setting-hide");
      }
    });
  };

  const handleConfigNameChange = (field, value) => {
    setConfig(prev => ({
      ...prev,
      [field]: value
    }));
  };

  const handleRuntimeChange = (field, value) => {
    setConfig(prev => ({
      ...prev,
      runtime_params: { ...prev.runtime_params, [field]: value }
    }));
  };

  const handleBaseChange = (field, value) => {
    setConfig(prev => ({
      ...prev,
      model_params: {
        ...prev.model_params,
        base: { ...prev.model_params.base, [field]: value }
      }
    }));
  };

  const handleBlackflyChange = (field, value) => {
    setConfig(prev => ({
      ...prev,
      model_params: {
        ...prev.model_params,
        blackfly: { ...prev.model_params.blackfly, [field]: value }
      }
    }));
  };

  const handleWormsChange = (field, value) => {
    setConfig(prev => ({
      ...prev,
      model_params: {
        ...prev.model_params,
        worms: { ...prev.model_params.worms, [field]: value }
      }
    }));
  };

  const handleMicrofilariaChange = (field, value) => {
    setConfig(prev => ({
      ...prev,
      model_params: {
        ...prev.model_params,
        microfilaria: { ...prev.model_params.microfilaria, [field]: value }
      }
    }));
  };

  const handleExposureChange = (field, value) => {
    setConfig(prev => ({
      ...prev,
      model_params: {
        ...prev.model_params,
        exposure: { ...prev.model_params.exposure, [field]: value }
      }
    }));
  };

  const handleHumanChange = (field, value) => {
    setConfig(prev => ({
      ...prev,
      model_params: {
        ...prev.model_params,
        human: { ...prev.model_params.human, [field]: value }
      }
    }));
  };

  const addTreatment = () => {
    if (!newTreatment.name.trim()) {
      alert('Treatment name is required');
      return;
    }
    setConfig(prev => ({
      ...prev,
      model_params: {
        ...prev.model_params,
        treatments: [...prev.model_params.treatments, { ...newTreatment }]
      }
    }));
    setNewTreatment(getDefaultTreatment());
  };

  const removeTreatment = (index) => {
    setConfig(prev => ({
      ...prev,
      model_params: {
        ...prev.model_params,
        treatments: prev.model_params.treatments.filter((_, i) => i !== index)
      }
    }));
  };

  const addVectorControl = () => {
    if (!newVectorControl.name.trim()) {
      alert('Vector control name is required');
      return;
    }
    setConfig(prev => ({
      ...prev,
      model_params: {
        ...prev.model_params,
        vector_control: [...prev.model_params.vector_control, { ...newVectorControl }]
      }
    }));
    setNewVectorControl(getDefaultVectorControl());
  };

  const removeVectorControl = (index) => {
    setConfig(prev => ({
      ...prev,
      model_params: {
        ...prev.model_params,
        vector_control: prev.model_params.vector_control.filter((_, i) => i !== index)
      }
    }));
  };

  const addOutput = () => {
    setConfig(prev => ({
      ...prev,
      output_params: [
        ...prev.output_params, { ...newOutput }
      ]
    }));
    setNewOutput(getDefaultOutput());
  };

  const removeOutput = (index) => {
    setConfig(prev => ({
      ...prev,
      output_params: prev.output_params.filter((_, i) => i !== index)
    }));
  };

  const generateJSON = () => {
    if (config.output_params.length === 0) {
      if (!window.confirm("No model output parameters defined. Do you want to continue?")) {
        return null;
      }
    }
    const json = JSON.stringify(config, null, 2);
    onConfigGenerated?.(json);
    setLastSavedConfig(json);
    return json;
  };

  const downloadJSON = () => {
    const json = generateJSON();
    if (!json) {
      return;
    }
    const blob = new Blob([json], { type: 'application/json' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = 'config.json';
    a.click();
    URL.revokeObjectURL(url);
  };

  const resetConfig = () => {
    setConfig(getDefaultConfig());
    setNewTreatment(getDefaultTreatment);
    setNewVectorControl(getDefaultVectorControl());
    setNewOutput(getDefaultOutput);
  }

  const [modelOutput, setModelOutput] = useState('');
  const runModel = async () => {
    setModelOutput('Running...');
    try {
      setIsModelRunning(true);
      const response = await fetch('http://localhost:3001/api/run-model', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify(config)
      });
      const result = await response.json();
      setModelOutput(result.output);
      if (result.success && result.output.includes('Successfully merged to:')) {
        onModelResultsGenerated?.(result);
      }
    } catch (error) {
      setModelOutput('Error: ' + error.message);
    }
    setIsModelRunning(false);
  };

  const allowModelRun = () => {
    return (lastSavedConfig === JSON.stringify(config, null, 2) && !isModelRunning);
  };

  return (
    <div className="config-form">
      <h1>Epioncho IBM Configuration</h1>

      <section className="section">
        <div className="form-group">
          <label>
            <input
              type="checkbox"
              checked={showAdvancedSettings.show_advanced_settings}
              onChange={(e) => handleToggleAdvanceSettings('show_advanced_settings', e.target.checked)}
            />
            Show Advanced Settings
          </label>
        </div>
      </section>

      <section className="section">
        <div className="form-group">
          <label>Configuration Name</label>
            <input
              type="text"
              value={config.name}
              onChange={(e) => handleConfigNameChange('name', e.target.value)}
              placeholder="Config"
            />
        </div>
      </section>

      <section className="section">
        <h2>Model Run Parameters</h2>
        <div className="form-group">
          <label>Total Years</label>
          <input
            type="number"
            value={config.runtime_params.total_years}
            onChange={(e) => handleRuntimeChange('total_years', parseInt(e.target.value))}
          />
        </div>
        <div className="form-group">
          <label>Number of Computer Cores</label>
          <input
            type="number"
            value={config.runtime_params.num_cores}
            onChange={(e) => handleRuntimeChange('num_cores', parseInt(e.target.value))}
          />
        </div>
        <div className="form-group">
          <label>Number of Model Repeats</label>
          <input
            type="number"
            value={config.runtime_params.num_repeats}
            onChange={(e) => handleRuntimeChange('num_repeats', parseInt(e.target.value))}
          />
        </div>
      </section>

      <section className="section">
        <h2>Base Parameters</h2>
        <div className="form-group">
          <label>Seed</label>
          <input
            type="number"
            value={config.model_params.base.seed}
            onChange={(e) => handleBaseChange('seed', parseInt(e.target.value))}
          />
        </div>
        <div className="form-group">
          <label>Population Size</label>
          <input
            type="number"
            value={config.model_params.base.n_people}
            onChange={(e) => handleBaseChange('n_people', parseInt(e.target.value))}
          />
        </div>
        <div className="form-group">
          <label>k_E (Exposure Heterogeneity)</label>
          <select value={config.model_params.base.k_E} onChange={(e) => handleBaseChange('k_E', parseFloat(e.target.value))}>
            <option value={0.2}>0.2</option>
            <option value={0.3}>0.3</option>
            <option value={0.4}>0.4</option>
            <option value={3.5}>3.5</option>
          </select>
        </div>
        <div className="form-group">
          <label>Delta Time (days)</label>
          <input
            type="number"
            step="0.5"
            value={config.model_params.base.delta_time_days}
            onChange={(e) => handleBaseChange('delta_time_days', parseFloat(e.target.value))}
          />
        </div>
      </section>

      <section className="section">
        <h2>Blackfly Parameters</h2>
        <div className="form-group">
          <label>Annual Biting Rate (ABR)</label>
          <input
            type="number"
            value={config.model_params.blackfly.bite_rate_per_person_per_year}
            onChange={(e) => handleBlackflyChange('bite_rate_per_person_per_year', parseFloat(e.target.value))}
          />
        </div>
        <div className="form-group">
          <label>
            <input
              type="checkbox"
              checked={config.model_params.blackfly.use_density_dependence}
              onChange={(e) => handleBlackflyChange('use_density_dependence', e.target.checked)}
            />
            Use Density Dependence
          </label>
        </div>
        <div className="form-group advanced-setting advanced-setting-hide">
          <label>k0 (intercept for relation between blackfly infection intensity and prevalence)</label>
          <input
            type="number"
            step="0.001"
            value={config.model_params.blackfly.k0}
            onChange={(e) => handleBlackflyChange('k0', parseFloat(e.target.value))}
          />
        </div>
        <div className="form-group advanced-setting advanced-setting-hide">
          <label>k1 (slope for relation between blackfly infection intensity and prevalence)</label>
          <input
            type="number"
            step="0.01"
            value={config.model_params.blackfly.k1}
            onChange={(e) => handleBlackflyChange('k1', parseFloat(e.target.value))}
          />
        </div>
        <div className="form-group advanced-setting advanced-setting-hide">
          <label>Human Blood Index</label>
          <input
            type="number"
            step="0.1"
            value={config.model_params.blackfly.human_blood_index}
            onChange={(e) => handleBlackflyChange('human_blood_index', parseFloat(e.target.value))}
          />
        </div>
        <div className="form-group advanced-setting advanced-setting-hide">
          <label>x1 (decrease in HBI in relation to ABR; set to -1 for no change)</label>
          <input
            type="number"
            step="0.00008"
            value={config.model_params.blackfly.x1}
            onChange={(e) => handleBlackflyChange('x1', parseFloat(e.target.value))}
          />
        </div>
        <div className="form-group advanced-setting advanced-setting-hide">
          <label>Lower-bound for HBI (only set if x1 is set)</label>
          <input
            type="number"
            step="0.01"
            value={config.model_params.blackfly.hbi_lb}
            onChange={(e) => handleBlackflyChange('hbi_lb', parseFloat(e.target.value))}
          />
        </div>
      </section>

      <section className="section">
        <h2>Microfilaria Parameters</h2>
        <div className="form-group">
          <label>
            <input
              type="checkbox"
              checked={config.model_params.microfilaria.use_kmf_const}
              onChange={(e) => handleMicrofilariaChange('use_kmf_const', e.target.checked)}
            />
            Use Constant kmf
          </label>
        </div>
        <div className="form-group advanced-setting advanced-setting-hide">
          <label>Rate of MF moving between compartments</label>
          <input
            type="number"
            step="0.1"
            value={config.model_params.microfilaria.mf_move_rate}
            onChange={(e) => handleMicrofilariaChange('mf_move_rate', parseFloat(e.target.value))}
          />
        </div>
        <div className="form-group advanced-setting advanced-setting-hide">
          <label>Number of age-stages for MF</label>
          <input
            type="number"
            step="0.1"
            value={config.model_params.microfilaria.mf_age_stages}
            onChange={(e) => handleMicrofilariaChange('mf_age_stages', parseInt(e.target.value))}
          />
        </div>
      </section>

      <section className="section">
        <h2>Exposure Parameters</h2>
        <div className="form-group advanced-setting advanced-setting-hide">
          <label>Relative Male to Female ratio of exposure to Blackfly Bites</label>
          <input
            type="number"
            step="0.1"
            value={config.model_params.exposure.Q}
            onChange={(e) => handleExposureChange('Q', parseFloat(e.target.value))}
          />
        </div>
        <div className="form-group">
          <label>
            <input
              type="checkbox"
              checked={config.model_params.exposure.use_onchosim_exposure}
              onChange={(e) => handleExposureChange('use_onchosim_exposure', e.target.checked)}
            />
            Use Onchosim Exposure
          </label>
        </div>
      </section>

      <section className="section">
        <h2>Human Parameters</h2>
        <div className="form-group">
          <label>Minimum skin-snip age</label>
          <input
            type="number"
            value={config.model_params.human.min_skinsnip_age}
            onChange={(e) => handleHumanChange('min_skinsnip_age', parseInt(e.target.value))}
          />
        </div>
        <div className="form-group advanced-setting advanced-setting-hide">
          <label>Weight of skin-snip (milligrams)</label>
          <input
            type="number"
            value={config.model_params.human.skin_snip_weight}
            onChange={(e) => handleHumanChange('skin_snip_weight', parseFloat(e.target.value))}
          />
        </div>
        <div className="form-group">
          <label>Number of skin-snips</label>
          <input
            type="number"
            value={config.model_params.human.skin_snip_number}
            onChange={(e) => handleHumanChange('skin_snip_number', parseInt(e.target.value))}
          />
        </div>
        <div className="form-group">
          <label>Proportion of Individuals who Serorevert Fast</label>
          <input
            type="number"
            step="0.01"
            min="0"
            max="1"
            value={config.model_params.human.prop_serorevert_fast}
            onChange={(e) => handleHumanChange('prop_serorevert_fast', parseFloat(e.target.value))}
          />
        </div>
      </section>

      <section className="section">
        <h2>Treatments (MDA)</h2>
        <div className="subsection">
          <h3>Add Treatment</h3>
          <div className="form-group">
            <label>Treatment Name</label>
            <input
              type="text"
              value={newTreatment.name}
              onChange={(e) => setNewTreatment({ ...newTreatment, name: e.target.value })}
              placeholder="e.g., IVM Annual"
            />
          </div>
          <div className="form-row">
            <div className="form-group">
              <label>Start Year</label>
              <input
                type="number"
                value={newTreatment.start_year}
                onChange={(e) => setNewTreatment({ ...newTreatment, start_year: parseInt(e.target.value) })}
              />
            </div>
            <div className="form-group">
              <label>End Year</label>
              <input
                type="number"
                value={newTreatment.end_year}
                onChange={(e) => setNewTreatment({ ...newTreatment, end_year: parseInt(e.target.value) })}
              />
            </div>
            <div className="form-group">
              <label>Interval (years)</label>
              <input
                type="number"
                step="0.5"
                value={newTreatment.interval_years}
                onChange={(e) => setNewTreatment({ ...newTreatment, interval_years: parseFloat(e.target.value) })}
              />
            </div>
          </div>
          <div className="form-row">
            <div className="form-group">
              <label>Coverage</label>
              <input
                type="number"
                step="0.01"
                min="0"
                max="1"
                value={newTreatment.coverage}
                onChange={(e) => setNewTreatment({ ...newTreatment, coverage: parseFloat(e.target.value) })}
              />
            </div>
            <div className="form-group">
              <label>Rho</label>
              <input
                type="number"
                step="0.01"
                min="0"
                max="1"
                value={newTreatment.rho}
                onChange={(e) => setNewTreatment({ ...newTreatment, rho: parseFloat(e.target.value) })}
              />
            </div>
          </div>

          <div className="form-group">
            <label>Drugs</label>
            <div className="drugs-section">
              {newTreatment.drugs.map((drug, idx) => (
                <div key={idx} className="drug-item">
                  <select 
                    value={drug.name}
                    onChange={(e) => {
                      const updatedDrugs = [...newTreatment.drugs];
                      updatedDrugs[idx].name = e.target.value;
                      setNewTreatment({ ...newTreatment, drugs: updatedDrugs });
                    }}
                  >
                    <option value="IVM">IVM</option>
                    <option value="MOX">MOX</option>
                  </select>
                  <input
                    type="number"
                    step="0.1"
                    min="0"
                    max="1"
                    placeholder="Allocation"
                    value={drug.allocation_proportion}
                    onChange={(e) => {
                      const updatedDrugs = [...newTreatment.drugs];
                      updatedDrugs[idx].allocation_proportion = parseFloat(e.target.value);
                      setNewTreatment({ ...newTreatment, drugs: updatedDrugs });
                    }}
                  />
                  <button onClick={() => {
                    const updatedDrugs = newTreatment.drugs.filter((_, i) => i !== idx);
                    setNewTreatment({ ...newTreatment, drugs: updatedDrugs });
                  }} className="btn-remove-small">Remove</button>
                </div>
              ))}
              <button onClick={() => {
                setNewTreatment({
                  ...newTreatment,
                  drugs: [...newTreatment.drugs, { name: 'IVM', allocation_proportion: 0.5 }]
                });
              }} className="btn-add-drug">Add Drug</button>
            </div>
          </div>
          <button onClick={addTreatment} className="btn-add">Add Treatment</button>
        </div>

        <div className="treatment-list">
          {config.model_params.treatments.map((treatment, idx) => (
            <div key={idx} className="treatment-item">
              <h4>{treatment.name}</h4>
              <p>Years {treatment.start_year}-{treatment.end_year}, Coverage: {treatment.coverage}</p>
              <button onClick={() => removeTreatment(idx)} className="btn-remove">Remove</button>
            </div>
          ))}
        </div>
      </section>

      <section className="section">
        <h2>Vector Control</h2>
        <div className="subsection">
          <h3>Add Vector Control</h3>
          <div className="form-group">
            <label>Control Name</label>
            <input
              type="text"
              value={newVectorControl.name}
              onChange={(e) => setNewVectorControl({ ...newVectorControl, name: e.target.value })}
              placeholder="e.g., Larviciding"
            />
          </div>
          <div className="form-row">
            <div className="form-group">
              <label>Start Year</label>
              <input
                type="number"
                value={newVectorControl.start_year}
                onChange={(e) => setNewVectorControl({ ...newVectorControl, start_year: parseInt(e.target.value) })}
              />
            </div>
            <div className="form-group">
              <label>End Year</label>
              <input
                type="number"
                value={newVectorControl.end_year}
                onChange={(e) => setNewVectorControl({ ...newVectorControl, end_year: parseInt(e.target.value) })}
              />
            </div>
            <div className="form-group">
              <label>Interval (years)</label>
              <input
                type="number"
                step="0.5"
                value={newVectorControl.interval_years}
                onChange={(e) => setNewVectorControl({ ...newVectorControl, interval_years: parseFloat(e.target.value) })}
              />
            </div>
          </div>
          <div className="form-group">
            <label>Efficacy</label>
            <input
              type="number"
              step="0.01"
              min="0"
              max="1"
              value={newVectorControl.efficacies[0]}
              onChange={(e) => setNewVectorControl({ ...newVectorControl, efficacies: [parseFloat(e.target.value)] })}
            />
          </div>
          <div className="form-group">
            <label>Years before ABR bounces back to pre-vector control value</label>
            <input
              type="number"
              step="1"
              min="0"
              value={newVectorControl.bounce_back_interval}
              onChange={(e) => setNewVectorControl({ ...newVectorControl, bounce_back_interval: parseInt(e.target.value) })}
            />
          </div>
          <button onClick={addVectorControl} className="btn-add">Add Vector Control</button>
        </div>

        <div className="vector-control-list">
          {config.model_params.vector_control.map((vc, idx) => (
            <div key={idx} className="vc-item">
              <h4>{vc.name}</h4>
              <p>Years {vc.start_year}-{vc.end_year}, Efficacy: {vc.efficacies[0]}</p>
              <button onClick={() => removeVectorControl(idx)} className="btn-remove">Remove</button>
            </div>
          ))}
        </div>
      </section>

      <section className="section">
        <h2>Model Output</h2>
        <div className="subsection">
          <h3>Add Model Output</h3>
          <div className="form-row">
            <div className="form-group">
              <label>Model Output Start Year</label>
              <input
                type="number"
                value={newOutput.start_time_years}
                min="0"
                onChange={(e) => setNewOutput({ ...newOutput, start_time_years: parseInt(e.target.value) })}
              />
            </div>
            <div className="form-group">
              <label>Model Output End Year</label>
              <input
                type="number"
                value={config.runtime_params.total_years}
                min="1"
                onChange={(e) => setNewOutput({ ...newOutput, end_time_years: parseInt(e.target.value) })}
              />
            </div>
            <div className="form-group">
              <label>Output Interval (years)</label>
              <input
                type="number"
                step="0.125"
                value={newOutput.interval_years}
                onChange={(e) => setNewOutput({ ...newOutput, interval_years: parseFloat(e.target.value) })}
              />
            </div>
          </div>
          <div className="form-row">
            <div className="form-group">
              <label>Model Output Start Age (inclusive)</label>
              <input
                type="number"
                value={newOutput.start_age}
                min="0"
                onChange={(e) => setNewOutput({ ...newOutput, start_age: parseInt(e.target.value) })}
              />
            </div>
            <div className="form-group">
              <label>Model Output End Age (exclusive)</label>
              <input
                type="number"
                min="1"
                value={newOutput.end_age}
                onChange={(e) => setNewOutput({ ...newOutput, end_age: parseInt(e.target.value) })}
              />
            </div>
          </div>
          <div className="form-group">
            <label>Year Label (starts at time = 0)</label>
            <input
              type="number"
              step="10"
              min="0"
              value={newOutput.year_label_start}
              onChange={(e) => setNewOutput({ ...newOutput, year_label_start: parseInt(e.target.value) })}
            />
          </div>
          <div className="form-row">
            <div className="form-group">
              <label>anti-Ov16 Test Sensitivity (%)</label>
              <input
                type="number"
                step="10"
                min="0"
                max="100"
                value={newOutput.anti_ov16_test_sens * 100}
                onChange={(e) => setNewOutput({ ...newOutput, anti_ov16_test_sens: parseFloat(e.target.value) / 100 })}
              />
            </div>
            <div className="form-group">
              <label>anti-Ov16 Test Specificity (%)</label>
              <input
                type="number"
                step="10"
                min="0"
                max="100"
                value={newOutput.anti_ov16_test_spec * 100}
                onChange={(e) => setNewOutput({ ...newOutput, anti_ov16_test_spec: parseFloat(e.target.value) / 100 })}
              />
            </div>
          </div>
          <div className="form-group">
            <label>Output Metrics to Track. Ctrl+click (Windows) or  CMD+click (Mac) to select multiple</label>
            <select 
              name="model-outputs"
              multiple
              defaultValue={[
                "mf_intensity", "mf_prevalence", "population_size", 
                "true_ov16_seroprevalence", "adjusted_ov16_seroprevalence", 
                "l3_per_blackfly", "l3_prevalence_blackflies"
              ]}
              onChange={(e) => setNewOutput({ ...newOutput, outputs_to_track: Array.from(e.target.selectedOptions).map(option => option.value) })}
            >
              <option value="mf_prevalence">MF Prevalence</option>
              <option value="mf_intensity">MF Intensity</option>
              <option value="population_size">Population Size</option>
              <option value="true_ov16_seroprevalence">Un-adjusted Ov16 Seroprevalence</option>
              <option value="adjusted_ov16_seroprevalence">Adjusted Ov16 Seroprevalence</option>
              <option value="l3_per_blackfly">L3 Intensity per Blackfly</option>
              <option value="l3_prevalence_blackflies">L3 Prevalence in Blackflies</option>
              <option value="mean_worm_load">Mean Worm Load</option>
              <option value="mean_male_worm_load">Mean Male Worm Load</option>
              <option value="mean_female_worm_load">Mean Female Worm Load</option>
              <option value="mean_fertile_female_worm_load">Mean Fertile Female Worm Load</option>
              <option value="mean_infertile_female_worm_load">Mean Infertile Female Worm Load</option>
              <option value="mean_perm_sterile_female_worm_load">Mean Permanently Sterile Female Worm Load</option>
              <option value="never_treated_percent">Compliance Percent</option>
              <option value="severe_itch_prevalence">Severe Itch Prevalence</option>
              <option value="rsd_prevalence">RSD Prevalence</option>
              <option value="atrophy_prevalence">Atrophy Prevalence</option>
              <option value="hanging_groin_prevalence">Hanging Groin Prevalence</option>
              <option value="depigmentation_prevalence">Depigmentation Prevalence</option>
              <option value="blindness_prevalence">Blindness Prevalence</option>
              <option value="visual_impairment_prevalence">Visual Impairment Prevalence</option>
              <option value="oae_prevalence">OAE Prevalence</option>
            </select>
          </div>

          <button onClick={addOutput} className="btn-add">Add Model Output</button>
        </div>

        <div className="model-output-list">
          {config.output_params.map((mo, idx) => (
            <div key={idx} className="mo-item">
              <h4>Age Group {mo.start_age}-{mo.end_age}</h4>
              <p>Years {mo.start_time_years}-{mo.end_time_years}</p>
              <button onClick={() => removeOutput(idx)} className="btn-remove">Remove</button>
            </div>
          ))}
        </div>
      </section>

      <section className="section actions">
        <section className="section actions json-attr">
          <section className="section actions">
            <button onClick={downloadJSON} className="btn-primary">Download JSON</button>
            <button onClick={generateJSON} className="btn-primary">Generate JSON</button>
          </section>
          <button onClick={resetConfig} className="btn-reset">Reset to Defaults</button>
        </section>
        <button onClick={runModel} id="run-model" className="btn-primary" disabled={!allowModelRun()}>Run Model</button>
      </section>
      <div className="model-output">
        <h3>Model Output</h3>
        <pre>{modelOutput}</pre>
      </div>
    </div>
  );
};

export default ConfigForm;
