import React, { useState } from 'react';
import ConfigForm from './ConfigForm';
import ResultsForm from './ResultsForm';
import './App.css';

function App() {
  const [generatedJSON, setGeneratedJSON] = useState(null);

  const handleConfigGenerated = (json) => {
    setGeneratedJSON(json);
  };

  const [modelResults, setModelResults] = useState(null);

  const handleModelResults = (results) => {
    setModelResults(results);
  };

  return (
    <div className="app">
      <script src="https://d3js.org/d3.v4.min.js"></script>
      <div className="container">
        <div className="setup">
          <div className="form-section">
            <ConfigForm 
              onConfigGenerated={handleConfigGenerated}
              onModelResultsGenerated={handleModelResults} 
            />
          </div>
          
          {generatedJSON && (
            <div className="json-section">
              <h2>Generated Configuration</h2>
              <pre>{generatedJSON}</pre>
            </div>
          )}
        </div>
        <div className="results">
          <ResultsForm results={modelResults} />
        </div>
      </div>
    </div>
  );
}

export default App;