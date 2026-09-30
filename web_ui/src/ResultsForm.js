import React, { useState, useMemo } from 'react';
import { LineChart, Line, Area, AreaChart, XAxis, YAxis, CartesianGrid, Tooltip, Legend, ResponsiveContainer, ComposedChart } from 'recharts';
import { quantile } from 'd3';
import './ResultsForm.css';

const ResultsForm = ({ results }) => {
  const [csvData, setCsvData] = useState(null);
  const [selectedMetric, setSelectedMetric] = useState('mf_prevalence');
  const [ageFilter, setAgeFilter] = useState({ start: 5, end: 81 });
  const [uniqueAgeGroups, setUniqueAgeGroups] = useState([]);

  React.useEffect(() => {
    if (results?.output) {
      const match = results.output.match(/Successfully merged to: (.+)/);
      if (match) {
        fetchCSV(match[1]);
      }
    }
  }, [results]);

  const fetchCSV = async (filePath) => {
    try {
      const response = await fetch(`http://localhost:3001/api/get-csv?path=${filePath}`);
      const data = await response.json();
      setCsvData(parseCSV(data.data));
    } catch (error) {
      console.error('Error fetching CSV:', error);
    }
  };

  const parseCSV = (csvText) => {
    const lines = csvText.split('\n');
    const headers = lines[0].split(',');
    const processed_data =  lines.slice(1).map(line => {
      const values = line.split(',');
      const obj = {};
      headers.forEach((header, idx) => {
        obj[header.trim()] = isNaN(values[idx]) ? values[idx] : parseFloat(values[idx]);
      });
      return obj;
    }).filter(obj => obj.output_year);

    const age_groups = [...new Set(processed_data.map(row => JSON.stringify({ start: row.age_start, end: row.age_end })))];
    setUniqueAgeGroups(age_groups.map(a => JSON.parse(a)));
    return processed_data;
  };

  const filteredData = useMemo(() => {
    if (!csvData) return [];
    
    const grouped = csvData
        .filter(row => row.age_start === ageFilter.start && row.age_end === ageFilter.end)
        .reduce((acc, row) => {
        const year = row.output_year;
        if (!acc[year]) acc[year] = [];
        acc[year].push(row[selectedMetric]);
        return acc;
        }, {});

    return Object.entries(grouped).map(([year, values]) => {
            const sorted = values.sort((a, b) => a - b);
            const mean = sorted.reduce((a, b) => a + b) / sorted.length;
            const lb = quantile(sorted, 0.025);//sorted[Math.floor(sorted.length * 0.025)];
            const ub = quantile(sorted, 0.975);//sorted[Math.floor(sorted.length * 0.975)];
            
            return {
            output_year: parseFloat(year),
            [`${selectedMetric}_mean`]: mean,
            [`${selectedMetric}_lb`]: lb,
            [`${selectedMetric}_ub`]: ub
            };
        }).sort((a, b) => a.output_year - b.output_year);
    }, [csvData, selectedMetric, ageFilter]);

  const metrics = [
    'mf_prevalence',
    'mf_intensity',
    'population_size',
    'true_ov16_seroprevalence',
    'adjusted_ov16_seroprevalence',
    'l3_per_blackfly',
    'l3_prevalence_blackflies'
  ];

  return (
    <div className="results-form">
      <h2>Results</h2>

      {!csvData ? (
        <div className="message">
            <p>Run the model to see results</p>
        </div>
        ) : (
        <>
          <div className="controls">
            <div className="control-group">
              <label>Metric</label>
              <select value={selectedMetric} onChange={(e) => setSelectedMetric(e.target.value)}>
                {metrics.map(m => (
                  <option key={m} value={m}>{m}</option>
                ))}
              </select>
            </div>
            <div className="control-group">
              <label>Age Range</label>
              <select value={JSON.stringify(ageFilter)} onChange={(e) => setAgeFilter(JSON.parse(e.target.value))}>
                {uniqueAgeGroups.map((age, idx) => (
                    <option key={idx} value={JSON.stringify(age)}>
                        {age.start} - {age.end}
                    </option>
                ))}
              </select>
            </div>
          </div>

            <ResponsiveContainer width="100%" height={400}>
                <ComposedChart data={filteredData}>
                    <defs>
                        <linearGradient id="colorUv" x1="0" y1="0" x2="0" y2="1">
                            <stop offset="5%" stopColor="#8884d8" stopOpacity={0.3}/>
                            <stop offset="95%" stopColor="#8884d8" stopOpacity={0}/>
                        </linearGradient>
                    </defs>
                    <CartesianGrid strokeDasharray="3 3" />
                    <XAxis dataKey="output_year" />
                    <YAxis domain={[0, 1]} />
                    <Tooltip />
                    <Line 
                        type="monotone" 
                        dataKey={`${selectedMetric}_lb`} 
                        stroke="url(#colorUv)"
                        strokeWidth={2}
                        dot={false}
                        strokeDasharray="5 5"
                        isAnimationActive={false}
                    />
                    <Line 
                        type="monotone" 
                        dataKey={`${selectedMetric}_mean`} 
                        stroke="#8884d8" 
                        strokeWidth={2}
                        fill="none"
                        dot={false}
                        isAnimationActive={false}
                    />
                    <Line 
                        type="monotone" 
                        dataKey={`${selectedMetric}_ub`} 
                        stroke="url(#colorUv)"
                        strokeWidth={2}
                        strokeDasharray="5 5"
                        dot={false}
                        isAnimationActive={false}
                    />
                </ComposedChart>
            </ResponsiveContainer>
        </>
      )}
    </div>
  );
};

export default ResultsForm;