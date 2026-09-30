const express = require('express');
const cors = require('cors');
const { spawn } = require('child_process');
const fs = require('fs');

const app = express();

app.use(express.json());
app.use(cors());

app.post('/api/run-model', (req, res) => {
  const config = req.body;
  
  fs.writeFileSync('config.json', JSON.stringify(config, null, 2));
  
  const model = spawn('../c-epioncho-ibm', ['--config', 'config.json', '--output-folder', '../model_output/']);
  let output = '';
  
  model.stdout.on('data', (data) => {
    output += data.toString();
    console.log(data.toString());
  });
  
  model.stderr.on('data', (data) => {
    output += data.toString();
    console.error(data.toString());
  });
  
  model.on('close', (code) => {
    res.json({ 
      success: code === 0, 
      output: output,
      message: code === 0 ? 'Model completed' : 'Model failed'
    });
  });
});

app.get('/api/get-csv', (req, res) => {
  const filePath = req.query.path;
  fs.readFile(filePath, 'utf8', (err, data) => {
    if (err) res.status(500).json({ error: err.message });
    else res.json({ data: data });
  });
});

app.listen(3001, () => {
  console.log('Server running on http://localhost:3001');
});