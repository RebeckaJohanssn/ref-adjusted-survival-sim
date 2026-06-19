// master.do
// Master script for running the full simulation study.
// Pipeline: generate truth -> generate data -> analyse data -> summarise results

// Edit ROOT to point to the local path of this repository

global ROOT "PATH_TO_REPOSITORY"
cd "${ROOT}"

// Load global settings
do "${ROOT}/globals.do"

// Step 1: Generate the true survival
do "${ROOT}/Do/generate_truth.do"

// Step 2: Generate simulated datasets (parallelised version)
do "${ROOT}/Do/generate_data_parallel.do"
  
// Non-parallel alternative
// do "${ROOT}/Do/generate_data.do"

// Step 3: Analyse simulated datasets (parallelised version)
do "${ROOT}/Do/analyse_data_parallel.do"
  
// Non-parallel alternative
// do "${ROOT}/Do/analyse_data.do"

// Step 4: Summarise results across simulations
do "${ROOT}/Do/summarise_results.do"
