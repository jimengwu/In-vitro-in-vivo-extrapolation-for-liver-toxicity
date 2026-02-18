# In vitro–in vivo extrapolation (IVIVE) for liver toxicity

## Code and models

* **`helper_function.R`** : shared helper functions used across scripts.
* **Stan models**:

  * **`hill_single.stan`** : Hill model (3-parameter form)
  * **`hill_bc_single.stan`**: Hill model (4-parameter form)

## Analysis workflow (script numbering)

The numbered scripts reflect the analysis sequence:

1. **`0.read_process_data.R`**: preprocesses compiled toxicity data (e.g., harmonises concentration units, substance names, assay names).
2. **Dose–response curve fitting** is split by endpoint class:

   * **`1.1...`** cytotoxicity
   * **`1.2...`** metabolic cell stress
   * **`1.3...`** cytokine
3. **Post-modelling filtering and screening** of fitted BMD results:

   * Outputs are saved under `results/` in separate subfolders for **human** and **mouse**.
4. **IVIVE translation**

   * After in vitro BMDs are extracted and filtered, results are then translated to corresponding **in vivo equivalent** doses using the different PBPK models for mouse and human system. 
  

### Results Folder structure

* **`mouse_dose_response_curve/`** contains fitted dose–response plots for the mouse dataset, organised by endpoint class:

  * **`cytotoxicity/`**
  * **`cytokine/`**
  * **`metabolic_cell_stress/`**
* Within each endpoint-class folder:

  * The **plot files** are linked to the corresponding **input data + metadata** saved as **`.Rdata`** files via **`CASE_ID`**.
  * Cases excluded after post-modelling filtering are retained in subfolders named by the relevant **screening-criteria step** within each endpoint class.
* **`human_dose_response_curve/`** mirrors the same structure and file format as the mouse dataset.

