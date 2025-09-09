# swot_swh_fxns
Code accompanying analysis of SWOT KaRIn SWH product. This repository contains two main directories:

- **`swh_fxns/matlab/`**  
  Utility functions for working with SWOT SWH KaRIn observations. Handles pre-processing (loading and deconflicting versions, despiking, patching, masking flagged bad-quality observations, calculating and applying empirical correction, etc.) and post-processing (wrapper for co-location via linear interpolation):  
    -- **`load_swot.m`**  
    -- **`correct_swotswh.m`**  
    -- **`process_swot.m`**  
    -- **`colocate_swot.m`**  
![example_processing_20230419T045829](https://github.com/user-attachments/assets/fe0dbfc3-f614-4195-b71d-4a12e8e1a81c)

- **`fig_code/`**  
  Code for generating the figures in publication.
