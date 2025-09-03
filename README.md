# swot_swh_fxns
Code accompanying analysis of SWOT KaRIn SWH product. This repository contains two main directories:

- **`swh_fxns/matlab/`**  
  Utility functions for working with SWOT SWH KaRIn observations. Handles pre-processing (loading and deconflicting versions, despiking, patching, masking flagged bad-quality observations, calculating and applying empirical correction, etc.) and post-processing (wrapper for co-location via linear interpolation):  
    -- **`load_swot.m`**  
    -- **`correct_swotswh.m`**  
    -- **`process_swot.m`**  
    -- **`colocate_swot.m`**  
![example_processing_20230419T045829](https://github.com/user-attachments/assets/abff04f8-c198-4445-ae9d-22cf9045a30d)
- **`swh_fxns/python/`**  
  Mirrored functions in Python; TBD
- **`fig_code/`**  
  Code for generating the figures in publication.
