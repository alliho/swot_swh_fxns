# swot_swh_fxns
Code accompanying analysis of SWOT KaRIn SWH product. This repository contains two main directories:

- **`swh_fxns/`**  
  Utility functions for working with SWOT SWH KaRIn observations. Handles pre-processing (loading and deconflicting versions, despiking, patching, masking flagged bad-quality observations, calculating and applying empirical correction, etc.) and post-processing (wrapper for co-location via linear interpolation):  
    -- **`load_swot.m`**  
    -- **`correct_swotswh.m`**  
    -- **`process_swot.m`**  
    -- **`colocate_swot.m`**  


- **`fig_code/`**  
  Code for generating the figures in publication.
