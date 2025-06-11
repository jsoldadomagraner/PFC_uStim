# PFC_uStim

-------------------
GENERAL INFORMATION
-------------------

Code package for the manuscript: 
"Robustness of working memory to prefrontal cortex microstimulation", Soldado-Magraner et al. 2025 
JNeuroscience, https://doi.org/10.1523/JNEUROSCI.2197-24.2025

Version: v1.0
Date Created: 2025/06/09

Author Information
    Name: Joana Soldado-Magraner
    ORCiD: 0000-0003-3607-7264
    Institution: Carnegie Mellon University
    Email: jsoldadomagraner@cmu.edu
    GitHub Username: jsoldadomagraner

---------------------
FILES & USAGE
---------------------

Main Files:
   A. Filename: reactiontimes.m       
      Purpose:  Computes behavioral stats (RTs)      
      Type: script
        
   B. Filename: saccadeprecision.m       
      Purpose:  Computes behavioral stats (saccade precision)    
      Type: script
      
   C. Filename: errors.m       
      Purpose:  Computes behavioral stats (error fraction)    
      Type: script
      
   D. Filename: NBdecoder_fit.m      
      Purpose:  Fits Naive Bayes decoder
      Type: script
      
   E. Filename: FRs_stats.m      
      Purpose:  Computes neural activity statistics
      Type: script
      
   F. Filename: fitFAtouStimdata.m      
      Purpose:  Fits factor analysis model to neural activity
      Type: script
      
   G. Filename: FAstats.m      
      Purpose:  Computes neural activity statistics in the dominant subspace of the data, defined using FA.
      Type: script
      
   H. Filename: dPCAstats.m       
      Purpose:  Computes neural activity statistics in the memory subspace of the data, defined using dPCA.
      Type: script


Directory Structure:

PFC_uStim/
├── README.md
├── license.md
└── code/
    ├── addpaths.m
    ├── behavior/
    ├── decoder/
    ├── FRs/
    ├── subspaces/
    └── utils/


System Requirements:
   Programming Language and version: Matlab R2019b
   Key Dependencies: The code package consists of customly written functions. It also uses functions to fit Factor Analysis models, written by Byron Yu, Version 1.01 2011, and the dPCA code package written by Kobak et al. 2016. These functions are included in the package (within the utils folder).
   

Usage Instructions:
   How to Run: Run first the main functions listed above to compute statistics based on neural and behavioral data, and then the rest of the functions can be used to recreate the different plots and analyses that appear on the manuscript.


--------------------------
REPRODUCIBILITY & ACCESS
--------------------------

License: CC BY-NC-SA 4.0
Related Publication: https://doi.org/10.1523/JNEUROSCI.2197-24.2025

Contact for Issues: Joana Soldado-Magraner, jsoldadomagraner@cmu.edu

