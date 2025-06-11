# PFC_uStim

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

   Filename: reactiontimes.m       
   Purpose:  Computes behavioral stats (RTs)      
        
   Filename: saccadeprecision.m       
   Purpose:  Computes behavioral stats (saccade precision)    
      
   Filename: errors.m       
   Purpose:  Computes behavioral stats (error fraction)    
      
   Filename: NBdecoder_fit.m      
   Purpose:  Fits Naive Bayes decoder
      
   Filename: FRs_stats.m      
   Purpose:  Computes neural activity statistics
      
   Filename: fitFAtouStimdata.m      
   Purpose:  Fits factor analysis model to neural activity
      
   Filename: FAstats.m      
   Purpose:  Computes neural activity statistics in the dominant subspace of the data, defined using FA.
      
   Filename: dPCAstats.m       
   Purpose:  Computes neural activity statistics in the memory subspace of the data, defined using dPCA.


Directory Structure:

PFC_uStim/  
├── README.md  
├── license.md  
└── code/  
&nbsp; &nbsp; &nbsp; &nbsp; ├── addpaths.m  
&nbsp; &nbsp; &nbsp; &nbsp; ├── behavior/  
&nbsp; &nbsp; &nbsp; &nbsp; ├── decoder/  
&nbsp; &nbsp; &nbsp; &nbsp; ├── FRs/  
&nbsp; &nbsp; &nbsp; &nbsp; ├── subspaces/  
&nbsp; &nbsp; &nbsp; &nbsp; └── utils/  


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

Data accessibility: The data package "PFC_uStim_data" can be downloaded from Zenodo 10.5281/zenodo.15640851  

Contact for Issues: Joana Soldado-Magraner, jsoldadomagraner@cmu.edu

