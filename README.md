# BEC2-microbes
Biogeochemical modifications to the BEC2 biogeochemical model as part of the Regional Ocean Modelling System.

## Overview
The code included in this repository includes both the physical and biogeochemical fortran files to compile the ROMS-BEC2 coupled physical-biogeochemical model in a 0.1 degree (latitude by longitude) configuration in the South East Tropical Pacific. With this model, we performed a 50-year control simulation after 
making major developments to the code, building on the developments made by McCoy et al. (2023). These developments were the addition of explicit bacterial and archaeal biomasses (functional types). They represented the processes of aerobic and anaerobic heterotrophy, aerobic ammonia and nitrite oxidation, 
and anaerobic ammonium oxidation (anammox). The microbial functional types performing these metabolisms were facultative nitrate-reducing heterotrophs (NAR), facultative nitrite-reducing heterotrophs (NIR), aerobic ammonia oxidising archaea (AOA), aerobic nitrite oxidising bacteria (NOB) and anammox bacteria (AOX).
With these microbes added to the model, we performed the 50-year control simulation from initial conditions. Additional idealised experiments that removed NOB via different methods were undertaken in parallel, branched at year 30. 

Output for these simulations are available at https://doi.org/10.5281/zenodo.7653737![image](https://github.com/pearseb/BEC2-microbes/assets/4089112/7d9bcffe-4e04-4f3a-a128-1fc8979b019a)

Descriptions of the BEC biogeochemical model can be found in:

Moore, J. K., Doney, S. C. & Lindsay, K. Upper ocean ecosystem dynamics and iron cycling in a global three-dimensional model. Global Biogeochem Cycles 18, 1–21 (2004).![image](https://github.com/pearseb/BEC2-microbes/assets/4089112/679e7594-e232-4949-8400-e194f8d139e7)
McCoy, D., Damien, P., Clements, D., Yang, S. & Bianchi, D. Pathways of Nitrous Oxide Production in the Eastern Tropical South Pacific Oxygen Minimum Zone. Global Biogeochem Cycles (2023) doi:10.1029/2022GB007670.![image](https://github.com/pearseb/BEC2-microbes/assets/4089112/e4ddb728-2761-48c4-bdbe-047f501e6d97)

Descriptions of the ROMS physical setting can be found in:

Shchepetkin, A. F. An adaptive, Courant-number-dependent implicit scheme for vertical advection in oceanic modeling. Ocean Model (Oxf) 91, 38–69 (2015).![image](https://github.com/pearseb/BEC2-microbes/assets/4089112/62ea9a2f-097f-409e-b0dc-194dd2eb34f4)
Shchepetkin, A. F. & McWilliams, J. C. The regional oceanic modeling system (ROMS): a split-explicit, free-surface, topography-following-coordinate oceanic model. Ocean Model (Oxf) 9, 347–404 (2005).
![image](https://github.com/pearseb/BEC2-microbes/assets/4089112/817abf7e-afae-4a4e-bb07-75133a260854)
