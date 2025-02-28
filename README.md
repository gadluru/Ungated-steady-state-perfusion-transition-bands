# Improving Ungated Steady-State Perfusion using Transition Bands

This code reproduces some of the results from the MRM publication, Mendes JK, Le JV, Arai AE, et al. Magn Reson Med 2025. This code was developed and tested on a RockLinux 8.6 operating system, with an AMD EPYC Milan 7543 32 2.8GHz 256 MB cache, 512 GB RAM and Nvidia A100 gpus. Systems with ~64 GB RAM may encounter issues with insufficient memory that require code changes to reduce memory requirements during certain function calls. 

Run recon_CRIMP_ADMM.m to perform the reconstruction for ungated, free-breathing, cardiac phase-resolved myocardial perfusion MRI using Continuous Radial Interleaved simultaneous Multi-slice acquisition at sPoiled steady-state (CRIMP). This reconstruction code builds upon the CRIMP framework to perform reconstruction for the same acquisition with the addition of transition bands to improve ungated steady-state cardiac perfusion. 

Run processing/quantification_code.m to perform blood flow quantification for the reconstruction generated from recon_CRIMP_ADMM.m. This quantification code builds upon the CRIMP quantification framework to perform quantification of the same acquisition with the addition of transition bands to improve ungated steady-state cardiac perfusion. This code uses the pTV registration toolbox for motion compensation and Bullseye Plot code from the Matlab File Exchange.

Reference:
[1] Mendes JK, Le JV, Arai AE, Ranjan R, DiBella EVR, Adluru G. Improving ungated steady-state cardiac perfusion using transition bands. Magn Reson Med 2025.
[2] 

      Johnathan Le
      E-mail: le.johnv@outlook.com

