# Improving Ungated Steady-State Perfusion using Transition Bands

This code reproduces some of the results from the MRM publication, Mendes JK, Le JV, Arai AE, et al. Magn Reson Med 2025. This code was developed and tested on a RockLinux 8.6 operating system, with an AMD EPYC Milan 7543 32 2.8GHz 256 MB cache, 512 GB RAM and Nvidia A100 gpus. Systems with ~64 GB RAM may encounter issues with insufficient memory that require code changes to reduce memory requirements during certain function calls. 

Run recon_CRIMP_ADMM.m to perform the reconstruction for ungated, free-breathing, cardiac phase-resolved myocardial perfusion MRI using Continuous Radial Interleaved simultaneous Multi-slice acquisition at sPoiled steady-state (CRIMP). This reconstruction code builds upon the CRIMP framework to perform reconstruction for the same acquisition with the addition of transition bands to improve ungated steady-state cardiac perfusion. 

Run processing/quantification_code.m to perform blood flow quantification for the reconstruction generated from recon_CRIMP_ADMM.m. This quantification code builds upon the CRIMP quantification framework to perform quantification of the same acquisition with the addition of transition bands to improve ungated steady-state cardiac perfusion. This code uses the pTV registration toolbox for motion compensation and Bullseye Plot code from the Matlab File Exchange.

<br />
<br />

![Bullseye](https://github.com/user-attachments/assets/572854ca-0445-4f4e-94f8-36b925992923)
Figure 7. Resting myocardial blood flow results (reported as Ktrans) for a normal human subject (HR=56 bpm) using transition bands.

<br />
<br />

![example_human_dataset](https://github.com/user-attachments/assets/30c5bb3d-c494-4ff7-86cb-96d8555afa41)
Supporting Information Figure 3. Movie demonstrating the proposed ungated steady-state perfusion reconstruction used to generate myocardial blood flow maps (reported as Ktranrs) of the normal human subject shown in Figure 7.

<br />
<br />

Reference:

[1] Mendes JK, Le JV, Arai AE, Ranjan R, DiBella EVR, Adluru G. Improving ungated steady-state cardiac perfusion using transition bands. Magn Reson Med 2025.

[2] Tian Y, Mendes J, Wilson B, Ross A, Ranjan R, DiBella E, Adluru G. Whole-heart, ungated, free-breathing, cardiac-phase-resolved myocardial perfusion MRI by using Continuous Radial Interleaved simultaneous Multi-slice acquisitions at sPoiled steady-state (CRIMP). Magn Reson Med 2020;84(6):3071-3087

[3] Vishnevskiy V, Gass T, Szekely G, Tanner C, Goksel O. Isotropic total variation regularization of displacements in parametric image registration. IEEE Trans Med Imaging 2017;36(2):385-395

[4] Adrian (2025). Bullseye Plot.zip (https://www.mathworks.com/matlabcentral/fileexchange/47454-bullseye-plot-zip), MATLAB Central File Exchange. Retrieved February 26, 2025. 

<br />
<br />

Contact:

Johnathan Le

le.johnv@outlook.com

Jason Mendes

Jason.Menes@hsc.utah.edu

Edward Dibella

Edward.DiBella@hsc.utah.edu

Ganesh Adluru

Ganeshsharma.Adluru@hsc.utah.edu
