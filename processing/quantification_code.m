
%--------------------------------------------------------------------------
%   Quantification for ungated, free-breathing, cardiac phase-resolved
%   myocardial perfusion MRI using Continuous Radial Interleaved
%   simultaneous Multi-slice acquisition at sPoiled steady-state (CRIMP). 
%   This quantification code builds upon the CRIMP framework
%   to perform quantification for the same acquisition with the addition of
%   transition bands to improve ungated steady-state cardiac perfusion.
%   This code uses the pTV registration toolbox for motion compensation 
%   and Bullseye Plot code from the Matlab File Exchange.
%--------------------------------------------------------------------------
%   Reference:
%       [1] Mendes JK, Le JV, Arai AE, Ranjan R, DiBella EVR, Adluru G.
%           Improving ungated steady-state cardiac perfusion using transition
%           bands. Magn Reson Med 2025.
%       [2] Tian Y, Mendes J, Wilson B, Ross A, Ranjan R, DiBella E, Adluru G.
%           Whole-heart, ungated, free-breathing, cardiac-phase-resolved 
%           myocardial perfusion MRI by using Continuous Radial Interleaved
%           simultaneous Multi-slice acquisitions at sPoiled steady-state 
%           (CRIMP). Magn Reson Med 2020;84(6):3071-3087
%       [3] Vishnevskiy V, Gass T, Szekely G, Tanner C, Goksel O. Isotropic
%           total variation regularization of displacements in parametric image
%           registration. IEEE Trans Med Imaging 2017;36(2):385-395
%       [4] Adrian (2025). Bullseye Plot.zip (https://www.mathworks.com/
%           matlabcentral/fileexchange/47454-bullseye-plot-zip),
%           MATLAB Central File Exchange. Retrieved February 26, 2025. 
%--------------------------------------------------------------------------
%   Author:
%       Johnathan Le
%       E-mail: le.johnv@outlook.com
%--------------------------------------------------------------------------

addpath(genpath('../mfiles/'))

dirs = dir('../ReconData/example_human_dataset.mat');
Q_init = Quant.init(dirs);
Q_curves = Quant.get_MISS_tissue_2_SMS_3(Q_init);
