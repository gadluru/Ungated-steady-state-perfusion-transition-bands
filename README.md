%% Improving Ungated Steady-State Perfusion using Transition Bands
% This code performs the reconstruction from the MRM publication,
% Mendes JK, Le JV, Arai AE, et al. Magn Reson Med 2025. This code
% was developed and tested on a RockLinux 8.6 operating system, with an AMD
% EPYC Milan 7543 32 2.8GHz 256 MB cache, 512 GB RAM and Nvidia A100 gpus.
% Systems with ~64 GB RAM may encounter issues with insufficient memory 
% that require code changes to reduce memory requirements during certain 
% function calls.
%--------------------------------------------------------------------------
%   Reference:
%       [1] Mendes JK, Le JV, Arai AE, Ranjan R, DiBella EVR, Adluru G.
%           Improving ungated steady-state cardiac perfusion using transition
%           bands. Magn Reson Med 2025.
%--------------------------------------------------------------------------
%   Author:
%       Johnathan Le
%       E-mail: le.johnv@outlook.com
%--------------------------------------------------------------------------
