function [par,fit,cnum] = muraseFit(time,plasma,tissue,num_par)
%--------------------------------------------------------------------------
%   [par,fit,cnum] = muraseFit(time,plasma,tissue,num_par)
%--------------------------------------------------------------------------
%   This function fits the AIF and tissue curves to the linearized two
%   compartment model to estimate the parameters of the two-compartment
%   model. The estimated parameters are then used to produce model-based
%   estimates of the tissue curves.
%--------------------------------------------------------------------------
%   Inputs:      
%       - time: time duration of the AIF and tissue curves [nt,1]
%       - plasma: AIF curve [nt,1]
%       - tissue: tissue curve [nt,1]
%       - num_par: number of model parameters [scalar]
%           -nt: number of perfusion time frames
%--------------------------------------------------------------------------
%   Outputs:
%       - par: estimated two-compartment model parameters [1,num_par]
%       - fit: model-based tissue curve from estimated parameters [nt,1]
%       - cnum: condition number for inversion [scalar]
%           -nt: number of perfusion time frames
%--------------------------------------------------------------------------
%   Reference:
%       [1] Murase K, Efficient method for calculating kinetic parameters
%       using T1-weighted dynamic contrast-enhanced magnetic resonance
%       imaging. Magn Reson Med 2004;51(4):858-862
%--------------------------------------------------------------------------

% default to fitting with k1,k2,Vp
if (nargin < 4); num_par = 3; end

if (num_par < 2 || num_par > 3); error('muraseFit :: number of parameters must be 2 or 3'); end
if (min(size(time)) ~= 1); error('muraseFit :: time and plasma curves must be vectors'); end

if (size(time,1) == 1); time = time'; end
if (size(plasma,1) == 1); plasma = plasma'; end
if (size(tissue,1) == 1); tissue = tissue'; end

if (size(time) ~= size(plasma)); error('muraseFit :: time and plasma curves are inconsistent'); end
if (size(time) ~= size(tissue)); error('muraseFit :: time and tissue curves are inconsistent'); end

A = zeros([length(time) num_par]);

if (num_par == 2)
    A(:,1) = cumtrapz(time,plasma);
    A(:,2) = -cumtrapz(time,tissue);
elseif (num_par == 3)
    A(:,1) = cumtrapz(time,plasma);
    A(:,2) = -cumtrapz(time,tissue);
    A(:,3) = plasma;  
end


warning off all

par = A\tissue;

fit = A*par;

if (num_par == 3)
    par(1) = par(1) - par(2)*par(3);
end

par = par';

if (nargout == 3); cnum = cond(A); end

end


