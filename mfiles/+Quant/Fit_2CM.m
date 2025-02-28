function param = Fit_2CM(AIF,Tissue,Q)
%--------------------------------------------------------------------------
%   param = Fit_2CM(AIF,Tissue,Q)
%--------------------------------------------------------------------------
%   This function estimates the model parameters from the two-compartment
%   model by performing a non-linear least squares curve fitting with the
%   AIF and tissue curves
%--------------------------------------------------------------------------
%   Inputs:      
%       - AIF: arterial input function for blood flow quantification [1,nt]
%           - nt: number of perfusion time frames
%       - Tissue: tissue curves for blood flow quantification [nt,nb]
%           - nt: number of perfusion time frames
%           - nb: number of tissue curves
%       - Q
%           - Tissue
%               -Mask_heart: Mask for blood flow quantification [nx,ny,nset*nsl]
%                   - nx: spatial x-dimension
%                   - ny: spatial y-dimension
%                   - nt: number of perfusion time frames
%                   - nset: number of slice groups
%                   - nsl: number of sms slices (SMS, multiband=3)
%           -param
%               -TD: time spacing between perfusion frames after interpolation [scalar]
%               -TimeStamps: time duration of perfusion after interpolation [1,nt]
%                   - nt: number of perfusion time frames
%--------------------------------------------------------------------------
%   Outputs:
%       - param
%           - estimated parameters from the two-compartment model (ktrans, kep, vb, dt) [4,nb]
%               -nb: number of tissue curves
%--------------------------------------------------------------------------
%   Reference:
%       [1] Tofts PS, Brix G, Buckley DL, Evelhoch JL, Henderson E, Knopp
%       MV, Larsson HB, Lee TY, Mayr NA, Parker GJ, Port RE, Taylor J,
%       Weisskoff RM. Estimating kinetic parameters from dynamic
%       contrast-enhanced T1-weighted MRI of a diffusable tracer:
%       Standardized quantities and symbols. J Magn Reson Imaging
%       1999;10(3):223-232.
%--------------------------------------------------------------------------

para_low = Q.TCM.para_low;
para_up = Q.TCM.para_up;
para_init = Q.TCM.para_init;

TimeStamps = Q.param.TimeStamps;
TD = Q.param.TD;

Tissue = double(Tissue');
t = TimeStamps/60;

in = {double(AIF),double(t),TD};

options = optimoptions('lsqcurvefit','Display','off');

parfor i=1:size(Tissue,1)
    param(:,i) = lsqcurvefit(@Two_component_model,para_init,in,Tissue(i,:),para_low,para_up,options);
end

end

function c = Two_component_model(x,y)

AIF = y{1};
t = y{2};
TD = y{3};

K_trans = x(1);
K_ep = x(2);
vb = x(3);
dt = x(4);

AIF = interp1(t,AIF,t-dt); AIF(isnan(AIF)) = 0;

c = conv(AIF,(TD./60).*K_trans.*exp(-K_ep*t));

c = c(1:length(t)) + vb*AIF(1:length(t));

end


