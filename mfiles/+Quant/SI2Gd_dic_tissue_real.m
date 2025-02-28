function [Gd_curve,T1_all,T1_pre_all] = SI2Gd_dic_tissue_real(SI_ratio,SI_ratio_pre,dic_real,SI,PD_SI,Nprecontrast)
%--------------------------------------------------------------------------
%   [Gd_curve,T1_all,T1_pre_all] = SI2Gd_dic_tissue_real(SI_ratio,SI_ratio_pre,dic_real,SI,PD_SI,Nprecontrast)
%--------------------------------------------------------------------------
%   This function converts the tissue signal intensity to [Gd] using the
%   estimated diciontary
%--------------------------------------------------------------------------
%   Inputs:      
%       - SI_ratio: perfusion signal intensity normalized by proton density [ny,nt]
%           - ny: spatial y-dimension
%           - nt: number of perfusion time frames
%       - SI_ratio_pre: prefusion precntrast signal intensity normalized by proton density [1,ny]
%           - ny: spatial y-dimension
%       - dic_real: signal dictionary to T1 and [Gd] estimation [1,nT1]
%           - nT1: number of T1 values
%       - SI: raw dictionary signal intensity [1,nT1]
%           - nT1: number of T1 values
%       - PD_SI: raw dictionary proton density signal intensity [1,nT1]
%           - nT1: number of T1 values
%       - Nprecontrast: number of precontrast frames [vector]
%--------------------------------------------------------------------------
%   Outputs:
%       - Gd_curve: estimated [Gd] curve [ny,nt]
%           - ny: spatial y-dimension
%           - nt: number of perfusion time frames
%       - T1_curve: estimated T1 curve [ny,nt]
%           - ny: spatial y-dimension
%           - nt: number of perfusion time frames
%--------------------------------------------------------------------------

[ny,nt] = size(SI_ratio);

SI_ratio_pre = reshape(SI_ratio_pre,[ny,1]);
SI_ratio = reshape(SI_ratio,[ny,1,nt]);

d = abs(SI_ratio_pre - dic_real);
[~,T1_pre_all] = min(d,[],2);

dic_real = SI ./ PD_SI(:,T1_pre_all)';

d = abs(SI_ratio - dic_real);
[~,T1_all] = min(d,[],2);

T1_all = reshape(T1_all,[ny,1,nt]);
T1_pre_all = mean(T1_all(:,:,Nprecontrast),3);

r = 3.7;
Gd_curve = (1./T1_all - 1./T1_pre_all) / r * 1000;
Gd_curve = Gd_curve - mean(Gd_curve(:,:,Nprecontrast),3);

Gd_curve(Gd_curve>15) = 15;
Gd_curve(Gd_curve<-0.3) = -0.3;

end