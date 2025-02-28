function [Gd_curve,T1_curve] = SI2Gd_dic_real(SI_ratio,SI_ratio_pre,dic_real,SI,PD_SI,Nprecontrast)
%--------------------------------------------------------------------------
%   [Gd_curve,T1_curve] = SI2Gd_dic_real(SI_ratio,SI_ratio_pre,dic_real,SI,PD_SI,Nprecontrast)
%--------------------------------------------------------------------------
%   This function converts the signal intensity AIF to [Gd] using the
%   estimated diciontary
%--------------------------------------------------------------------------
%   Inputs:      
%       - SI_ratio: perfusion signal intensity normalized by proton density [1,nt]
%           - nt: number of perfusion time frames
%       - SI_ratio_pre: prefusion precntrast signal intensity normalized by proton density [scalar]
%       - dic_real: signal dictionary to T1 and [Gd] estimation [1,nT1]
%           - nT1: number of T1 values
%       - SI: raw dictionary signal intensity [1,nT1]
%           - nT1: number of T1 values
%       - PD_SI: raw dictionary proton density signal intensity [1,nT1]
%           - nT1: number of T1 values
%       - Nprecontrast: number of precontrast frames [vector]
%--------------------------------------------------------------------------
%   Outputs:
%       - Gd_curve: estimated [Gd] curve [1,nt]
%           - nt: number of perfusion time frames
%       - T1_curve: estimated T1 curve [1,nt]
%           - nt: number of perfusion time frames
%--------------------------------------------------------------------------

%% Convert Signal Intensity to T1

d = abs(SI_ratio_pre' - dic_real);

[~,T1_pre] = min(d,[],2);

dic_real = SI ./ PD_SI(:,T1_pre);

d = abs(SI_ratio' - dic_real);

[~,T1_curve] = min(d,[],2);

T1_pre = mean(T1_curve(Nprecontrast,:),1);

%% Convert T1 to [Gd]
r = 3.7;
Gd_curve = (1./T1_curve - 1./T1_pre)/r*1000;
Gd_curve = Gd_curve - mean(Gd_curve(Nprecontrast,:));

Gd_curve(Gd_curve>15) = 15;
Gd_curve(Gd_curve<-0.3) = -0.3;

end