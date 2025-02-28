function vare = murase_fitg(curves)
%--------------------------------------------------------------------------
%   vare = murase_fitg(curves)
%--------------------------------------------------------------------------
%   This function generates model-based images by fitting the estimated AIF
%   and tissue curves to the two-compartment model using a linear least
%   squares method. The two-compartment model is linearized by integrating
%   both sides of the equation in order to do the linear least squares
%   method.
%--------------------------------------------------------------------------
%   Inputs:      
%       - curves: time curves for AIF and tissue [nb,nt]
%           -nb: number of curves (first curve is AIF)
%           -nt: number of perfusion time frames
%--------------------------------------------------------------------------
%   Outputs:
%       - vare: model-based estimate of each tissue curve [nb-1,nt]
%           -nb: number of curves (first curve is AIF)
%           -nt: number of perfusion time frames
%--------------------------------------------------------------------------
%   Reference:
%       [1] Murase K, Efficient method for calculating kinetic parameters
%       using T1-weighted dynamic contrast-enhanced magnetic resonance
%       imaging. Magn Reson Med 2004;51(4):858-862
%--------------------------------------------------------------------------

bldcurve=curves(1,:)';

nRegs=size(curves,1)-1;
tisscurves=curves(2:nRegs+1,:);
tisscurve=tisscurves';
nTimes=length(bldcurve);

max_delay=round(size(curves,2))-1;

vare=zeros(size(curves,2),size(curves,1)-1);

num_par=3;
time=(0:(length(bldcurve)-1))/size(curves,2);

parfor ii=1:nRegs
        
    err=zeros(max_delay+1,1);    
    tissue=tisscurve(:,ii);
    
    for j=-max_delay:1:0
        
        t_delay=j;
        tmpbldcurve=bldcurve;
        tmpbldcurve(1-fix(t_delay):nTimes)=bldcurve(1:nTimes+fix(t_delay));
        tmpbldcurve(1:1-fix(t_delay))=0;
        try
            [~,fit,~] = Quant.muraseFit(time,tmpbldcurve,tissue,num_par);        
        catch
            fit=bldcurve;
        end
        err(j+max_delay+1)=sum((tissue-fit).^2);        
    end
    
    pos=find(err==min(err));
    t_delay=(pos)-(max_delay+1);
    
    tmpbldcurve=bldcurve;
    tmpbldcurve(1-fix(t_delay):nTimes)=bldcurve(1:nTimes+fix(t_delay));
    tmpbldcurve(1:1-fix(t_delay))=0;

    tissue=tisscurve(:,ii);    
    [~,fit,~] = Quant.muraseFit(time,tmpbldcurve,tissue,num_par);
    
    vare(:,ii)=fit;

end

end

