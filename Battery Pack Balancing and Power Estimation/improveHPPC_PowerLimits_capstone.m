addpath readonly
load ./readonly/CellModel.mat

default.z0 = 0.5;
default.T = 25;
default.dT = 10;
default.eta = 1;
default.ns = 1;
default.np = 1;
limits.zMin = 0.1;
limits.zMax = 0.9;
limits.vMin = 2.8;
limits.vMax = 4.3;
limits.iMin = -200;
limits.iMax = 350;
limits.pMin = -1000;
limits.pMax = 1000;
default.limits = limits;
[pChg,pDis] = HPPCpower(default.z0,default.T,default.dT,default.eta,default.ns,default.np,model,default.limits)

function [pChg, pDis] = HPPCpower(z0,T,dT,eta,ns,np,model,limits)
 
    zMin = limits.zMin; zMax = limits.zMax; % Retrieve SOC limits [unitless]
    vMin = limits.vMin; vMax = limits.vMax; % Retrieve voltage limits [V]
    iMin = limits.iMin; iMax = limits.iMax; % Retrieve current limits [A]
    pMin = limits.pMin; pMax = limits.pMax; % Retrieve design power limits [W]

    Q = getParamESC('QParam',T,model); 
    iChgPulse = 10*Q*[zeros(5,1); -ones(dT,1); zeros(5,1)];  % [A] charge pulse
    iDisPulse = 10*Q*[zeros(5,1);  ones(dT,1); zeros(5,1)];  % [A] discharge pulse
    [vk,~,~,~,~] = simCell(iChgPulse,T,model,1,z0,0,0);
    rChg  = abs((max(vk)-vk(1))/min(iChgPulse))
    [vk,~,~,~,~] = simCell(iDisPulse,T,model,1,z0,0,0);
    rDis  = abs((min(vk)-vk(1))/max(iDisPulse))

    % HPPC Power Estimation: Truth
    OCV      = OCVfromSOCtemp(z0,T,model);
    iDisMaxV = (OCV-vMin)/rDis;
    iDisMaxZ = (z0 - zMin)*3600*Q/dT;
    iDisMax  = np*max(0,min([iDisMaxV;iDisMaxZ;iMax*ones(size(z0))]));
    %   OCVDis_2 = OCVfromSOCtemp(z0 - iDisMax*dT/np/(Q*3600),T,model);
    pDisMax  = min(ns*np*vMin*iDisMax,ns*pMax*ones(size(z0)));
  
    iChgMinV = (OCV-vMax)/rChg;
    iChgMinZ = (z0 - zMax)*3600*Q/(dT*eta);
    iChgMin  = np*max([iChgMinV;iChgMinZ;iMin*ones(size(z0))]);
    %   OCVChg_2 = OCVfromSOCtemp(z0 - iChgMin*eta*dT/np/(3600*Q),T,model);
    pChgMin  = min(0,max(ns*np*vMax*iChgMin,ns*pMin*ones(size(z0))));

    pDis = pDisMax;
    pChg = pChgMin; 
end




