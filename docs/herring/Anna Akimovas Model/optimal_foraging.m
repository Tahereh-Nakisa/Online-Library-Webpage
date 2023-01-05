% %   optimal foraging  routine for a herring larve. Input variables are:
% fl - larval length (in mm)
% prey_l_bins - mean lenght of zooplankton within size-spectra bins (in mm)
% prey_dw_bins - mean dry weight of zooplankton within size-spectra bins (in mg)
% prey_c_bins - zooplankton concentration in each bin (ind/m3)

% % Output variables:
% cons - larval foraging rate( in mg/h)
% ops - is the optimal prey size (in mm)
% mposps - maximal possible prey size due to gape limitation (in mm)
% maxps - maximum targeted prey size (maximal prey size included in the
% diet) in mm
% minps - minimun targeted prey size (minimal prey size included in the
% diet) in mm

% % Akimova Anna, Thünen-Institute of Sea Fisheries, 04.01.2023


function [cons, ops, mposps, maxps, minps]=...
    optimal_foraging(fl, prey_l_bins, prey_dw_bins, prey_c_bins)

prey.l=prey_l_bins; 
prey.dw=prey_dw_bins; 
prey.co=prey_c_bins;
larvae.l=fl; 

if sum(isnan(prey.l)) || sum(isnan(prey.dw)) || sum(isnan(prey.co))
    disp('Error in the input data')
    return;
end


 % % here the variables are converted to the unites used in the optimal foraging routine in Hufnagl and Peck 2011
 prey.l=prey.l*1e3; % prey size is in µm in Hufnagl and Peck, 2011
 prey.dw=prey.dw*1e3; % % prey weight is in µg
 prey.co=prey.co*1e-6; % % prey concentration is in ind/mm3
               
 max_p_a= 2200;      %  2200
 max_p_b= -2;         % 2
 CS_b=1.1;   % % capture success parameter
 HT_a=0.264; % % handling time parameter1
 HT_b=20;    % % handling time parameter1
 PF=0.35;    % % pause frequency
 PD=1.3;     % % 1.3  pause duration mean from Myrons data > 9mm  < 13 mm
 w=1.3;      % turbulent velocity [mm/s] from Dower et al 2002
        
        
 prey.v=(prey.l/1000)*3; % % prey speed in mm/s
 larvae.v=127.21/(  1+exp(-(larvae.l-30.566)/3.67) )+6.7; %larval swimming speed [mm/s]
        
 alpha = 0.0167*exp(9.14-2.4*log(larvae.l)+0.229*(log(larvae.l))^2); % angle of visual acuity
 alpha = alpha * pi/180; % % in radians   
 RD= 1e-3*prey.l./(2*tan(alpha/2)); % % reactive distance (in mm)
 V=sqrt(prey.v.^2+2*w.^2+larvae.v.^2); % Velocity component of contact rate mm/s  
 % % actually to make it correct: larvae.l should be also in �m 
 PS_max= max_p_a./(1+(larvae.l/14)^(max_p_b)); % maximum prey size consumed in µm
        
 % % here is the prey.co is converted to ind/ml as in was used in Hufnagl
 % and Peck, 2011
ER= 0.5*pi.*(RD.^2).*(prey.co./1000)*(1-PF*PD).*V;
CS=max(1.1-CS_b*prey.l/PS_max,0);  % capture success in ind/s
HT=exp(HT_a.*10.^(HT_b.*(prey.l*1e-3/larvae.l))); % handling time [s])

% prey type ranking
rank_par=prey.dw.*CS./HT; 
if ~isempty(rank_par(rank_par<0))
    disp ('Alarm ranking!')
    return;
end

bb=rank_par(rank_par>=0); % ranking all prey classes with length<PS_max
[ranking ii]=sort(bb, 'descend');
% % cost-benefits ratio
par_bc=cumsum(prey.dw(ii).*CS(ii).*ER(ii))./(1+cumsum(ER(ii).*HT(ii)));
benefit=cumsum(prey.dw(ii).*CS(ii).*ER(ii));
cost=1+cumsum(ER(ii).*HT(ii));
[mm ind_mm]= max(par_bc); clear mm
cons=max(par_bc)*60*60; % % in µg per hour
cons=cons*1e-3; % % convert to mg/h
        
prey.dw_ranked=prey.dw(ii);
prey.l_ranked=prey.l(ii);
ops=prey.l_ranked(1)*1e-3; % % optimal prey size in mm   
mposps=PS_max*1e-3; % % max prey length in mm
maxps=max(prey.l_ranked(1:ind_mm))*1e-3;
minps=min(prey.l_ranked(1:ind_mm))*1e-3;


