% % function describing hourly energy budget of a herring larvae
% % input parameters are:
% length - larval length at the previous time step (in mm)
% weight - larval weight at the previous time step (in mug)
% light - light (0 during the night or 1 during the day)
% for_consump - foraging rate per hour(mug/h)
% t - temperature
% % gc_bfore - gut content at the previous time step

% % output parameters
% % energy_budget  - energy budget (in mug/h)
% % egain - energy gain
% % eloss - energy loss
% % gc_after - gut content at the end of the time step

% % Akimova Anna, Thünen-Institute of Sea Fisheries, 04.01.2023



function [energy_budget, egain, eloss, gc_after]=...
    energy_budget(length, weight, light, for_consump, gc_before, t)


if (length<0) || (isnan(length))
    disp('Error: length is negative of NaN');
    return;
end

if (weight<0) || (isnan(weight))
    disp('Error: weight is negative or NaN');
    return;
end

if (gc_before<0) || (isnan(gc_before))
    disp('Error: gut content is negative of NaN');
    return;
end

if (t<0) || (isnan(t))
    disp('Error: temperature is negative of NaN');
    return;
end


l=length;
dw=weight; %larval dry weight in mug (from Huebert and Peck, 2014)

Q10GER = 2.52; % from HP2011
SDA=0.10; % % specific dynamic action

dw_ref= 0.018521*l^(3.613981*exp(mean(t)*0.006266)); % reference weight(Hufnagl and Peck, 2011)

betta=0.60*(1-0.3*exp(-0.003*dw-dw_ref)); % % assimilation efficiency
gcmax=10^(1.72+1.02*log10(dw_ref/1000)); % max gut content in mug

      
GER_incr=2.5;
convers=(15.999*2)/1000; % convert  nmol O2 ind-1 h-1 into mug ind-1 h-1
Rs=0.0528* (dw^0.8859)* exp(t * 0.1046); %  Respiration rate is in nmol O2/h (from Moyano, 2016)
Rs=Rs*convers; % in mug/h
k_day=2; k_night=1;

GER    = (1.792 * l^(-0.828)*Q10GER ^((t-12)/10))*100; % gu evacuation rate (Hufnagl and Peck, 2011)
GER = GER*GER_incr; % % increment of the gut evacuation rate during active feeding

gc_curr = min(gc_before + for_consump*light, gcmax); % increase in the gut content %µg + µg/h * h
egain = min(gc_curr*GER/100, gc_curr); % the food processed in Hufnagl and Peck, 2011
if light==1
    eloss=egain*(1-betta+SDA*betta)+k_day*Rs;
else
    eloss=egain*(1-betta+SDA*betta)+k_night*Rs;
end

energy_budget=egain - eloss;
gc_after=max(0, gc_curr-(gc_curr * GER/100)); % % the reduction of the stomach content due to evacuation rate

