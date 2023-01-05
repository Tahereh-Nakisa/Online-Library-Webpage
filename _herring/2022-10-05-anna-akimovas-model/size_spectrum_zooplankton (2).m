% % create a size-spectrum of zooplankton
% % based on the input parameters:
% %
% % min_l - minimum size of prey (mm)
% % max_l - maximum size of prey (mm)
% % slope - the slope of the size-spectrum 
% % bin_num - the number of the size bins in the size-spectra
% % biomass - total biomass (mg/m3) of all zooplankton organisms

% % the output variables area: 
% % b_bins - zooppankton biomass in the size bins
% % c_bins - zooplankton abundance in the size bins
% % length_bins - bins bounds in length 
% % bin_mean_l - mean length in the bins
% % bin_mean_dw - mean dry weight in the bins

% % Akimova Anna, Thünen-Institute of Sea Fisheries, 04.01.2023

function [b_bins, c_bins, length_bins, bin_mean_l, bin_mean_dw]= ...
    size_spectrum_zooplankton(min_l, max_l, bin_num, biomass, slope)

l1=min_l; 
l2=max_l; 
tb=biomass; 
s=slope;

x=linspace(l1, l2, bin_num); % in mm
bin_mean_l=x(1:end-1)+diff(x)/2; % in mm 
length_bins=x;

dw=14*x.^2.5; % in mug see: Huebert et al, 2018
bin_mean_dw=14*bin_mean_l.^2.5; % in mug see: Huebert et al, 2018
dw=dw*1e-3; % % convert to mg
bin_mean_dw=bin_mean_dw*1e-3; % % convert to mg

% % Calculate normalized size-spectra beased on Blanco 1994 and corrected
% Huebert et al, 2018
if (s==-1) 
    b0_n=(log(dw(end)) - log(dw(1)));
else
    b0_n=(1/(s+1))*(dw(end)^(s+1) - dw(1)^(s+1));
end
        
        
for k=1:length(x)-1
    if (s==-1) 
        b_bin(k)=( log(dw(k+1)) - log(dw(k)) );
    else
        b_bin(k)=(1/(s+1))*( dw(k+1)^(s+1) - dw(k)^(s+1) );
    end
    tb_bin(k)=tb*b_bin(k)/b0_n; % % in mg/m3
end
        
b_bins=tb_bin; % % biomass distribution in mg/m3        
c_bins=tb_bin./bin_mean_dw; % % in [num/m3]



        
        

