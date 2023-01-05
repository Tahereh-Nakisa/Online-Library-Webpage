% % The main script for simulation of a daily growth of a herring larvae
% % of different initial length, at various temperatures, total zooplankton biomass and
% size-spectra slopes

% % Akimova Anna, Thünen-Institute of Sea Fisheries, 04.01.2023


clear all; close all;

%l_init=[8:1:13]; % initial larval length
l_init=13;
%t=[10:2:16]; % % temperature
t=12;

min_l=0.01; % % lower prey length (mm)
max_l=2.00; % % higher prey length (mm)
bin_num=200; % % number of prey size bins
%biomass=[0.01:0.2:1.9 2:1:59 60:25:1000 1010:100:5000 10000]; % total prey biomass in mg/m3
biomass=[10000]; % total prey biomass in mg/m3
%slope=[-2.0:0.1:0]; % % prey size spectra slope
slope=[-1.5];
exp_dur=48; % the duration of the simulation in hours (24 hours spin out and 24 hours of larval growth)
light_hours=12; % % daylength(in h)

a=length(l_init);
b=length(t);
c=length(biomass);
d=length(slope);

fish_l=NaN(a,b,c,d,exp_dur); % % initiate fish lenght array [in mm]
fish_w=NaN(a,b,c,d,exp_dur); % % initiate fish weight array [in mug]

light=zeros(exp_dur,1);
for n=1:round(exp_dur/24) 
    light(1+24*(n-1):11+24*(n-1))=1;
end

gc(1:a,1:b,1:c,1:d,1)=0; % % empty gut before the spin off simulations

for dd=1:d
    for cc=1:c
        [bbins, cbins, lbins, bin_mean_l, bin_mean_dw]=...
                size_spectrum_zooplankton(min_l, max_l, bin_num, biomass(cc), slope(dd));
        for bb=1:b
            fish_l(1:a,bb,cc,dd,1)=l_init'; % initial fish length
            for aa=1:a
                fish_w(aa,bb,cc,dd,1)=0.018521.*l_init(aa).^(3.613981.*exp(t(bb).*0.006266)); % initial fish weight
                gc=0; % % gut content is zero before the spin off simulations
                CI=1; % % conditional index
                for m=2:exp_dur 
                    [cons, ops, mposps, maxps, minps]=optimal_foraging(fish_l(aa,bb,cc,dd,m-1), bin_mean_l, bin_mean_dw, cbins);           
                     cons=cons*1000; % convert consumption in mug/h
                    [ebudget, egain, eloss, gc]=...
                     energy_budget(fish_l(aa,bb,cc,dd,m-1), fish_w(aa,bb,cc,dd,m-1), light(m-1), cons, gc, t(bb));
                    clear cons
                    if (m-(exp_dur-24))<=0 %  spin off to establish a realistic gut content (gc)
                        fish_w(aa,bb,cc,dd,m)=fish_w(aa,bb,cc,dd,m-1);
                        fish_l(aa,bb,cc,dd,m)=fish_l(aa,bb,cc,dd,m-1);
                    else
                        [fish_l(aa,bb,cc,dd,m), fish_w(aa,bb,cc,dd,m), CI]=...
                                        energy_allocation_and_growth(ebudget, fish_l(aa,bb,cc,dd,m-1), fish_w(aa,bb,cc,dd,m-1), t(bb));          
                    end              
                    clear energy_bunget egain eloss CI
                end
                %[m aa bb cc dd fish_w(aa,bb,cc,dd,1) fish_w(aa,bb,cc,dd,m)  fish_l(aa,bb,cc,dd,1) fish_l(aa,bb,cc,dd,m)]
            end
        end
    end
end


