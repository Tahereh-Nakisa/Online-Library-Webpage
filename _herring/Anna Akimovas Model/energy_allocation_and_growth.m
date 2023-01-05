% % % the growth routine of herring larvae
% % input variables: 
% % ebundget - energy gain from the enrgy budget routine ()
% % l_old - larval length at the previous step [in mm]
% % w_old - larval weight at the previous step [in mg]
% % t - temeprature
% % output variables: 
% % l_new - larval length  [in mm]
% % w_new - larval weight  [in mg]
% % CI - conditional index (ratio of w_new to the reference larval weight)


% % Akimova Anna, Thünen-Institute of Sea Fisheries, 04.01.2023


function [l_new, dw_new, CI]=...
    energy_allocation_and_growth(ebudget, l_old, dw_old, t)


    dw_ref = 0.018521  * l_old.^(3.613981*exp(t*0.006266)); % % reference dry weight following Hufnagl an Peck, 2011
    
    dw_new=dw_old+ebudget; % % ebudget can be positive or negative
    CI=dw_new/dw_ref;
    
    if dw_new > dw_ref  % (or CI>1)   
        l_new=(dw_new/0.018521)^(1/(3.613981*exp(t*0.006266))); % %energy allocation in growth in length
    else
        l_new=l_old; % % no growth in length
    end
    
    
 