function [m]=epr(ds,T,F,sp)
%calculate egg production rate (m) based on starting day (ds),
%time series of temperature (T), food (F) and species
 
%h: hatching time
%d1: n1-c1   (from beginning of n1 to begining of c1)
%d2: c1-c5
%d3: c5-c6

if (strcmp(sp,'cfin'))

    %temperature seems not control epr (Campbell et al., 2000)
    %need further test if temperature dependent
 
    kfood = 1.0;      %half saturation: 1 ug/l chla 
    m = 70.*F(ds)./(F(ds)+kfood);
end  
 
