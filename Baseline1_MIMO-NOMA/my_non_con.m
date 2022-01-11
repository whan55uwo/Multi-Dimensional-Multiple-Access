function [c,ceq] = my_non_con(p_vector)
%% Iniatilization
global BS; % type: structure
%%
[f_tot, u_fun, rate_ue, rate_ue_norm, SINR, SINR_SIC] = data_rate_cal(p_vector);
%% Bower bound of data rate
c = zeros(2*BS.num_of_UE, 1);
for i = 1 : BS.num_of_UE
     c(i) = (BS.rate_min - rate_ue(i))/BS.rate_max; 
     c(i + BS.num_of_UE) = (rate_ue(i) - BS.rate_max)/BS.rate_max; 
end % end of "for i = 1 : BS.num_of_UE"
ceq = [];
end % end of function