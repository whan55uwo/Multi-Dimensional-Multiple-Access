function [c,ceq] = my_non_con(x_vector)
%% Input: p_vector can be divided into two parts:
%           1) p_c_vector is the transmit power for the common message;
%           2) p_p_vector is the transmit power for the private message
% a_vector: the precentage of common message allocated to each user in each group
%% Iniatilization
global BS; % type: structure
p_c_vector = x_vector(1:BS.num_of_SCs);
p_p_vector = x_vector(BS.num_of_SCs + 1 : BS.num_of_SCs + BS.num_of_UE);
a_vector = x_vector(BS.num_of_SCs + BS.num_of_UE + 1: end);
c = zeros(2*BS.num_of_UE, 1);
[rate_ue_tot, rate_ue_com, rate_ue_private, rate_group_com_tot, SINR, SINR_SIC] = data_rate_cal(p_c_vector, p_p_vector, a_vector);
%% Bound of data rate
for i = 1 : BS.num_of_UE
     c(i) = (BS.rate_min - rate_ue_tot(i))/BS.rate_max; 
     c(i + BS.num_of_UE) = (rate_ue_tot(i) - BS.rate_max)/BS.rate_max; 
end % end of "for i = 1 : BS.num_of_UE"
ceq = [];
end % end of function