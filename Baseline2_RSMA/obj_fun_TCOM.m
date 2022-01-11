function [f] = obj_fun_TCOM(x_vector)
global BS; % type: structure
global Users;
global Groups;
%% One-layer RS structure
% Input: p_vector can be divided into two parts:
%           1) p_c_vector is the transmit power for the common message;
%           2) p_p_vector is the transmit power for the private message
% a_vector: the precentage of common message allocated to each user in each group
p_c_vector = x_vector(1:BS.num_of_SCs);
p_p_vector = x_vector(BS.num_of_SCs + 1 : BS.num_of_SCs + BS.num_of_UE);
a_vector = x_vector(BS.num_of_SCs + BS.num_of_UE + 1: end);
%% Calcualte utility function and utilization cost

f = 0;
[rate_ue_tot, rate_ue_com, rate_ue_private, rate_group_com_tot, SINR, SINR_SIC] = data_rate_cal(p_c_vector, p_p_vector, a_vector);
for group_idx = 1:BS.num_of_SCs
    SC_idx = group_idx;
    UE_group = Groups(group_idx).user_group;
    group_size = length(UE_group);
    if group_size == 1
        UE_idx = UE_group(1);
        f = f + min(rate_ue_tot(UE_idx)/BS.rate_max, 1) - Users(UE_idx).cost_weight*Users(UE_idx).utilization_cost;
    else
        SINR_SIC_group = SINR_SIC(UE_group);
        [cost_UEs] = g_cost_cal(UE_group, SC_idx, SINR_SIC_group, 0);
        for i = 1 : group_size
            UE_i_idx = UE_group(i);
            f = f + min(rate_ue_tot(UE_i_idx)/BS.rate_max, 1) - Users(UE_i_idx).cost_weight*cost_UEs(i);
        end % end of "for i = 1 : group_size"
    end % end of "if group_size == 1"
end % end of "for group_idx = 1:BS.num_of_SCs"
%% Output
f = -f;
if ~isreal(f)
    fprintf('error in function "obj_fun_TWC.m"!\n');
end
end % end of function