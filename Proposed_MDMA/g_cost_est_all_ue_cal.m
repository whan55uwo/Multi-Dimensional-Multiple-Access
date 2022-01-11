function [sum_cost, utilization_cost_est_UEs] = g_cost_est_all_ue_cal(temp_UE_coalitions, MODE)
%% INPUT:
% 1. User Coalitioon Set, it contains the idx of each UE
% 2. MODE = 0, only calculate ESTIMATED ue g_cost based on location
% infomation in Problem P1
%     MODE = 1, SET THE VALUE OF: 
%     1) Users.if_P_NOMA;
%     2) Users.if_Strong_UE_P_NOMA; 
%     3) Users.partner_P_NOMA and Users.partner_S_NOMA;
%     4) Users.PD_ND and Users.SD_ND
%% Import the user data
global Users;
global BS;
global Coalitions;
%% Inialization
utilization_cost_est_UEs = zeros(BS.num_of_UE, 1);
for Coalition_idx = 1 : BS.num_of_SCs
    coalition_ue_set = find(temp_UE_coalitions(Coalition_idx, :) == 1);
    [cost_UEs_est] = g_cost_est_cal(coalition_ue_set, MODE);
    utilization_cost_est_UEs(coalition_ue_set) = cost_UEs_est;
    if MODE == 1
        Coalitions(Coalition_idx).user_subset = coalition_ue_set;
        len_ = length(coalition_ue_set);
        alpha = 0;
        beta = 0;
        if len_ == 1
            Coalitions(Coalition_idx).MA_mode = 0; %OMA
        else
            for i = 1 : len_
                UE_idx = coalition_ue_set(i);
                if Users(UE_idx).if_P_NOMA == 1
                    alpha = alpha + 1;
                end
                if Users(UE_idx).if_S_NOMA == 1
                    beta = beta + 1;
                end
            end % end of for i = 1 : len_
            if alpha > 0 && beta == 0
                Coalitions(Coalition_idx).MA_mode = 1; %1: P-NOMA mode
            end
            if alpha == 0 && beta > 0
                Coalitions(Coalition_idx).MA_mode = 2; %2: S-NOMA mode
            end
            if alpha >0 && beta > 0
                Coalitions(Coalition_idx).MA_mode = 3; %3: Hybrid NOMA mode
            end
        end % if len_ == 1
    end % if MODE == 1
end % for Coalition_idx = 1 : BS.num_of_SCs
sum_cost = sum(utilization_cost_est_UEs.*([Users(1:BS.num_of_UE).cost_weight])');
end