function [cost_UEs] = g_cost_cal_est(UE_group, SC_index)
%% Import the user data
global Users;
global BS;
%% Initialization
group_size = length(UE_group);
cost_UEs = zeros(group_size, 1);
%% Calculation
SD_ND_of_each_UE = zeros(group_size, 1);
%PD_ND_est_of_each_UE = zeros(group_size, 1);
if group_size > 1
    for i = 1 : group_size
        UE_idx = UE_group(i);
        CH_UE_i = Users(UE_idx).H_matrix(:, SC_index);
%         CH_gain_max = norm(CH_UE_i)^2;
%         CH_gain_min = norm(CH_UE_i)^2;
        for j = 1 : group_size
            UE_partner_idx = UE_group(j);
            CH_UE_j = Users(UE_partner_idx).H_matrix(:, SC_index);
            if UE_partner_idx ==  UE_idx
                continue;
            end % end of if UE_idx == UE_partner_idx
%             CH_gain_max = max(norm(CH_UE_i)^2, CH_gain_max);
%             CH_gain_min = min(norm(CH_UE_i)^2, CH_gain_min);
            Partner_idx_SD_NOMA = [];
            %% SD-NOMA
            Partner_idx_SD_NOMA = [Partner_idx_SD_NOMA; UE_partner_idx];
            SD_ND_of_each_UE(i) = SD_ND_of_each_UE(i) +  norm(CH_UE_i'*CH_UE_j)/(norm(CH_UE_i)*norm(CH_UE_j));
            %% PD
%             PD_ND_est_of_each_UE(i) = ...
%                         max(BS.constant_G + BS.rho1*log10(CH_gain_min/CH_gain_max), 0);
        end % end of "for j = 1 : group_size"
        cost_UEs(i) = SD_ND_of_each_UE(i);
    end % end of "for i = 1 : group_size"
end % end of "if group_size > 1"
end
