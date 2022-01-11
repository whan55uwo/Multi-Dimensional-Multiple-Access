function [cost_UEs] = g_cost_cal(UE_group, SC_idx, SINR_SIC, MODE)
%% INPUT:
% 1. User Group Set, it contains the idx of each UE
% 2. The idx of SC
% 3. SINR at SIC receiver
%% Import the user data
global Users;
global BS;
%% INATALIZATION
group_size = length(UE_group);
cost_UEs = zeros(group_size, 1);
if group_size > 1
    for i = 1 : group_size
        UE_i_idx = UE_group(i);
        h = Users(UE_i_idx).H_matrix(:, SC_idx);
        %% PD-NOMA FOR solving common message
        PD_ND = 0;
        if Users(UE_i_idx).if_SIC_capability == 1 %if this UE has SIC capability
            PD_ND =  BS.constant_G - BS.rho1*log10(min(SINR_SIC(i), BS.sinr_max));
        end
        %% SD-NOMA
        SD_ND = 0;
        for j = 1 : group_size
            if j == i
                continue;
            end
            UE_j_idx = UE_group(j);
            h_j = Users(UE_j_idx).H_matrix(:, SC_idx);
            SD_ND = SD_ND + BS.rho2*(norm(h'*h_j)/(norm(h)*norm(h_j)));
        end % end of "for j = 1 : Partner_num_SD_NOMA"
        cost_UEs(i) = PD_ND + SD_ND;
        if MODE == 1
            Users(UE_i_idx).PD_ND = PD_ND;
            Users(UE_i_idx).SD_ND = SD_ND;
            Users(UE_i_idx).utilization_cost = cost_UEs(i);
        end
    end % end of "for i = 1 : group_size"
else
    UE_i_idx = UE_group(1);
    if MODE == 1
        Users(UE_i_idx).PD_ND = 0;
        Users(UE_i_idx).SD_ND = 0;
        Users(UE_i_idx).utilization_cost = 0;
    end
end % end of "if group_size > 1"
end % end of function