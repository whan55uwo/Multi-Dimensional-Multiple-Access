function [cost_UE] = g_cost_cal(UE_idx, SC_idx, SINR_SIC)
%% INPUT:
% 1. UE index
% 2. Index of SC
% 3. If possible, SINR at SIC receiver is required
%% Import the user data
global Users;
global BS;
%% Processing
if Users(UE_idx).if_P_NOMA == 0 && Users(UE_idx).if_S_NOMA == 0
    Users(UE_idx).PD_ND = 0;
    Users(UE_idx).SD_ND = 0;
    cost_UE = 0;
    Users(UE_idx).utilization_cost = 0;
    return
end
%% INATALIZATION
h = Users(UE_idx).H_matrix(:, SC_idx);
%% PD-NOMA
PD_ND = 0;
if Users(UE_idx).if_P_NOMA == 1 && Users(UE_idx).if_Strong_UE_P_NOMA == 1 %if this UE is a strong UE
    PD_ND =  BS.constant_G - BS.rho1*log10(min(SINR_SIC, BS.sinr_max));
end
Users(UE_idx).PD_ND = PD_ND;
%% SD-NOMA
SD_ND = 0;
if Users(UE_idx).if_S_NOMA == 1
    Partner_idx_SD_NOMA = Users(UE_idx).partner_S_NOMA;
    Partner_num_SD_NOMA = length(Partner_idx_SD_NOMA);
    for j = 1 : Partner_num_SD_NOMA
        UE_j_idx = Partner_idx_SD_NOMA(j);
        h_j = Users(UE_j_idx).H_matrix(:, SC_idx);
        SD_ND = SD_ND + BS.rho2*(norm(h'*h_j)/(norm(h)*norm(h_j)));
    end % end of "for j = 1 : Partner_num_SD_NOMA"
end
Users(UE_idx).SD_ND = SD_ND;
cost_UE = PD_ND + SD_ND;
Users(UE_idx).utilization_cost = cost_UE;
end % end of function