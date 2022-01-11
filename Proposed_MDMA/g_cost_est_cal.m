function [cost_UEs] = g_cost_est_cal(UE_coalition, MODE)
%% INPUT:
% 1. User Coalitioon Set, it contains the idx of each UE
% 2. MODE = 0, only calculate ESTIMATED ue g_cost based on location
% infomation in Problem P1
%     MODE = 1, SET THE VALUE OF: 1) Users.if_P_NOMA;
%     2) Users.if_Strong_UE_P_NOMA; 3) Users.partner_P_NOMA and Users.partner_S_NOMA;
%     4) Users.PD_ND and Users.SD_ND
%% Import the user data
global Users;
global BS;
%% Initialization
coalition_size = length(UE_coalition);
beam_idx_of_each_UE = zeros(coalition_size, 1);
cost_UEs = zeros(coalition_size, 1);
for i = 1 : coalition_size
    UE_idx = UE_coalition(i);
    beam_idx_of_each_UE(i) = Users(UE_idx).beam_idx;
end
%% Caluation
PD_ND_est_of_each_UE = zeros(coalition_size, 1);
SD_ND_est_of_each_UE = zeros(coalition_size, 1);
if coalition_size > 1
    for i = 1 : coalition_size
        UE_idx = UE_coalition(i);
        Partner_idx_PD_NOMA = [];
        Partner_idx_SD_NOMA = [];
        if_PD_NOMA = 0;
        if_strong_UE = 0;
        for j = 1 : coalition_size
            UE_partner_idx = UE_coalition(j);
            if UE_partner_idx ==  UE_idx
                continue;
            end % end of if UE_idx == UE_partner_idx
            if beam_idx_of_each_UE(i) == beam_idx_of_each_UE(j)
                %% PD-NOMA
                if_PD_NOMA = 1;
                Partner_idx_PD_NOMA = [Partner_idx_PD_NOMA; UE_partner_idx];
                if Users(UE_partner_idx).dist > Users(UE_idx).dist
                    if_strong_UE = 1;
                    PD_ND_est_of_each_UE(i) = ...
                        BS.constant_G + BS.rho1*log10(max((Users(UE_partner_idx).pathloss + BS.noise_power/BS.tx_power)/Users(UE_idx).pathloss, 1/BS.sinr_max));
                end % determine if it is strong user
            else
                %% SD-NOMA
                Partner_idx_SD_NOMA = [Partner_idx_SD_NOMA; UE_partner_idx];
                SD_ND_est_of_each_UE(i) = SD_ND_est_of_each_UE(i) +  ...
                    BS.rho2 * norm((Users(UE_idx).a_vector)'*Users(UE_partner_idx).a_vector);
            end
        end % end of "for j = 1 : group_size"
        cost_UEs(i) = PD_ND_est_of_each_UE(i) + SD_ND_est_of_each_UE(i);
        if  BS.rho1*PD_ND_est_of_each_UE(i) + BS.rho2*SD_ND_est_of_each_UE(i) == 0 && if_strong_UE == 1
            fprintf('Check in function g_cost_est_cal.m. \n');
            pause(inf);
        end
        if length(Partner_idx_PD_NOMA) > 2
            fprintf('Check in function g_cost_est_cal.m. \n');
            pause(inf);
        end
        if MODE == 1
            Users(UE_idx).if_P_NOMA = if_PD_NOMA;
            Users(UE_idx).if_S_NOMA = (~isempty(Partner_idx_SD_NOMA));
            Users(UE_idx).if_Strong_UE_P_NOMA = if_strong_UE;
            Users(UE_idx).partner_P_NOMA = Partner_idx_PD_NOMA;
            Users(UE_idx).partner_S_NOMA = Partner_idx_SD_NOMA;
            Users(UE_idx).PD_ND_est = PD_ND_est_of_each_UE(i);
            Users(UE_idx).SD_ND_est = SD_ND_est_of_each_UE(i);
            Users(UE_idx).utilization_cost_est = cost_UEs(i);
        end
    end % end of "for i = 1 : group_size"
else
    % Group_size = 1, that is, only one UE in this SC
    UE_idx = UE_coalition(1);
    if MODE == 1
        Users(UE_idx).if_P_NOMA = 0;
        Users(UE_idx).if_S_NOMA = 0;
        Users(UE_idx).if_Strong_UE_P_NOMA = 0;
        Users(UE_idx).partner_P_NOMA = [];
        Users(UE_idx).partner_S_NOMA = [];
        Users(UE_idx).PD_ND_est = 0;
        Users(UE_idx).SD_ND_est = 0;
        Users(UE_idx).utilization_cost_est = 0;
    end
end % end of "if group_size > 1"
end