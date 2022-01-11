function [f_tot, u_fun, rate_ue, rate_ue_norm, SINR, SINR_SIC] = data_rate_cal(p_vector)
%%
global BS; % type: structure
global Users;
%%
SINR = zeros(BS.num_of_UE, 1);
SINR_SIC = zeros(BS.num_of_UE, 1);
rate_ue = zeros(BS.num_of_UE, 1);
u_fun = zeros(BS.num_of_UE, 1);
f_tot = 0;
%% SINR and SINR at SIC recievier
for loop_UE_idx = 1 : BS.num_of_UE
    SC_index = Users(loop_UE_idx).SC_idx;
    h_k = Users(loop_UE_idx).H_matrix(:, SC_index);
    w_k = Users(loop_UE_idx).beam_vector;
    S_pow = p_vector(loop_UE_idx) * norm(h_k'*w_k)^2;
    %% if it is using SD-NOMA
    I_SD_NOMA = 0;
    if Users(loop_UE_idx).if_S_NOMA == 1
        Partner_idx_SD_NOMA = Users(loop_UE_idx).partner_S_NOMA;
        Partner_num_SD_NOMA = length(Partner_idx_SD_NOMA);
        for j = 1 : Partner_num_SD_NOMA
            UE_j_idx = Partner_idx_SD_NOMA(j);
            w_j = Users(UE_j_idx).beam_vector;
            I_SD_NOMA = I_SD_NOMA + p_vector(UE_j_idx) * norm(h_k'*w_j)^2;
        end % end of "for j = 1 : Partner_num_SD_NOMA"
    end
    %% if it is using PD-NOMA and it is a FAR UE
    I_PD_NOMA = 0;
    if Users(loop_UE_idx).if_P_NOMA == 1 && Users(loop_UE_idx).if_Strong_UE_P_NOMA == 0
        % this UE is a weak signal user
        UE_strong_idx = Users(loop_UE_idx).partner_P_NOMA;
        w_UE_strong = Users(UE_strong_idx).beam_vector;
        I_PD_NOMA = p_vector(UE_strong_idx)*norm(h_k'*w_UE_strong);
    end
    %% if it is using PD-NOMA and it is a NEAR UE
    if Users(loop_UE_idx).if_P_NOMA == 1 && Users(loop_UE_idx).if_Strong_UE_P_NOMA == 1
        % this UE is a strong signal user
        UE_weak_idx = Users(loop_UE_idx).partner_P_NOMA;
        w_UE_weak = Users(UE_weak_idx).beam_vector;
        pow_UE_weak = p_vector(UE_weak_idx);
        SIC_S_pow = pow_UE_weak*norm(h_k'*w_UE_weak)^2;
        SINR_SIC(loop_UE_idx) = SIC_S_pow/(S_pow + I_SD_NOMA + BS.noise_power);
    end
    %% SINR caluculation
    SINR(loop_UE_idx) = S_pow/(I_SD_NOMA + I_PD_NOMA + BS.noise_power);
    %% Data Rate of UE
    rate_ue(loop_UE_idx) = 0.5*BS.sys_SC_bandwidth*log2(1+ SINR(loop_UE_idx));
    %% Utylity function calculation of each UE
    if Users(loop_UE_idx).if_P_NOMA == 1 && Users(loop_UE_idx).if_Strong_UE_P_NOMA == 1
        u_fun(loop_UE_idx) = min(rate_ue(loop_UE_idx)/BS.rate_max, 1) - ...
            Users(loop_UE_idx).cost_weight*g_cost_cal(loop_UE_idx, SC_index, SINR_SIC(loop_UE_idx), 0);
    else
        u_fun(loop_UE_idx) = min(rate_ue(loop_UE_idx)/BS.rate_max, 1) - ...
            Users(loop_UE_idx).cost_weight*g_cost_cal(loop_UE_idx, SC_index, 0, 0);
    end
    f_tot = f_tot + u_fun(loop_UE_idx);
end % end of "for UE_idx = 1 : BS.num_of_UE"
rate_ue_norm = rate_ue/BS.rate_max;
end % end of function