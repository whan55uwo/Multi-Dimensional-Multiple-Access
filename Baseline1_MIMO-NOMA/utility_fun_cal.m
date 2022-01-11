function [] = utility_fun_cal()
global BS; % type: structure
global Users;
%% Assign pow for PD-NOMA PAIRS
for loop_UE_idx = 1 : BS.num_of_UE
    if Users(loop_UE_idx).if_P_NOMA == 1 && Users(loop_UE_idx).if_Strong_UE_P_NOMA == 0
        % this UE is a weak signal user
        Users(loop_UE_idx).pow = Users(loop_UE_idx).pow*1.6;
    end
    if Users(loop_UE_idx).if_P_NOMA == 1 && Users(loop_UE_idx).if_Strong_UE_P_NOMA == 1
        % this UE is a strong signal user
        Users(loop_UE_idx).pow = Users(loop_UE_idx).pow*0.4;
    end
end
%% Calculation
while 1
    for loop_UE_idx = 1 : BS.num_of_UE
        SC_index = Users(loop_UE_idx).SC_idx;
        h_k = Users(loop_UE_idx).H_matrix(:, SC_index);
        w_k = Users(loop_UE_idx).beam_vector;
        S_pow = Users(loop_UE_idx).pow * norm(h_k'*w_k)^2;
        %% if it is using SD-NOMA
        I_SD_NOMA = 0;
        if Users(loop_UE_idx).if_S_NOMA == 1
            Partner_idx_SD_NOMA = Users(loop_UE_idx).partner_S_NOMA;
            Partner_num_SD_NOMA = length(Partner_idx_SD_NOMA);
            for j = 1 : Partner_num_SD_NOMA
                UE_j_idx = Partner_idx_SD_NOMA(j);
                w_j = Users(UE_j_idx).beam_vector;
                I_SD_NOMA = I_SD_NOMA + Users(UE_j_idx).pow * norm(h_k'*w_j)^2;
            end % end of "for j = 1 : Partner_num_SD_NOMA"
        end
        %% if it is using PD-NOMA and it is a FAR UE
        I_PD_NOMA = 0;
        if Users(loop_UE_idx).if_P_NOMA == 1 && Users(loop_UE_idx).if_Strong_UE_P_NOMA == 0
            % this UE is a weak signal user
            UE_strong_idx = Users(loop_UE_idx).partner_P_NOMA;
            w_UE_strong = Users(UE_strong_idx).beam_vector;
            I_PD_NOMA = Users(UE_strong_idx).pow * norm(h_k'*w_UE_strong);
        end
        %% if it is using PD-NOMA and it is a NEAR UE
        if Users(loop_UE_idx).if_P_NOMA == 1 && Users(loop_UE_idx).if_Strong_UE_P_NOMA == 1
            % this UE is a strong signal user
            UE_weak_idx = Users(loop_UE_idx).partner_P_NOMA;
            w_UE_weak = Users(UE_weak_idx).beam_vector;
            pow_UE_weak = Users(UE_weak_idx).pow;
            SIC_S_pow = pow_UE_weak*norm(h_k'*w_UE_weak)^2;
            Users(loop_UE_idx).SINR_SIC = SIC_S_pow/(S_pow + I_SD_NOMA + BS.noise_power);
        end
        %% SINR caluculation
        Users(loop_UE_idx).SINR = S_pow/(I_SD_NOMA + I_PD_NOMA + BS.noise_power);
        Users(loop_UE_idx).SINR_dB = 10*log10(Users(loop_UE_idx).SINR);
        %% Data Rate of UE
        Users(loop_UE_idx).rate = 0.5*BS.sys_SC_bandwidth*log2(1+Users(loop_UE_idx).SINR);
        %% Utylity function calculation of each UE
        if Users(loop_UE_idx).if_P_NOMA == 1 && Users(loop_UE_idx).if_Strong_UE_P_NOMA == 1
            Users(loop_UE_idx).u_fun = Users(loop_UE_idx).rate/BS.rate_max - ...
                Users(loop_UE_idx).cost_weight*g_cost_cal(loop_UE_idx, SC_index, Users(loop_UE_idx).SINR_SIC, 1);
        else
            Users(loop_UE_idx).u_fun = Users(loop_UE_idx).rate/BS.rate_max - ...
                Users(loop_UE_idx).cost_weight*g_cost_cal(loop_UE_idx, SC_index, 0, 1);
        end
    end % end of "for UE_idx = 1 : BS.num_of_UE"
    %% Re-allocate pow
    count = 0;
    for loop_UE_idx = 1 : BS.num_of_UE
        if Users(loop_UE_idx).SINR > BS.sinr_max
            Users(loop_UE_idx).pow = Users(loop_UE_idx).pow*(BS.sinr_max/Users(loop_UE_idx).SINR);
            count =  count + 1;
        end % end of "if SINR_vecotr(k) > BS.sinr_max"
    end
    %% Break condition
    if count == 0
        break
    end
end % end of while
end % end of function