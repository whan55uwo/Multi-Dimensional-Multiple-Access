function [rate_ue_tot, rate_ue_com, rate_ue_private, rate_group_com_tot, SINR, SINR_SIC] = data_rate_cal(p_c_vector, p_p_vector, a_vector)
%%
global BS; % type: structure
global Users;
global Groups;
%%
SINR = zeros(BS.num_of_UE, 1);
SINR_SIC = zeros(BS.num_of_UE, 1);
rate_ue_tot = zeros(BS.num_of_UE, 1);
rate_ue_com = zeros(BS.num_of_UE, 1);
rate_ue_private = zeros(BS.num_of_UE, 1);
rate_group_com_tot = zeros(BS.num_of_SCs, 1);
%% SINR and SINR at SIC recievier
for group_idx = 1:BS.num_of_SCs
    SC_idx = group_idx;
    UE_group = Groups(group_idx).user_group;
    group_size = length(UE_group);
    w_c = Groups(group_idx).beam_vector_common;
    %% COMMON MESSAGE RECEIVEING
    % Single UE case
    if group_size == 1
        UE_idx = UE_group(1);
        h_vector =  Users(UE_idx).H_matrix(:, SC_idx);
        SINR(UE_idx) = p_c_vector(group_idx)*norm(h_vector'*w_c)^2/BS.noise_power;
        rate_group_com_tot(group_idx) = 0.5*BS.sys_SC_bandwidth*log2(1 + SINR(UE_idx));
        rate_ue_tot(UE_idx) = rate_group_com_tot(group_idx);
        rate_ue_com(UE_idx) = rate_group_com_tot(group_idx);
        continue;
    end
    % Multiple UEs case
    I_pow_p = zeros(group_size, group_size); % interference between private UEs
    count = 1;
    for i = 1 : group_size
        UE_i_idx = UE_group(i);
        h_i_vector =  Users(UE_i_idx).H_matrix(:, SC_idx);
        S_pow_c = p_c_vector(group_idx)*norm(h_i_vector'*w_c)^2;
        I_pow_c = 0;
        for j = 1 : group_size
            UE_j_idx = UE_group(j);
            if Users(UE_j_idx).if_SIC_capability == 1
                w_p = Users(UE_j_idx).beam_vector_private;
                I_pow_p(i, j) = p_p_vector(UE_j_idx)*norm(h_i_vector'*w_p)^2;
                I_pow_c = I_pow_c + I_pow_p(i, j);
            end
        end % end of "for j = 1 : group_size"
        SINR_SIC(UE_i_idx) = S_pow_c/(I_pow_c + BS.noise_power);
        if Users(UE_i_idx).if_SIC_capability == 0
            SINR(UE_i_idx) = SINR_SIC(UE_i_idx);
        end
        if count == 1
            rate_group_com_tot(group_idx) = 0.5*BS.sys_SC_bandwidth*log2(1 + SINR_SIC(UE_i_idx));
        else
            rate_group_com_tot(group_idx) = min(rate_group_com_tot(group_idx), 0.5*BS.sys_SC_bandwidth*log2(1 + SINR_SIC(UE_i_idx)));
        end
        count = count + 1;
    end % end of "for i = 1 : group_size"
    %% COMBINE COMMON and PRIVATE MESSAGE from UEs
    for i = 1 : group_size
        UE_i_idx = UE_group(i);
        rate_ue_com(UE_i_idx) = rate_group_com_tot(group_idx)*a_vector(UE_i_idx);
        if Users(UE_i_idx).if_SIC_capability == 1
            Users(UE_i_idx).pow_private = p_p_vector(UE_i_idx);
            S_pow_p = I_pow_p(i,i);
            I_pow = sum(I_pow_p(i,:)) - S_pow_p;
            SINR(UE_i_idx) = S_pow_p/(I_pow + BS.noise_power);
            rate_ue_private(UE_i_idx) = 0.5*BS.sys_SC_bandwidth*log2(1 + SINR(UE_i_idx));
        else
            rate_ue_private(UE_i_idx) = 0;
        end % end of "if Users(UE_i_idx).if_SIC_capability == 1"
        rate_ue_tot(UE_i_idx) = min(rate_ue_com(UE_i_idx) + rate_ue_private(UE_i_idx), BS.rate_max);
    end % end of "for i = 1 : group_size"
end % end of "for group_idx = 1:BS.num_of_SCs"
end % end of function