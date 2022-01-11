function [power_vector, J_nm_value, Rate, Rate_low_bound, SINR_vector, SINR_SIC_vector, iteration_num] = ...
    fix_point_equation(Eta_vector,Mu_vector,Nu,Coalition_idx,SC_idx)
%% Initialzation
global BS; % type: structure
global Coalitions;
global Users
%% Pre-Processing
UE_set = Coalitions(Coalition_idx).user_subset;
UE_set_size = length(UE_set);
power_vector = ones(UE_set_size, 1);
SINR_vector = zeros(UE_set_size, 1);
SINR_SIC_vector = zeros(UE_set_size, 1);
Rate = zeros(UE_set_size, 1);
Rate_low_bound = zeros(UE_set_size, 1);
J_nm_value = 0;
iteration_num = 1;
iteration_num_max = 31;
%% Coalitions(Coalition_idx).MA_mode == 0
if UE_set_size == 1
    UE_idx = UE_set(1);
    h = Users(UE_idx).H_matrix(:, SC_idx);
    D = BS.D0(UE_idx)*(1/BS.rate_max + Eta_vector(UE_idx) - Mu_vector(UE_idx));
    pow_ = D/Nu;
    SINR_ = (norm(h'*(Users(UE_idx).beam_vector_matrix(:,SC_idx)))^2 * pow_)/BS.noise_power;
    if SINR_ > BS.sinr_max
        power_vector(1) = BS.sinr_max*BS.noise_power/(norm(h'*(Users(UE_idx).beam_vector_matrix(:,SC_idx)))^2);
        SINR_vector(1) = BS.sinr_max;
    elseif SINR_ < BS.sinr_th
        power_vector(1) = BS.sinr_th*BS.noise_power/(norm(h'*(Users(UE_idx).beam_vector_matrix(:,SC_idx)))^2);
        SINR_vector(1) = BS.sinr_th;
    else % SINR_ >= BS.sinr_th
        power_vector(1) = pow_;
        SINR_vector(1) = SINR_;
    end
end % if UE_set_size == 1
%% Coalitions(Coalition_idx).MA_mode == 1, i.e., Power-domain NOMA mode
if Coalitions(Coalition_idx).MA_mode == 1
    %% Strong signal (Near) UE
    for k = 1 : UE_set_size
        UE_k_idx = UE_set(k);
        if Users(UE_k_idx).if_Strong_UE_P_NOMA == 1
            h_k = Users(UE_k_idx).H_matrix(:, SC_idx);
            D_k = BS.D0(UE_k_idx)*(1/BS.rate_max + Eta_vector(UE_k_idx) - Mu_vector(UE_k_idx));
            pow_k_ = D_k/Nu;
            SINR_k_ = (norm(h_k'*(Users(UE_k_idx).beam_vector_matrix(:,SC_idx)))^2 * pow_k_)/BS.noise_power;
            if SINR_k_ > BS.sinr_max
                power_vector(k) = BS.sinr_max*BS.noise_power/(norm(h_k'*(Users(UE_k_idx).beam_vector_matrix(:,SC_idx)))^2);
                SINR_vector(k) = BS.sinr_max;
            elseif SINR_k_ < BS.sinr_th
                power_vector(k) = BS.sinr_th*BS.noise_power/(norm(h_k'*(Users(UE_k_idx).beam_vector_matrix(:,SC_idx)))^2);
                SINR_vector(k) = BS.sinr_th;
            else % SINR_ >= BS.sinr_th
                power_vector(k) = pow_k_;
                SINR_vector(k) = SINR_k_;
            end
        end % end of if strong signal UE
    end % end of "for k = 1 : UE_set_size"
    %% Weak signal (Far) UE
    for k = 1 : UE_set_size
        UE_k_idx = UE_set(k);
        if Users(UE_k_idx).if_Strong_UE_P_NOMA == 0
            h_k = Users(UE_k_idx).H_matrix(:, SC_idx);
            SINR_vector(k) = BS.sinr_max;
            % Interference from strong UE
            for i = 1 : UE_set_size
                if k == i
                    continue
                end
                UE_i_idx = UE_set(i);
                h_i = Users(UE_i_idx).H_matrix(:, SC_idx);
                w_i = Users(UE_i_idx).beam_vector_matrix(:,SC_idx);
                I_PD_NOMA = power_vector(i) * norm(h_k'*w_i)^2;
                I_SIC = power_vector(i) * norm(h_i'*w_i)^2;
                power_vector(k) = max(BS.sinr_max*(BS.noise_power + I_PD_NOMA)...
                    /(norm(h_k'*(Users(UE_k_idx).beam_vector_matrix(:,SC_idx)))^2), ...
                    BS.sinr_max*(BS.noise_power + I_SIC)...
                    /(norm(h_i'*(Users(UE_k_idx).beam_vector_matrix(:,SC_idx)))^2));
                SINR_SIC_vector(i) = power_vector(k)*(norm(h_i'*(Users(UE_k_idx).beam_vector_matrix(:,SC_idx)))^2)/...
                    (BS.noise_power + I_SIC);
            end
        end % end of if weak signal UE
    end % end of "for k = 1 : UE_set_size"
end % end of "if Coalitions(Coalition_idx).MA_mode == 1"
%% Coalitions(Coalition_idx).MA_mode == 2, Spatial-domain NOMA
if Coalitions(Coalition_idx).MA_mode == 2
    D_vector = zeros(UE_set_size, 1);
    G_mat =  zeros(UE_set_size, UE_set_size);
    if_break = 0;
    %% fix-point equation iteration
    while 1
        power_vector_pre = power_vector;
        %% SINR calculation
        for i = 1 : UE_set_size
            UE_i_idx = UE_set(i);
            D_vector(i) = BS.D0(UE_i_idx)*(1/BS.rate_max + Eta_vector(UE_i_idx) - Mu_vector(UE_i_idx));
            h_i = Users(UE_i_idx).H_matrix(:, SC_idx);
            S_pow = power_vector(i)*norm(h_i'*(Users(UE_i_idx).beam_vector_matrix(:,SC_idx)))^2;
            %% Interference: if it is using SD-NOMA
            I_SD_NOMA = 0;
            if Users(UE_i_idx).if_S_NOMA == 1
                Partner_idx_SD_NOMA = Users(UE_i_idx).partner_S_NOMA;
            end
            for k = 1 : UE_set_size
                if k == i
                    continue
                end
                UE_k_idx = UE_set(k);
                % SD-NOMA
                if Users(UE_i_idx).if_S_NOMA == 1
                    if ismember(UE_k_idx, Partner_idx_SD_NOMA)
                        w_k = Users(UE_k_idx).beam_vector_matrix(:,SC_idx);
                        I_SD_NOMA = I_SD_NOMA + power_vector(k) * norm(h_i'*w_k)^2;
                        G_mat(i, k) = BS.G0(UE_i_idx, UE_k_idx, SC_idx)*...
                            (1/BS.rate_max + Eta_vector(UE_i_idx) - Mu_vector(UE_i_idx));
                    end
                end % end of "if Users(UE_i_idx).if_S_NOMA == 1"
            end % end of "for i = 1 : UE_set_size"
            % SINR calculation
            SINR_vector(i) = S_pow/(I_SD_NOMA + BS.noise_power);
        end % end of "for k = 1 : UE_set_size"
        %% Renew Power allocation
        if if_break == 0
            for k = 1 : UE_set_size
                numerator = D_vector(k);
                denominator = Nu;
                for i = 1 : UE_set_size
                    if k == i
                        continue
                    end
                    denominator = denominator + G_mat(i,k)*SINR_vector(i)/(power_vector(i));
                end % end of "for i = 1 : UE_set_size"
                if denominator == 0
                    fprintf('Check the function "fix_point_equation.m".\n');
                    pause(inf);
                end
                power_vector(k) = numerator/denominator;
            end % end of "for k = 1 : UE_set_size"
        end % end of "if if_break == 0"
        %% Termination Conditons
        if sum(abs(power_vector - power_vector_pre))/sum(power_vector) <= 0.1 && if_break == 0
            if_break = 1;
        end
        if iteration_num > iteration_num_max && if_break == 0
            if_break = 1;
        end
        if if_break == 1 % re-adjust power level of UE, who exceeds the BS.rate_max (BS.sinr_max)
            count  = 0;
            for k = 1 : UE_set_size
                if SINR_vector(k) > BS.sinr_max
                    power_vector(k) = power_vector(k)*(BS.sinr_max/SINR_vector(k));
                    count =  count + 1;
                end % end of "if SINR_vecotr(k) > BS.sinr_max"
            end % "for k = 1 : UE_set_size"
            if count == 0
                break
            end
        end % end of "if if_break == 1"
        %% Renew iteration number
        iteration_num = iteration_num + 1;
    end % end of while
end
%% Coalitions(Coalition_idx).MA_mode == 3, Hybrid NOMA
if Coalitions(Coalition_idx).MA_mode == 3
    D_vector = zeros(UE_set_size, 1);
    E_vector = zeros(UE_set_size, 1);
    F_mat =  zeros(UE_set_size, UE_set_size);
    G_mat =  zeros(UE_set_size, UE_set_size);
    V_mat = zeros(UE_set_size, UE_set_size, UE_set_size);
    if_break = 0;
    %% fix-point equation iteration
    while 1
        power_vector_pre = power_vector;
        %% SINR calculation
        for i = 1 : UE_set_size
            UE_i_idx = UE_set(i);
            D_vector(i) = BS.D0(UE_i_idx)*(1/BS.rate_max + Eta_vector(UE_i_idx) - Mu_vector(UE_i_idx));
            if Users(UE_i_idx).if_P_NOMA == 1 && Users(UE_i_idx).if_Strong_UE_P_NOMA == 0
                % UE i is a weak signal (FAR) user and UE k is a strong signal (NEAR) user
                E_vector(i) = BS.rho1*Users(UE_i_idx).cost_weight/log(10);
            end
            h_i = Users(UE_i_idx).H_matrix(:, SC_idx);
            w_i = Users(UE_i_idx).beam_vector_matrix(:,SC_idx);
            S_pow = power_vector(i)*norm(h_i'*w_i)^2;
            %% Interference: if it is using SD-NOMA/PD-NOMA
            I_SD_NOMA = 0;
            if Users(UE_i_idx).if_S_NOMA == 1
                Partner_idx_SD_NOMA = Users(UE_i_idx).partner_S_NOMA;
            end
            for k = 1 : UE_set_size
                if k == i
                    continue
                end
                UE_k_idx = UE_set(k);
                % SD-NOMA
                if Users(UE_i_idx).if_S_NOMA == 1
                    if ismember(UE_k_idx, Partner_idx_SD_NOMA)
                        w_k = Users(UE_k_idx).beam_vector_matrix(:,SC_idx);
                        I_SD_NOMA = I_SD_NOMA + power_vector(k) * norm(h_i'*w_k)^2;
                        G_mat(i, k) = BS.G0(UE_i_idx, UE_k_idx, SC_idx)*...
                            (1/BS.rate_max + Eta_vector(UE_i_idx) - Mu_vector(UE_i_idx));
                    end
                end % end of "if Users(UE_i_idx).if_S_NOMA == 1"
                % PD-NOMA
                I_PD_NOMA = 0;
                if Users(UE_i_idx).if_P_NOMA == 1 && Users(UE_i_idx).if_Strong_UE_P_NOMA == 0 ...
                        && UE_k_idx == Users(UE_i_idx).partner_P_NOMA
                    % UE i is a weak signal (FAR) user and UE k is a strong signal (NEAR) user
                    h_k = Users(UE_k_idx).H_matrix(:, SC_idx);
                    w_k = Users(UE_k_idx).beam_vector_matrix(:,SC_idx);
                    I_PD_NOMA = power_vector(k) * norm(h_i'*w_k)^2;
                    G_mat(i,k) = BS.G0(UE_i_idx,UE_k_idx, SC_idx)*...
                        (1/BS.rate_max + Eta_vector(UE_i_idx) - Mu_vector(UE_i_idx));
                    F_mat(i, k) = BS.rho1*Users(UE_i_idx).cost_weight/log(10)*...
                        (norm(h_k'*w_k)^2/norm(h_k'*w_i)^2);
                end
                % PD-NOMA SIC
                if Users(UE_i_idx).if_P_NOMA == 1 && Users(UE_i_idx).if_Strong_UE_P_NOMA == 1 ...
                        && UE_k_idx == Users(UE_i_idx).partner_P_NOMA
                    % UE i is a strong signal (NEAR) user and UE k is a weak signal (FAR) user
                    w_k = Users(UE_k_idx).beam_vector_matrix(:,SC_idx);
                    SIC_S_pow = power_vector(k) * norm(h_i'*w_k)^2; % THE receive power at SIC receiver
                end
            end % end of "for k = 1 : UE_set_size"
            % SINR calculation
            SINR_vector(i) = S_pow/(I_PD_NOMA + I_SD_NOMA + BS.noise_power);
            % SINR_SIC calculation
            if Users(UE_i_idx).if_P_NOMA == 1 && Users(UE_i_idx).if_Strong_UE_P_NOMA == 1
                % UE i is a strong signal (NEAR) user
                I_SIC = S_pow;
                SINR_SIC_vector(i) = SIC_S_pow/(I_SIC + I_SD_NOMA + BS.noise_power);
            end
        end % end of "for i = 1 : UE_set_size"
        %% Renew Power allocation
        if if_break == 0
            for k = 1 : UE_set_size
                UE_k_idx = UE_set(k);
                numerator = D_vector(k) + E_vector(k);
                denominator = Nu;
                for i = 1 : UE_set_size
                    UE_i_idx = UE_set(i);
                    if k == i
                        continue
                    end
                    denominator = denominator + G_mat(i,k)*SINR_vector(i)/(power_vector(i)) ...
                        + F_mat(i, k)*SINR_SIC_vector(k)/(power_vector(i));
                    %                     if Users(UE_k_idx).if_P_NOMA == 1 && Users(UE_k_idx).if_Strong_UE_P_NOMA == 1 ...
                    %                         && UE_i_idx == Users(UE_k_idx).partner_P_NOMA
                    %                         fprintf('Check in function "fix_point_equation.m". \n');
                    %                     end
                    for ii = 1 : UE_set_size
                        if k == ii || i == ii
                            continue;
                        end
                        UE_ii_idx = UE_set(ii);
                        if Users(UE_i_idx).if_P_NOMA == 1 && Users(UE_i_idx).if_Strong_UE_P_NOMA == 1 ...
                                && UE_ii_idx == Users(UE_i_idx).partner_P_NOMA
                            h_i = Users(UE_i_idx).H_matrix(:, SC_idx);
                            w_k = Users(UE_k_idx).beam_vector_matrix(:,SC_idx);
                            w_ii = Users(UE_ii_idx).beam_vector_matrix(:,SC_idx);
                            V_mat(i,ii,k) = BS.rho1*Users(UE_i_idx).cost_weight/log(10)*...
                                (norm(h_i'*w_k)^2/norm(h_i'*w_ii)^2);
                            denominator = denominator +  V_mat(i,ii,k)*SINR_SIC_vector(i)/(power_vector(ii));
                        end
                    end % end of "for ii = 1 : UE_set_size"
                end % end of "for i = 1 : UE_set_size"
                if denominator == 0
                    fprintf('Check the function "fix_point_equation.m".\n');
                    pause(inf);
                end
                power_vector(k) = numerator/denominator;
            end % end of "for k = 1 : UE_set_size"
        end % end of "if if_break == 0"
        %% Termination Conditons
        if sum(abs(power_vector - power_vector_pre))/sum(power_vector) <= 0.1 && if_break == 0
            if_break = 1;
        end
        if iteration_num > iteration_num_max && if_break == 0
            if_break = 1;
        end
        if if_break == 1 % re-adjust power level of UE, who exceeds the BS.rate_max (BS.sinr_max)
            count  = 0;
            for k = 1 : UE_set_size
                if SINR_vector(k) > BS.sinr_max
                    power_vector(k) = power_vector(k)*(BS.sinr_max/SINR_vector(k));
                    count =  count + 1;
                end % end of "if SINR_vecotr(k) > BS.sinr_max"
            end % "for k = 1 : UE_set_size"
            if count == 0
                break
            end
        end % end of "if if_break == 1"
        %% Renew iteration number
        iteration_num = iteration_num + 1;
    end % end of while
end
%% Calculation of Function J_{n,m}
for k = 1 : UE_set_size
    UE_k_idx = UE_set(k);
    if Users(UE_k_idx).if_P_NOMA == 1 && Users(UE_k_idx).if_Strong_UE_P_NOMA == 1
        SINR_SIC = SINR_SIC_vector(k);
    else
        SINR_SIC = 0;
    end
    J_nm_value = J_nm_value + (1/BS.rate_max + Eta_vector(UE_k_idx) - Mu_vector(UE_k_idx))*BS.sys_SC_bandwidth*BS.a_k*log2(SINR_vector(k)) ...
        - Nu*power_vector(k) - Users(UE_k_idx).cost_weight*g_cost_cal(UE_k_idx, SC_idx, SINR_SIC);
    %     if (1/BS.rate_max + Eta_vector(UE_k_idx) - Mu_vector(UE_k_idx))*BS.sys_SC_bandwidth*BS.a_k*log2(SINR_vector(k)) ...
    %             - Nu*power_vector(k) - Users(UE_k_idx).cost_weight*g_cost_cal(UE_k_idx, SC_idx, SINR_SIC) > 1000
    %         fprintf('Check the function "fix_point_equation.m": Part exceeds 1000. \n');
    %     end
    Rate(k) = 0.5*BS.sys_SC_bandwidth*log2(1 + SINR_vector(k));
    Rate_low_bound(k) = 0.5*BS.sys_SC_bandwidth*(BS.a_k*log2(SINR_vector(k)) + BS.b_k);
end % end of "for k = 1 : UE_set_size"
if J_nm_value == -Inf
    fprintf('Check the function "fix_point_equation.m". \n');
    pause(inf);
end
end % end of function