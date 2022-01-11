function [p_c_vector_opt, p_p_vector_opt, a_vector_opt] = precoding_power_allocation()
global BS; % type: structure
global Users;
global Groups;
%% Initalization
[Groups(1:BS.num_of_SCs).H_mat] = deal(0); % channel condition of UEs in the group
%% Precoder Desgin
for group_idx = 1:BS.num_of_SCs
    SC_idx = group_idx;
    UE_group = Groups(group_idx).user_group;
    group_size = length(UE_group);
    if group_size == 1
        UE_idx = UE_group(1);
        h_vector =  Users(UE_idx).H_matrix(:, SC_idx);
        Groups(group_idx).H_mat = h_vector;
        Users(UE_idx).beam_vector_private = h_vector/norm(h_vector);
        Groups(group_idx).beam_vector_common = h_vector/norm(h_vector);
        continue;
    end
    %% Precoding for Common Message: Broadcast
    angle_set = zeros(group_size, 1);
    for j = 1 : group_size
        UE_idx = UE_group(j);
        angle_set(j) = Users(UE_idx).ang;
    end
    ang_brodcast = mean(angle_set);
    a_vector_temp = zeros(BS.Nt, 1);
    for i = 1:BS.Nt
        a_vector_temp(i) = exp(1i*sin(ang_brodcast - 2*pi*(i-1)/BS.Nt));
    end
    Groups(group_idx).beam_vector_common = a_vector_temp/norm(a_vector_temp);
    %% Precoding for Private Message
    H_mat = zeros(BS.Nt, group_size);
    for j = 1 : group_size
        UE_idx = UE_group(j);
        H_mat(:, j) = Users(UE_idx).H_matrix(:, SC_idx);
    end
    Groups(group_idx).H_mat = H_mat;
    W_ZF = H_mat*pinv(H_mat'*H_mat);
    W_ZF=W_ZF./repmat(sqrt(sum(W_ZF.*conj(W_ZF),1)),size(W_ZF,1),1); % normalization
    for j = 1 : group_size
        UE_idx = UE_group(j);
        Users(UE_idx).beam_vector_private = W_ZF(:, j);
    end
end % end of group_idx
%% Power allocation
% Input: optimization variable x_vector can be divided into two parts:
% (1) p_vector can be further divided into two parts:
%           1) p_c_vector is the transmit power for the common message;
%           2) p_p_vector is the transmit power for the private message
x0_vector = ones(BS.num_of_SCs + 2*BS.num_of_UE, 1)*(BS.tx_power/(BS.num_of_SCs + BS.num_of_UE));
% (2) a_vector: the precentage of common message allocated to each user in each group
x0_vector(BS.num_of_SCs + BS.num_of_UE + 1 : end) = 0;
for group_idx = 1:BS.num_of_SCs
    UE_group = Groups(group_idx).user_group;
    group_size = length(UE_group);
    for i = 1 : group_size
        UE_i_idx = UE_group(i);
        x0_vector(BS.num_of_SCs + BS.num_of_UE + UE_i_idx) = 1/group_size;
    end % for i = 1 : group_size
end
%% linear inequalities: maximum power constraint
A_mat = zeros(1, BS.num_of_SCs + 2*BS.num_of_UE);
a_ = ones(1, BS.num_of_SCs + BS.num_of_UE);
A_mat(1, 1 : BS.num_of_SCs + BS.num_of_UE) = a_;
b = BS.tx_power;
%% Linear equality constraints
A_eq_mat = zeros(BS.num_of_SCs, BS.num_of_SCs + 2*BS.num_of_UE);
b_eq_vector = zeros(BS.num_of_SCs, 1);
for group_idx = 1 : BS.num_of_SCs
    UE_group = Groups(group_idx).user_group;
    % Equality constraint of common message assinged to each UE in each group
    a_ = zeros(1, BS.num_of_UE);
    a_(UE_group) = 1;
    A_eq_mat(group_idx, BS.num_of_SCs + BS.num_of_UE + 1: end) = a_;
    b_eq_vector(group_idx) = 1;
end % end of "for group_idx = 1:BS.num_of_SCs"
% Equality constraint of private message: UE with NON-SIC capability cannot
% perform this stage, therefore the UE with no-SIC capability is given zero
% transmit power
for group_idx = 1 : BS.num_of_SCs
    UE_group = Groups(group_idx).user_group;
    group_size = length(UE_group);
    if group_size == 1 % single UE case
        UE_idx = UE_group(1);
        h_vector =  Users(UE_idx).H_matrix(:, SC_idx);
        w_c = Groups(group_idx).beam_vector_common;
        % power of single UE
        a_eq_vec = zeros(1, BS.num_of_SCs + 2*BS.num_of_UE);
        a_eq_vec(group_idx) = 1;
        A_eq_mat = [A_eq_mat; a_eq_vec];
        b_eq = (BS.sinr_max*BS.noise_power)/(norm(h_vector'*w_c)^2);
        b_eq_vector = [b_eq_vector; b_eq];
        a_eq_vec = zeros(1, BS.num_of_SCs + 2*BS.num_of_UE);
        a_eq_vec(BS.num_of_SCs + UE_idx) = 1;
        A_eq_mat = [A_eq_mat; a_eq_vec];
        b_eq_vector = [b_eq_vector; 0];
    else
        for i = 1 : group_size
            UE_i_idx = UE_group(i);
            if Users(UE_i_idx).if_SIC_capability == 0
                x0_vector(BS.num_of_SCs + UE_i_idx) = 0;
                a_eq_vec = zeros(1, BS.num_of_SCs + 2*BS.num_of_UE);
                a_eq_vec(BS.num_of_SCs + UE_i_idx) = 1;
                A_eq_mat = [A_eq_mat; a_eq_vec];
                b_eq_vector = [b_eq_vector; 0];
            end
        end % end of "for i = 1 : group_size"
    end
end % end of "for group_idx = 1:BS.num_of_SCs"
%% Lower bounds of variables
lb = ones(BS.num_of_SCs + 2*BS.num_of_UE, 1)*1e-2;
% avoid the fmincon give negetive pow
for UE_idx = 1 :  BS.num_of_UE
    if Users(UE_idx).if_SIC_capability == 0 || isempty(Users(UE_i_idx).partner)
        lb(UE_idx) = 0;
    end
end % end of "for i = 1 : BS.num_of_SCs + BS.num_of_UE"
%% Upper bounds of variables
ub = ones(BS.num_of_SCs + 2*BS.num_of_UE, 1);
for i = 1 : BS.num_of_SCs + BS.num_of_UE
    ub(i) = .9*BS.tx_power;
end % end of "for i = 1 : BS.num_of_SCs + BS.num_of_UE"
%% Nonlinear Constraint
nonlcon = @my_non_con;
%% Optimization Process
options = optimoptions('fmincon','Display','final-detailed','Algorithm','sqp', 'MaxFunctionEvaluations', 10000);
[x_vector_opt, fval] = fmincon(@obj_fun_TCOM,x0_vector, A_mat, b, A_eq_mat, b_eq_vector, lb, ub, nonlcon, options);
%% Output
p_c_vector_opt = x_vector_opt(1:BS.num_of_SCs);
p_p_vector_opt = x_vector_opt(BS.num_of_SCs + 1 : BS.num_of_SCs + BS.num_of_UE);
a_vector_opt = x_vector_opt(BS.num_of_SCs + BS.num_of_UE + 1: end);
[rate_ue_tot, rate_ue_com, rate_ue_private, rate_group_com_tot, SINR, SINR_SIC] = data_rate_cal(p_c_vector_opt, p_p_vector_opt, a_vector_opt);
for group_idx = 1:BS.num_of_SCs
    SC_idx = group_idx;
    UE_group = Groups(group_idx).user_group;
    group_size = length(UE_group);
    SINR_SIC_group = SINR_SIC(UE_group);
    [cost_UEs] = g_cost_cal(UE_group, SC_idx, SINR_SIC_group, 1);
    Groups(group_idx).pow_common = p_c_vector_opt(group_idx);
    Groups(group_idx).rate_common_total = rate_group_com_tot(group_idx);
    Groups(group_idx).rate_common_total_norm = Groups(group_idx).rate_common_total/BS.rate_max;
    Groups(group_idx).a_vector = a_vector_opt(UE_group);
    for i = 1 : group_size
        UE_i_idx = UE_group(i);
        Users(UE_i_idx).SINR = SINR(UE_i_idx);
        Users(UE_i_idx).SINR_dB = 10*log10(SINR(UE_i_idx));
        Users(UE_i_idx).SINR_SIC = SINR_SIC(UE_i_idx);
        if Users(UE_i_idx).SINR_SIC ~= 0
            Users(UE_i_idx).SINR_SIC_dB = 10*log10(SINR_SIC(UE_i_idx));
        end
        Users(UE_i_idx).rate_total = rate_ue_tot(UE_i_idx);
        Users(UE_i_idx).rate_total_norm = Users(UE_i_idx).rate_total/BS.rate_max;
        Users(UE_i_idx).u_fun = min(Users(UE_i_idx).rate_total_norm, 1) - Users(UE_i_idx).cost_weight*cost_UEs(i);
        Users(UE_i_idx).rate_common = rate_ue_com(UE_i_idx);
        Users(UE_i_idx).rate_private = rate_ue_private(UE_i_idx);
        Users(UE_i_idx).pow_private = p_p_vector_opt(UE_i_idx);
    end % end of "for i = 1 : group_size"
end % end of "for group_idx = 1:BS.num_of_SCs"
end % end of function