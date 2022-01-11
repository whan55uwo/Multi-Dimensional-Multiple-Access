function [] = User_Group_Formation_Process()
global BS; % type: structure
global Users;
global Groups;
%% Output:
% 1. the MA mode of each UE
% 2. The UEs share the same SC with UE in different radio resource domain
% 3. SC allocation
% from “Location-Based MIMO-NOMA: Multiple Access Regions and Low-Complexity User Pairing”------Algorithm 1
%% Theorem 1: Corollary 1-Corollary 2
BS.pow_allocate_coeff =  0.5;
rho =  100; % power to noise
dist_array = [Users(1:BS.num_of_UE).dist];
[Users(1:BS.num_of_UE).delta_beta] = deal(0);
[Users(1:BS.num_of_UE).delta_d] = deal(0);
[Users(1:BS.num_of_UE).R_N_set] = deal(0);
for u = 1 : BS.num_of_UE
    Users(u).delta_beta = 0.8*(1/BS.pow_allocate_coeff - 3+sqrt((3 - 1/BS.pow_allocate_coeff)^2 + 4*BS.Nt*rho*(Users(u).dist^(-3.76))))/...
        (2*BS.Nt*rho*(Users(u).dist^(-3.76)));
    %% Power-domain NOMA
    if Users(u).if_SIC_capability == 1
        Users(u).delta_d = 10^((10 + 37.6*log10(Users(u).dist/1000))/37.6)*1000 - Users(u).dist;
        R_N_set_ = find(dist_array >= Users(u).dist + Users(u).delta_d);
        if isempty(R_N_set_)
            Users(u).R_N_set = [];
        else
            distance_ = dist_array(R_N_set_);
            [B, index_] = sort(distance_,'descend');
            R_N_set_ = R_N_set_(index_);
            R_N_set = R_N_set_;
            len_ = length(R_N_set_);
            for i = 1 : len_
                ue_i1_idx = R_N_set_(i);
                ND_S_est = Spatial_Non_est_cal(u, ue_i1_idx);
                if ND_S_est <= Users(u).delta_beta
                    R_N_set = setdiff(R_N_set, [ue_i1_idx]);
                end
            end % end of for i = 1 : len_
            Users(u).R_N_set = R_N_set;
        end % end of if
    else
        Users(u).delta_d = Inf;
        R_N_set_ = [];
        Users(u).R_N_set = R_N_set_;
    end % end of "if Users(u).if_SIC_capability == 1"
end % end of "for u = 1 : BS.num_of_UE"
%% Processing
beta_TH = 0.1;
%% Spatial-domain NOMA
[Users(1:BS.num_of_UE).R_S_set] = deal(1);
[Users(1:BS.num_of_UE).R_S_value_set] = deal(0);
for u = 1 : BS.num_of_UE
    R_S_set_ = [];
    Users(u).R_S_value_set = zeros(BS.num_of_UE, 1);
    for j = 1 : BS.num_of_UE
        if  j == u
            continue
        end
        ND_S_est = Spatial_Non_est_cal(u, j);
        Users(u).R_S_value_set(j) = ND_S_est;
        if ND_S_est <= beta_TH
            R_S_set_ = union(R_S_set_, [j]);
        end
    end % end of "for j = 1 : BS.num_of_UE"
    Users(u).R_S_set = R_S_set_;
end
%% Processing, Algorithm 1
T1 = [1 : BS.num_of_UE];
matched_indicator = zeros(BS.num_of_UE, 1);
max_pair_size = 3; % 3 pairs of NOMA in each group
for group_idx = 1 : BS.num_of_SCs
    %% break condition
    if isempty(T1)
        break;
    end
    %% find a NOMA pair for current group
    group_ = zeros(2*max_pair_size, 1);
    pair_idx = 1;
    while 1
        % break condition
        if  pair_idx > max_pair_size || isempty(T1)
            break;
        end
        % add a NOMA pair
        ue_i1_idx = T1(randperm(length(T1),1)); % selected as strong UE in the NOMA pair
        matched_indicator(ue_i1_idx) = 1;
        T2 = intersect(T1, Users(ue_i1_idx).R_N_set);
        if isempty(T2)
            group_(2*(pair_idx - 1) + 1, 1) = ue_i1_idx;
            group_(2*pair_idx, 1) = 0;
            T1 = intersect(T1, Users(ue_i1_idx).R_S_set);
        else
            ue_i2_idx = T2(1);
            matched_indicator(ue_i2_idx) = 1;
            group_(2*pair_idx - 1, 1) = ue_i1_idx;
            group_(2*pair_idx, 1) = ue_i2_idx;
            T1 = intersect(T1, Users(ue_i1_idx).R_S_set);
            T1 = intersect(T1, Users(ue_i2_idx).R_S_set);
        end % end of "if isempty(T2)"
        % Renew pair idx in current group
        pair_idx = pair_idx + 1;
    end % inner iteration
    %% Renew T1 for outer iteration
    matched_UE_set = find(matched_indicator == 1);
    T1 = setdiff([1 : BS.num_of_UE],  matched_UE_set);
    %% Value assignment
    Groups(group_idx).user_subset = group_;
    %     Groups(group_idx).user_subset_prev = group_;
    Groups(group_idx).SC_idx = group_idx;
end % end of "for group_idx = 1 : BS.num_of_SCs"
%% Re-assign UE if it is not paired
if sum(matched_indicator) ~= BS.num_of_UE
    for ue_idx = 1 : BS.num_of_UE
        if matched_indicator(ue_idx) == 0
            %% if this UE has not been assgined SC
            canadidate_group = inf*ones(BS.num_of_SCs, 1);
            ue_single_idx =  ue_idx;
            for  i = 1 : BS.num_of_SCs
                if length(Groups(i).user_subset) ~= 6
                    continue;
                end
                if Groups(i).user_subset(3) == 0 || Groups(i).user_subset(5) == 0
                    canadidate_group(i) = 0;
                    for j = 1 : 6
                        group_ue_idx = Groups(i).user_subset(j);
                        if group_ue_idx == 0
                            continue;
                        end
                        canadidate_group(i) = canadidate_group(i) + ...
                            Users(ue_single_idx).R_S_value_set(group_ue_idx);
                    end % end of "for j = 1 : 6"
                end % end of "if Groups(i).user_subset(3) == 0 || Groups(i).user_subset(5) == 0"
            end % end of "for  i = 1 : group_idx - 1"
            [Min_value, opt_index] = min(canadidate_group);
            if Min_value < inf
                if Groups(opt_index).user_subset(3) == 0
                    Groups(opt_index).user_subset(3) = ue_single_idx;
                else
                    Groups(opt_index).user_subset(5) = ue_single_idx;
                end
            end % end of "if Min_value < inf"
        end % end of "if matched_indicator(ue_idx) == 1"
    end % end of "for ue_idx = 1 : BS.num_of_UE"
end % end of if
%% Determine MA mode of UE and its paried ue in the same SC
for group_idx = 1 : BS.num_of_SCs
    group_ = Groups(group_idx).user_subset;
    SC_index = Groups(group_idx).SC_idx;
    if sum(group_) == 0
        continue;
    end
    alpha = 0;
    beta = 0;
    for pair_idx = 1:3
        NOMA_pair =  group_(pair_idx*2-1:pair_idx*2);
        ue_strong_idx = NOMA_pair(1);
        if ue_strong_idx == 0
            continue;
        end
        ue_weak_idx = NOMA_pair(2);
        if ue_weak_idx == 0
            other_ue_set = setdiff(group_, [0]);
            other_ue_set = setdiff(other_ue_set, [ue_strong_idx]);
            % single UE case
            Users(ue_strong_idx).if_P_NOMA = 0;
            Users(ue_strong_idx).if_Strong_UE_P_NOMA = 0;
            Users(ue_strong_idx).partner_P_NOMA = [];
            Users(ue_strong_idx).SC_idx = SC_index;
            Users(ue_strong_idx).beam_vector = (Users(ue_strong_idx).H_matrix(:, SC_index))/...
                norm(Users(ue_strong_idx).H_matrix(:, SC_index));
            Users(ue_strong_idx).if_S_NOMA = (~isempty(other_ue_set));
            Users(ue_strong_idx).partner_S_NOMA = other_ue_set;
        else
            other_ue_set = setdiff(group_, [0]);
            other_ue_set = setdiff(other_ue_set, [ue_strong_idx; ue_weak_idx]);
            % Non-single UE case
            % Strong signal UE
            Users(ue_strong_idx).if_P_NOMA = 1;
            Users(ue_strong_idx).if_Strong_UE_P_NOMA = 1;
            Users(ue_strong_idx).partner_P_NOMA = [ue_weak_idx];
            Users(ue_strong_idx).SC_idx = SC_index;
            Users(ue_strong_idx).if_S_NOMA = (~isempty(other_ue_set));
            Users(ue_strong_idx).partner_S_NOMA = other_ue_set;
            % Weak Signal UE
            Users(ue_weak_idx).if_P_NOMA = 1;
            Users(ue_weak_idx).if_Strong_UE_P_NOMA = 0;
            Users(ue_weak_idx).partner_P_NOMA = [ue_strong_idx];
            Users(ue_weak_idx).SC_idx = SC_index;
            Users(ue_weak_idx).if_S_NOMA = (~isempty(other_ue_set));
            Users(ue_weak_idx).partner_S_NOMA = other_ue_set;
            % beam vector
            angle_set_strong_ue = Users(ue_strong_idx).ang;
            angle_set_weak_ue = Users(ue_weak_idx).ang;
            ang_ = (angle_set_strong_ue + angle_set_weak_ue)/2;
            a_vector_temp = zeros(BS.Nt, 1);
            for i = 1:BS.Nt
                a_vector_temp(i) = exp(1i*sin(ang_ - 2*pi*(i-1)/BS.Nt));
            end
            Users(ue_weak_idx).beam_vector = a_vector_temp/norm(a_vector_temp);
            Users(ue_strong_idx).beam_vector = a_vector_temp/norm(a_vector_temp);
        end
        alpha = alpha + Users(ue_strong_idx).if_P_NOMA;
        beta = beta + Users(ue_strong_idx).if_S_NOMA;
    end % end of "for pair_idx = 1:3"
    if alpha > 0 && beta == 0
        Groups(group_idx).MA_mode = 1; %1: P-NOMA mode
    end
    if alpha == 0 && beta > 0
        Groups(group_idx).MA_mode = 2; %2: S-NOMA mode
    end
    if alpha >0 && beta > 0
        Groups(group_idx).MA_mode = 3; %3: Hybrid NOMA mode
    end
    Groups(group_idx).user_subset = setdiff(Groups(group_idx).user_subset, [0]);
end% end of "for group_idx = 1 : BS.num_of_SCs"
if sum(matched_indicator) == BS.num_of_UE
    for i = 1 : BS.num_of_UE
        if Users(i).SC_idx == 0
            fprintf('Check.\n');
        end
    end
end
end
% a_1N = BS.pow_allocate_coeff*rho*(Users(u).dist^(-3.76));% formula (12), page 2296, RHS
% beta = 0.01;
% C1 = (1 + BS.Nt*a_1N)*(1 + BS.Nt*a_1N*beta*(1/BS.pow_allocate_coeff - 1)); % Constant C1 in page 2299, LHS
% C2 = 1 + BS.Nt*a_1N + BS.Nt*a_1N*beta*(1/BS.pow_allocate_coeff - 1);
% eta = (C1 - C2)/(BS.Nt*a_1N*(C2*(beta + 1/BS.pow_allocate_coeff - 1) - C1*(beta/BS.pow_allocate_coeff)));
% %eta_max = 1/(BS.Rician_fading_factor*BS.Nt - rho*(Users(u).dist^(-3.76))/(1+BS.Rician_fading_factor));
% Users(u).delta_d = Users(u).dist*(eta^(-1/3.76) - 1);