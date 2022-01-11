function [SC_allocation] = User_Grouping_SC_Allocation()

global BS; % type: structure
global Users;
global Beamspaces;
global Groups;

%% Initialization
temp_SC_allocation = zeros(BS.num_of_SCs, BS.num_of_UE); % the output of this part
UEs_unmatched = ones(BS.num_of_UE, 1); % 1: unmatched; and 0: matched
SC_occpuied_flag = zeros(BS.num_of_SCs, 1); % see if this SC has been occpied, 0-idle, 1-occpuied
UE_SC_Prference_List = cell(BS.num_of_UE, 2); % it is a K-size cell matrix
[Users(1:BS.num_of_UE).channel_gain]    = deal(0); % which beamspace this UE belongs to
for loop_UE_idx = 1 : BS.num_of_UE
    % caluculation of channel conditions: w_c and h
    UE_SC_Prference_List{loop_UE_idx, 1} = zeros(BS.num_of_SCs, 1);
    for loop_SC = 1 : BS.num_of_SCs
        h = Users(loop_UE_idx).H_matrix(:, loop_SC);
        UE_SC_Prference_List{loop_UE_idx, 1}(loop_SC,1) = norm(h);
    end % end of loop_SC
    % gloabl variable
    Users(loop_UE_idx).channel_gain = UE_SC_Prference_List{loop_UE_idx, 1};
    % calculate the preference list
    [B, I] = sort(UE_SC_Prference_List{loop_UE_idx, 1}, 'descend');
    UE_SC_Prference_List{loop_UE_idx, 2} = I;
end % end of loop_UE_idx
%% Phase 1: Obtain a  Feasiable Matching
%Step 1
rand_permutation_list = randperm(BS.num_of_UE);
for loop_UE = 1 : BS.num_of_UE
    % try to occpuy one idle SC with best channel codntion
    loop_UE_idx = rand_permutation_list(loop_UE);
    if sum(SC_occpuied_flag) <  BS.num_of_SCs % if there are still some idle SCs
        idle_SC_set = find(SC_occpuied_flag == 0);
        channel_condition_set = UE_SC_Prference_List{loop_UE_idx, 1}(SC_occpuied_flag == 0,1);
        [Y, index] = max(channel_condition_set);
        temp_SC_allocation(idle_SC_set(index), loop_UE_idx) = 1; % choose the best channel conditon SC in the idle SC set
        SC_occpuied_flag(idle_SC_set(index)) = 1;
        UEs_unmatched(loop_UE_idx) = 0; % means that this UE has been allocated one SC
    else
        continue;
    end % end if sum(...)
end % end of loop_UE_idx
%Step 2: Accepted or Rejected Process
rand_permutation_list = randperm(BS.num_of_UE);
for loop_UE = 1 : BS.num_of_UE
    loop_UE_idx = rand_permutation_list(loop_UE);
    if UEs_unmatched(loop_UE_idx) == 0 % means that this UE has been allocated one SC
        continue
    end
    temp_SC_allocation(:, loop_UE_idx) = 0;
    SC_index = UE_SC_Prference_List{loop_UE_idx, 2}(1);
    UE_group_origin = find(temp_SC_allocation(SC_index, :) == 1);
    group_size = length(UE_group_origin);
    % Exclude forbidden pair
    UE_group = [UE_group_origin, loop_UE_idx];
    % Determine Accepted or Rjected
    if group_size < BS.L_max && rand > 0.5
        % Accepted
        UEs_unmatched(loop_UE_idx) = 0;
        temp_SC_allocation(SC_index, loop_UE_idx) = 1; % choose the best channel conditon SC in the idle SC set
    else
        % Rejected
        UE_group_candidates = nchoosek(UE_group, BS.L_max);
        len_ = size(UE_group_candidates, 1);
        g_cost_group = zeros(len_, 1);
        for  loop_UE_group_idx = 1 : len_
            UE_group_candidate = UE_group_candidates(loop_UE_group_idx, :);
            g_cost_group(loop_UE_group_idx, 1) = sum(g_cost_cal_est(UE_group_candidate, SC_index));
        end % for  loop_UE_group_idx = 1 : len_
        [Y, index] = min(g_cost_group);
        UE_group_new = UE_group_candidates(index, :);
        Rejected_UE = setdiff(UE_group_origin, UE_group_new);
        if isempty(Rejected_UE)
            UEs_unmatched(loop_UE_idx) = 1;
        else
            UEs_unmatched(Rejected_UE) = 1; % Rejected_UE is unmatched
            UEs_unmatched(loop_UE_idx) = 0; % this UE is accpeted on this SC
            temp_SC_allocation(:, Rejected_UE) = 0; % Withdraw the allocation of this SC for Rejected_UE
            temp_SC_allocation(SC_index, loop_UE_idx) = 1;
        end % if isempty(Rejected_UE)
    end % if group_size < BS.L_max
end % end of loop_UE_idx
%Step 3
rand_permutation_list = randperm(BS.num_of_UE);
for loop_UE = 1 : BS.num_of_UE
    loop_UE_idx = rand_permutation_list(loop_UE);
    if UEs_unmatched(loop_UE_idx) == 0
        continue
    end
    beam_space_idx = Users(loop_UE_idx).beam_idx;
    beam_user_set = Beamspaces(beam_space_idx).users_set;
    temp_SC_allocation(:, loop_UE_idx) = 0;
    for i = 1 : BS.num_of_SCs
        SC_index = UE_SC_Prference_List{loop_UE_idx, 2}(i);
        UE_group_origin = find(temp_SC_allocation(SC_index, :) == 1);
        group_size = length(UE_group_origin);
        if group_size == BS.L_max
            continue;
        end
        UE_group = [UE_group_origin, loop_UE_idx];
        if group_size < BS.L_max
            % Accepted
            UEs_unmatched(loop_UE_idx) = 0;
            temp_SC_allocation(SC_index, loop_UE_idx) = 1; % choose the suitable SC for this unmatched UE
            break;
        end
    end % for i = 1 : BS.num_of_SCs
end
%% Phase 2: Rotation
iteration_rotation_idx = 1;
iteration_rotation_min = 100;
iteration_rotation_max = 400;
continus_failed_num = 0;
while (iteration_rotation_idx<=iteration_rotation_max)
    % Rotation and Switch
    continus_failed_num = continus_failed_num + 1;
    % Step 1: Rnadomly choose two SCs
    while 1
        x = randi(BS.num_of_SCs);
        y = randi(BS.num_of_SCs);
        if x ~= y
            Rotation_SC_set2 = [x, y];
            break;
        end % end of if x ~= y
    end % end of while
    % Step 2: Swtich UE pair
    % Slect two UEs
    SC_index_x = Rotation_SC_set2(1);
    UE_group_x_origin = find(temp_SC_allocation(SC_index_x, :) == 1);
    UE_index_x = randi(length(UE_group_x_origin)); % the index in the UE_GROUP_X
    UE_idx_x = UE_group_x_origin(UE_index_x);
    SC_index_y = Rotation_SC_set2(2);
    UE_group_y_origin = find(temp_SC_allocation(SC_index_y, :) == 1);
    UE_idnex_y = randi(length(UE_group_y_origin)); % the index in the UE_GROUP_Y
    UE_idx_y = UE_group_y_origin(UE_idnex_y);
    % Switch
    UE_group_x_new = UE_group_x_origin;
    UE_group_x_new(UE_index_x) = UE_idx_y;
    UE_group_y_new = UE_group_y_origin;
    UE_group_y_new(UE_idnex_y) = UE_idx_x;
    % compare g_cost: output: variable if_switch, 1-switch, 0-deny
    g_cost_SCx_origin = sum(g_cost_cal_est(UE_group_x_origin, SC_index_x));
    g_cost_SCy_origin = sum(g_cost_cal_est(UE_group_y_origin, SC_index_y));
    g_cost_SCx_new = sum(g_cost_cal_est(UE_group_x_new, SC_index_x));
    g_cost_SCy_new = sum(g_cost_cal_est(UE_group_y_new, SC_index_y));
    if  (g_cost_SCx_new+g_cost_SCy_new) < (g_cost_SCx_origin+g_cost_SCy_origin)% swith user in each SC UE set
        temp_SC_allocation(:, UE_idx_x) = 0;
        temp_SC_allocation(SC_index_y, UE_idx_x) = 1;
        temp_SC_allocation(:, UE_idx_y) = 0;
        temp_SC_allocation(SC_index_x, UE_idx_y) = 1;
        continus_failed_num = 0;
    end
    % terminate condition
    if  (iteration_rotation_idx > iteration_rotation_min) && (continus_failed_num > 10)
        break
    end
    iteration_rotation_idx = iteration_rotation_idx + 1;
end % end of while
%% Output
SC_allocation = temp_SC_allocation; % the output of this part
for loop_SC = 1 : BS.num_of_SCs
    UE_group = find(SC_allocation(loop_SC, :) == 1);
    UE_group_size = length(UE_group);
    if UE_group_size > 0
        Groups(loop_SC).user_group = UE_group;
        if UE_group_size > 1
            for i = 1 : UE_group_size
                UE_i_idx = UE_group(i);
                Users(UE_i_idx).partner = setdiff(UE_group, UE_i_idx);
                Users(UE_i_idx).SC_idx = loop_SC;
            end % end of "for i = 1 : UE_group_size"
        else
            UE_i_idx = UE_group(1);
            Users(UE_i_idx).partner = [];
            Users(UE_i_idx).SC_idx = loop_SC;
        end % end of "if UE_group_size > 1"
    else
        fprintf('SC index = %d is idel status!\n', loop_SC);
        Groups(loop_SC).user_group = [];
    end
end % end of loop_SC
%% Check
for loop_UE = 1 : BS.num_of_UE
    if  sum(SC_allocation(:,loop_UE)) == 0
        fprintf('error!\n');
    end
end % end of loop_UE
end