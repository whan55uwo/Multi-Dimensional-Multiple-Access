function [sum_cost, utilization_cost_UE_array, iteration_rotation_idx] =  Many_to_one_matching_Yannan_TWC()

global BS; % type: structure
global Users;
global Beamspaces;

%% Initialization
temp_UE_coalitions = zeros(BS.num_of_SCs, BS.num_of_UE);
Coalition_empty_flag = zeros(BS.num_of_SCs, 1); % see if this SC has been occpied, 0-empty, 1-occpuied
UEs_unmatched = ones(BS.num_of_UE, 1); % 1: unmatched; and 0: matched
UE_Coalition_Prference_List = cell(BS.num_of_UE, 2); % it is a K-size cell matrix
%% Phase 1: Obtain a  Feasiable Matching
iteration_rotation_idx =  1;
rand_permutation_list = randperm(BS.num_of_UE);
for loop_UE = 1 : BS.num_of_UE
    % try to occpuy one empty coalition with best channel codntion
    loop_UE_idx = rand_permutation_list(loop_UE);
    UE_Coalition_Prference_List{loop_UE_idx, 1} = rand(BS.num_of_SCs, 1);
    if sum(Coalition_empty_flag) <  BS.num_of_SCs % if there are still some empty coalition
        empty_coalition_set = find(Coalition_empty_flag == 0);
        random_set =UE_Coalition_Prference_List{loop_UE_idx, 1}(Coalition_empty_flag == 0, 1);
        [Y, index] = max(random_set);
        temp_UE_coalitions(empty_coalition_set(index), loop_UE_idx) = 1;
        % random choose a UE coalition
        Coalition_empty_flag(empty_coalition_set(index)) = 1;
        UEs_unmatched(loop_UE_idx) = 0; % means that this UE has been added into this coalition
    else
        continue;
    end % end if sum(...)
end % end of loop_UE_idx
%% Phase 2: Accepted or Rejected Process
for loop_UE = 1 : BS.num_of_UE
    [Coalition_Prference_List, C_index] = cal_coalition_preference(loop_UE, temp_UE_coalitions);
    UE_Coalition_Prference_List{loop_UE, 1} = Coalition_Prference_List;
    UE_Coalition_Prference_List{loop_UE, 2} = C_index ;
end
% iteration
rand_permutation_list = randperm(BS.num_of_UE);
while 1
    for loop_UE = 1 : BS.num_of_UE
        loop_UE_idx = rand_permutation_list(loop_UE);
        if UEs_unmatched(loop_UE_idx) == 0 % means that this UE has been allocated one SC
            continue
        end
        iteration_rotation_idx =  iteration_rotation_idx +  1;
        temp_UE_coalitions(:, loop_UE_idx) = 0;
        beam_space_idx = Users(loop_UE_idx).beam_idx;
        beam_user_set = Beamspaces(beam_space_idx).users_set;
        coalition_index = UE_Coalition_Prference_List{loop_UE_idx, 2}(1);
        UE_coalition_origin = find(temp_UE_coalitions(coalition_index, :) == 1);
        coalition_size = length(UE_coalition_origin);
        UE_coalition_new = [UE_coalition_origin, loop_UE_idx];
        % Exclude forbidden pair
        [flag] = is_forbidden_coalition(UE_coalition_origin, loop_UE_idx);
        if  flag == 1
            continue;
        end
        % Determine Accepted or Rjected
        if coalition_size < BS.L_max
            % Accepted
            UEs_unmatched(loop_UE_idx) = 0;
            temp_UE_coalitions(coalition_index, loop_UE_idx) = 1; % choose the best channel conditon SC in the idle SC set
        else
            % Rejected
            UE_coalition_candidates = nchoosek(UE_coalition_new, BS.L_max);
            len_ = size(UE_coalition_candidates, 1);
            g_cost_coalition = zeros(len_, 1);
            for  loop_coalition_idx = 1 : len_
                UE_coalition_candidate = UE_coalition_candidates(loop_coalition_idx, :);
                set_ = intersect(UE_coalition_candidate, beam_user_set);
                if length(set_) > 2
                    g_cost_coalition(loop_coalition_idx, 1) = 10000; % this is a forbbien pair
                elseif  length(set_) == 2 && Users(loop_UE_idx).if_SIC_capability == 0% Power-domain NOMA
                    % if this UE does not have SIC capability
                    distance = Users(loop_UE_idx).dist;
                    pair_UE_idx = setdiff(set_, [loop_UE_idx]);
                    distance_ = Users(pair_UE_idx).dist;
                    if distance <= distance_
                        g_cost_coalition(loop_coalition_idx, 1) = 10000; % this is a forbbien pair
                    end
                    if distance_ <= distance && Users(pair_UE_idx).if_SIC_capability == 0
                        % if partner has no SIC capability
                        g_cost_coalition(loop_coalition_idx, 1) = 10000; % this is a forbbien pair
                        continue;
                    end
                else
                    g_cost_coalition(loop_coalition_idx, 1) = sum(g_cost_est_cal(UE_coalition_candidate, 0));
                end
            end % for  loop_UE_group_idx = 1 : len_
            [Y, index] = min(g_cost_coalition);
            UE_coalition_new = UE_coalition_candidates(index, :);
            Rejected_UE = setdiff(UE_coalition_origin, UE_coalition_new);
            if isempty(Rejected_UE)
                UEs_unmatched(loop_UE_idx) = 1;
            else
                UEs_unmatched(Rejected_UE) = 1; % Rejected_UE is unmatched
                % Re-calculate Preference List: Update the preference list of the rejected UEs
                [Coalition_Prference_List, C_index] = cal_coalition_preference(Rejected_UE, temp_UE_coalitions);
                UE_Coalition_Prference_List{Rejected_UE, 1} = Coalition_Prference_List;
                UE_Coalition_Prference_List{Rejected_UE, 2} = C_index ;
                UEs_unmatched(loop_UE_idx) = 0; % this UE is accpeted on this SC
                temp_UE_coalitions(:, Rejected_UE) = 0; % Withdraw the allocation of this SC for Rejected_UE
                temp_UE_coalitions(coalition_index, loop_UE_idx) = 1;
            end % if isempty(Rejected_UE)
        end % if group_size < BS.L_max
    end % end of loop_UE_idx
    if sum(UEs_unmatched) == 0
        break;
    end
end % end of while
%% Output
best_UE_coalitions = temp_UE_coalitions;
[sum_cost, utilization_cost_UE_array] = g_cost_est_all_ue_cal(best_UE_coalitions, 1);
end



function [Coalition_Prference_List, C_index] = cal_coalition_preference(UE_idx, temp_UE_coalitions)
%% inistalization
global BS; % type: structure
Coalition_Prference_List = zeros(BS.num_of_SCs, 1);
%% cal
for Coalition_idx = 1 : BS.num_of_SCs
    % L_max
    UE_coalition_origin = find(temp_UE_coalitions(Coalition_idx, :) == 1);
    coalition_origin_size = length(UE_coalition_origin);
    if coalition_origin_size == BS.L_max
        Coalition_Prference_List(Coalition_idx,1) = Inf;
        continue;
    end
    % Exclude forbidden pair: Constraint C1
    UE_coalition_new = [UE_coalition_origin, UE_idx];
    [flag] = is_forbidden_coalition(UE_coalition_new, UE_idx);
    if  flag == 1
        Coalition_Prference_List(Coalition_idx,1) = Inf;
        continue;
    end
    % Calculate
    utilization_cost_array = g_cost_est_cal(UE_coalition_new, 0);
    Coalition_Prference_List(Coalition_idx,1) = ...
        utilization_cost_array(end);
end % for loop_C = 1 : BS.num_of_SCs
[min_value, C_index] =  min(Coalition_Prference_List);
end