function [sum_cost, utilization_cost_UE_array, iteration_rotation_idx] = User_Coalition_Formation_Process()
global BS; % type: structure
%% Initialization for outer iteration
best_performance = inf;
best_UE_coalitions = [];
avg_u_cost_record_all = cell(10, 1);
iteration_all = 1;
for outer_loop_idx = 1 : 10
    %% Initialization for inner iteration
    temp_UE_coalitions = zeros(BS.num_of_SCs, BS.num_of_UE);
    UE_Coalition_Prference_List = cell(BS.num_of_UE, 1); % it is a K-size cell matrix
    UEs_unmatched = ones(BS.num_of_UE, 1); % 1: unmatched; and 0: matched
    %% Phase I ; establish a feasible coalition
    %Step 1
    rand_permutation_list = randperm(BS.num_of_UE);
    for Coalition_idx = 1 : BS.num_of_SCs
        % try to occpuy one idle SC
        loop_UE_idx = rand_permutation_list(Coalition_idx);
        temp_UE_coalitions(Coalition_idx, loop_UE_idx) = 1;
        UEs_unmatched(loop_UE_idx) = 0;
    end % end of loop_UE_idx
    %Step 2
    for loop_UE = 1 : BS.num_of_UE
        loop_UE_idx = rand_permutation_list(loop_UE);
        if UEs_unmatched(loop_UE_idx) == 0
            % means that this UE has been allocated one Coalition
            continue
        end
        UE_Coalition_Prference_List{loop_UE_idx, 1} = zeros(BS.num_of_SCs, 1);
        for Coalition_idx = 1 : BS.num_of_SCs
            % L_max
            UE_coalition_origin = find(temp_UE_coalitions(Coalition_idx, :) == 1);
            coalition_origin_size = length(UE_coalition_origin);
            if coalition_origin_size == BS.L_max
                UE_Coalition_Prference_List{loop_UE_idx, 1}(Coalition_idx,1) = Inf;
                continue;
            end
            % Exclude forbidden pair: Constraint C1
            UE_coalition_new = [UE_coalition_origin, loop_UE_idx];
            [flag] = is_forbidden_coalition(UE_coalition_new, loop_UE_idx);
            if  flag == 1
                UE_Coalition_Prference_List{loop_UE_idx, 1}(Coalition_idx,1) = Inf;
                continue;
            end
            % Calculate
            utilization_cost_array = g_cost_est_cal(UE_coalition_new, 0);
            UE_Coalition_Prference_List{loop_UE_idx, 1}(Coalition_idx,1) = ...
                utilization_cost_array(end);
        end % for loop_C = 1 : BS.num_of_SCs
        [min_value, C_index] =  min(UE_Coalition_Prference_List{loop_UE_idx, 1});
        temp_UE_coalitions(C_index, loop_UE_idx) = 1;
        UEs_unmatched(loop_UE_idx) = 0;
    end % end of loop_UE_idx
    % the output of this part
    %% Phase II: Rotation for optimality
    iteration_rotation_idx = 1;
    iteration_rotation_min = 100;
    iteration_rotation_max = 1000;
    weighted_sum_cost_record = zeros(iteration_rotation_max, 1);
    [sum_cost, utilization_cost_UE_array] = g_cost_est_all_ue_cal(temp_UE_coalitions, 0);
    cost_record_idx = 1;
    weighted_sum_cost_record(cost_record_idx) = sum_cost;
    continus_failed_num = 0;
    while 1
        % Rotation and Switch
        continus_failed_num = continus_failed_num + 1;
        rotation_cases_S2 = nchoosek(1 : BS.num_of_UE, 2);
        len_rotation_cases_S2 = size(rotation_cases_S2, 1);
        rotation_cases_S3 = nchoosek(1 : BS.num_of_UE, 3);
        len_rotation_cases_S3 = size(rotation_cases_S3, 1);
        if iteration_rotation_idx <= len_rotation_cases_S2
            %------------------Size =2-------------------------------
            rotation_case = rotation_cases_S2(iteration_rotation_idx, :);
            % Exclude the forbidden cases
            ue_x_idx = rotation_case(1);
            ue_y_idx = rotation_case(2);
            [flag_switch, sum_cost_new, utilization_cost_UE_array_new, temp_UE_coalitions_new] =...
                Rotation_S2(ue_x_idx, ue_y_idx, temp_UE_coalitions, utilization_cost_UE_array);
        else
            %------------------Size = 3---------------------------------
            temp_ = iteration_rotation_idx - len_rotation_cases_S2;
            S3_rotation_idx = ceil(temp_/2);
            rotation_case = rotation_cases_S3(S3_rotation_idx, :);
            ue_x_idx = rotation_case(1);
            ue_y_idx = rotation_case(2);
            ue_z_idx = rotation_case(3);
            if mod(temp_, 2) == 1
                [flag_switch, sum_cost_new, utilization_cost_UE_array_new, temp_UE_coalitions_new] = ...
                    Rotation_S3(ue_x_idx, ue_y_idx, ue_z_idx, 1, temp_UE_coalitions, utilization_cost_UE_array);
            else
                [flag_switch, sum_cost_new, utilization_cost_UE_array_new, temp_UE_coalitions_new] = ...
                    Rotation_S3(ue_x_idx, ue_y_idx, ue_z_idx, 2, temp_UE_coalitions, utilization_cost_UE_array);
            end
        end
        % determin if "sum_cost_new" is better than former cases
        if flag_switch == 1 && weighted_sum_cost_record(cost_record_idx) > sum_cost_new
            cost_record_idx = cost_record_idx + 1;
            weighted_sum_cost_record(cost_record_idx) = sum_cost_new;
            temp_UE_coalitions = temp_UE_coalitions_new;
            utilization_cost_UE_array = utilization_cost_UE_array_new;
            continus_failed_num = 0;
        end
        % Break Conditions
        if iteration_rotation_idx >= iteration_rotation_max || iteration_rotation_idx >= (len_rotation_cases_S2 + len_rotation_cases_S3)
            break;
        end
        if  (iteration_rotation_idx > iteration_rotation_min) && (continus_failed_num > 70)
            break
        end
        % Renew iteration_rotation_idx
        iteration_rotation_idx = iteration_rotation_idx + 1;
        iteration_all = iteration_all + 1;
    end % end of inner iteration
    avg_u_cost_record_all{outer_loop_idx, 1} = weighted_sum_cost_record;
    if weighted_sum_cost_record(cost_record_idx) < best_performance
        best_UE_coalitions = temp_UE_coalitions;
        best_performance = weighted_sum_cost_record(cost_record_idx);
    end 
end % end of outer iteration
%% Output
[sum_cost, utilization_cost_UE_array] = g_cost_est_all_ue_cal(best_UE_coalitions, 1);
% DETERMINE THE VALUE OF GLOBAL VARIABLES
end
