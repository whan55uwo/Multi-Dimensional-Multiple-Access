function [flag_switch, sum_cost_new, utilization_cost_UE_array_new, temp_UE_coalitions_new] = ...
    Rotation_S3(ue_x_idx, ue_y_idx, ue_z_idx, step_size, temp_UE_coalitions, utilization_cost_UE_array)
%% INPUT:
% 1. Switched UE index, the new UE add inthe coalition
% 2. step_size = 1 OR 2
%% OUTPUT
% 1. FLAG = 1, if this switch is successed
% 2. Renew the sum_cost_new, if flag = 1
% 3. Renew matrix "temp_UE_coalitions", if flag = 1
global BS;
global Users;
%%
temp_UE_coalitions_new = [];
sum_cost_new = [];
utilization_cost_UE_array_new = [];
flag_switch = 0;
coalition_A_idx = find(temp_UE_coalitions(:, ue_x_idx) == 1);
coalition_B_idx = find(temp_UE_coalitions(:, ue_y_idx) == 1);
coalition_C_idx = find(temp_UE_coalitions(:, ue_z_idx) == 1);
if coalition_A_idx == coalition_B_idx || coalition_B_idx == coalition_C_idx ...
        || coalition_A_idx == coalition_C_idx
    return;
end
coalition_A_origin = find(temp_UE_coalitions(coalition_A_idx, :) == 1);
coalition_B_origin = find(temp_UE_coalitions(coalition_B_idx, :) == 1);
coalition_C_origin = find(temp_UE_coalitions(coalition_C_idx, :) == 1);
if step_size == 1
    coalition_A_new = [setdiff(coalition_A_origin, [ue_x_idx]), ue_z_idx];
    [flag1] = is_forbidden_coalition(coalition_A_new, ue_z_idx);
    if flag1 == 1
        return;
    end % Exclude forbidden pair
    coalition_B_new = [setdiff(coalition_B_origin, [ue_y_idx]), ue_x_idx];
    [flag2] = is_forbidden_coalition(coalition_B_new, ue_x_idx);
    if flag2 == 1
        return;
    end % Exclude forbidden pair
    coalition_C_new = [setdiff(coalition_C_origin, [ue_z_idx]), ue_y_idx];
    [flag3] = is_forbidden_coalition(coalition_C_new, ue_y_idx);
    if flag3 == 1
        return;
    end % Exclude forbidden pair
else
    coalition_A_new = [setdiff(coalition_A_origin, [ue_x_idx]), ue_y_idx];
    [flag1] = is_forbidden_coalition(coalition_A_new, ue_y_idx);
    if flag1 == 1
        return;
    end % Exclude forbidden pair
    coalition_B_new = [setdiff(coalition_B_origin, [ue_y_idx]), ue_z_idx];
    [flag2] = is_forbidden_coalition(coalition_B_new, ue_z_idx);
    if flag2 == 1
        return;
    end % Exclude forbidden pair
    coalition_C_new = [setdiff(coalition_C_origin, [ue_z_idx]), ue_x_idx];
    [flag3] = is_forbidden_coalition(coalition_C_new, ue_x_idx);
    if flag3 == 1
        return;
    end % Exclude forbidden pair
end
[cost_coalition_A_new] = g_cost_est_cal(coalition_A_new, 0);
[cost_coalition_B_new] = g_cost_est_cal(coalition_B_new, 0);
[cost_coalition_C_new] = g_cost_est_cal(coalition_C_new, 0);
utilization_cost_UE_array_new = utilization_cost_UE_array;
utilization_cost_UE_array_new(coalition_A_new) = cost_coalition_A_new;
utilization_cost_UE_array_new(coalition_B_new) = cost_coalition_B_new;
utilization_cost_UE_array_new(coalition_C_new) = cost_coalition_C_new;
sum_cost_new = sum(utilization_cost_UE_array_new.*([Users(1:BS.num_of_UE).cost_weight])');
if sum_cost_new < sum(utilization_cost_UE_array)
    flag_switch = 1;
    %% Renew matrix "temp_UE_coalitions"
    if step_size == 1
        temp_UE_coalitions_new = temp_UE_coalitions;
        temp_UE_coalitions_new(:, ue_x_idx) = 0;
        temp_UE_coalitions_new(coalition_B_idx, ue_x_idx) = 1;
        temp_UE_coalitions_new(:, ue_y_idx) = 0;
        temp_UE_coalitions_new(coalition_C_idx, ue_y_idx) = 1;
        temp_UE_coalitions_new(:, ue_z_idx) = 0;
        temp_UE_coalitions_new(coalition_A_idx, ue_z_idx) = 1;
    else
        temp_UE_coalitions_new = temp_UE_coalitions;
        temp_UE_coalitions_new(:, ue_x_idx) = 0;
        temp_UE_coalitions_new(coalition_C_idx, ue_x_idx) = 1;
        temp_UE_coalitions_new(:, ue_y_idx) = 0;
        temp_UE_coalitions_new(coalition_A_idx, ue_y_idx) = 1;
        temp_UE_coalitions_new(:, ue_z_idx) = 0;
        temp_UE_coalitions_new(coalition_B_idx, ue_z_idx) = 1;
    end % if step_size == 1
end
end