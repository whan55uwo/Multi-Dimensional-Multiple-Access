function [flag_switch, sum_cost_new, utilization_cost_UE_array_new, temp_UE_coalitions_new] =...
    Rotation_S2(ue_x_idx, ue_y_idx, temp_UE_coalitions, utilization_cost_UE_array)
%% INPUT:
% Switched UE index, the new UE add inthe coalition
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
ue_x_coalition_idx_origin = find(temp_UE_coalitions(:, ue_x_idx) == 1);
ue_y_coalition_idx_origin = find(temp_UE_coalitions(:, ue_y_idx) == 1);
if ue_x_coalition_idx_origin == ue_y_coalition_idx_origin
    return;
end
ue_x_coalition_origin = find(temp_UE_coalitions(ue_x_coalition_idx_origin, :) == 1);
ue_y_coalition_origin = find(temp_UE_coalitions(ue_y_coalition_idx_origin, :) == 1);
ue_x_coalition_new = [setdiff(ue_x_coalition_origin, [ue_x_idx]), ue_y_idx];
ue_y_coalition_new = [setdiff(ue_y_coalition_origin, [ue_y_idx]), ue_x_idx];
[flag1] = is_forbidden_coalition(ue_x_coalition_new, ue_y_idx);
[flag2] = is_forbidden_coalition(ue_y_coalition_new, ue_x_idx);
if flag1 == 1 || flag2 == 1
    return;
end % Exclude forbidden pair
[cost_coalition_x_new] = g_cost_est_cal(ue_x_coalition_new, 0);
[cost_coalition_y_new] = g_cost_est_cal(ue_y_coalition_new, 0);
utilization_cost_UE_array_new = utilization_cost_UE_array;
utilization_cost_UE_array_new(ue_x_coalition_new) = cost_coalition_x_new;
utilization_cost_UE_array_new(ue_y_coalition_new) = cost_coalition_y_new;
sum_cost_new = sum(utilization_cost_UE_array_new.*([Users(1:BS.num_of_UE).cost_weight])');
if sum_cost_new < sum(utilization_cost_UE_array)
    flag_switch = 1;
    % Renew matrix "temp_UE_coalitions"
    temp_UE_coalitions_new = temp_UE_coalitions;
    temp_UE_coalitions_new(:, ue_x_idx) = 0;
    temp_UE_coalitions_new(ue_y_coalition_idx_origin, ue_x_idx) = 1;
    temp_UE_coalitions_new(:, ue_y_idx) = 0;
    temp_UE_coalitions_new(ue_x_coalition_idx_origin, ue_y_idx) = 1;
end
end
