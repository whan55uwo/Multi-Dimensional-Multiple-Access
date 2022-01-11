function [iteration_idx] = SC_power_allocation()
global BS; % type: structure
global Users;
global Coalitions;
%% Initialization
ZF_precoder_at_SC();
%Constant Values D_k and G_{i,k,m}
D_0 = zeros(BS.num_of_UE, 1);
G_0 = zeros(BS.num_of_UE, BS.num_of_UE, BS.num_of_SCs);
for k = 1 : BS.num_of_UE
    D_0(k) = (1/log(2))*BS.sys_SC_bandwidth*BS.a_k;
    for i = 1 : BS.num_of_UE
        for m = 1 : BS.num_of_SCs
            temp_numerator = norm((Users(i).H_matrix(:, m))'*(Users(k).beam_vector_matrix(:,m)))^2;
            temp_denominator = norm((Users(i).H_matrix(:, m))'*(Users(i).beam_vector_matrix(:,m)))^2;
            G_0(i,k,m) = (1/log(2))*BS.sys_SC_bandwidth*BS.a_k*(temp_numerator/temp_denominator);
        end % loop SC
    end % loop_i
end % loop_k
BS.D0 = D_0; % make matrix D_0 GLOBAL VARIABLE
BS.G0 = G_0; % make matrix G_0 GLOBAL VARIABLE
%% Iteration Processing
% Largrange Mulitpler
iteration_idx = 1;
Mu_vector = zeros(BS.num_of_UE, 1); % realted to right side of C4
Nu = 1e-39;
Nu_max = 20*mean(D_0)*(1/BS.rate_max);
Nu_min = 0;
outer_iteration_index = 1; % the index of outer iteration
max_outer_iter_num = 24;
Total_power_record = zeros(max_outer_iter_num, 1);
while 1 % Iteration related to Nu by bi-searching method
    %% Inialitzation of the inner iteration
    inner_iteration_index = 1; % the index of inner iteration
    max_inner_iter_num = 200;
    min_inner_iter_num = 10; % the objective function has no significant decrease within min_iter_num iterations
    epsilon = 1e-8;
    window_size_grad = 20;
    lambda_grad = 1-1/window_size_grad;
    window_size_deacy = 20;
    lambda_deacy = 1-1/window_size_deacy;
    J_fun_record = zeros(max_inner_iter_num,1); % the record window of current value of obejetive function
    J_fun = Inf; % the record of best value of obejetive function
    J_indiviual = zeros(BS.num_of_SCs, BS.num_of_SCs); %row: Colaition index, column: SC index
    temp_power_vector = ones(BS.num_of_UE, 1)*(BS.tx_power/BS.num_of_UE);
    temp_SINR_vector = zeros(BS.num_of_UE, 1);
    temp_SINR_SIC_vector = zeros(BS.num_of_UE, 1);
    temp_Rate_vector = zeros(BS.num_of_UE, 1);
    Eta_vector = zeros(BS.num_of_UE, 1); % realted to left side of C4
    decay_ave_grad_pre = zeros(BS.num_of_UE, 1); %  the decaying average of gradient at previous iteration
    decay_ave_step_size_pre = zeros(BS.num_of_UE, 1); %  the decaying average of step size at previous iteration
    EMS_step_size_pre = zeros(BS.num_of_UE, 1); %  the root mean squared (RMS) error criterion of step size
    while 1
        %% Solve Problem (13a) based on Lemma 1
        power_mat = cell(BS.num_of_SCs, BS.num_of_SCs); %row: Colaition index, column: SC index
        SINR_mat = cell(BS.num_of_SCs, BS.num_of_SCs); %row: Colaition index, column: SC index
        SINR_SIC_mat = cell(BS.num_of_SCs, BS.num_of_SCs); %row: Colaition index, column: SC index
        Rate_mat = cell(BS.num_of_SCs, BS.num_of_SCs); %row: Colaition index, column: SC index
        Rate_low_bound_mat = cell(BS.num_of_SCs, BS.num_of_SCs); %row: Colaition index, column: SC index
        for Coalition_idx = 1 : BS.num_of_SCs % index of coalition
            for SC_idx = 1 : BS.num_of_SCs % index of SC
                [power_, J_nm_value, Rate, Rate_low_bound, SINR_, SINR_SIC_, iteration_] = fix_point_equation(Eta_vector,Mu_vector,Nu,Coalition_idx,SC_idx);
                J_indiviual(Coalition_idx, SC_idx) = J_nm_value;
                SINR_mat{Coalition_idx, SC_idx} = SINR_;
                SINR_SIC_mat{Coalition_idx, SC_idx} = SINR_SIC_;
                power_mat{Coalition_idx, SC_idx} = power_;
                Rate_mat{Coalition_idx, SC_idx} = Rate;
                Rate_low_bound_mat{Coalition_idx, SC_idx} = Rate_low_bound;
            end % end of "for SC_idx = 1 : BS.num_of_SCs"
        end % end of "for Coalition_idx = 1 : BS.num_of_SCs"
        %% SC assignment
        % biparte matching
        if outer_iteration_index == 1
            [Temp_SC_allocation_mat,J_value,iteration_num_]=bipartite_graph_match(-J_indiviual); %row: Colaition index, column: SC index
            J_fun_pre = J_fun;
            J_fun = -J_value;
            J_fun_record(inner_iteration_index) = -J_value;
        else
            J_fun_pre = J_fun;
            J_fun = sum(sum(J_indiviual(Temp_SC_allocation_mat == 1)));
        end
        % determine power allocation
        for Coalition_idx = 1 : BS.num_of_SCs % index of coalition
            index = find(Temp_SC_allocation_mat(Coalition_idx,:) == 1);
            index = index(1);
            UE_set = Coalitions(Coalition_idx).user_subset;
            UE_set_size = length(UE_set);
            power_vector = power_mat{Coalition_idx, index};
            Rate_vecotor = Rate_mat{Coalition_idx, index};
            SINR_vector = SINR_mat{Coalition_idx, index};
            SINR_SIC_vector = SINR_SIC_mat{Coalition_idx, index};
            for loop_UE = 1 : UE_set_size
                UE_idx = UE_set(loop_UE);
                temp_power_vector(UE_idx) = power_vector(loop_UE);
                temp_Rate_vector(UE_idx) = Rate_vecotor(loop_UE);
                temp_SINR_vector(UE_idx) = SINR_vector(loop_UE);
                temp_SINR_SIC_vector(UE_idx) = SINR_SIC_vector(loop_UE);
            end % loop_UE = 1 : UE_set_size
        end % end of "for Coalition_idx = 1 : BS.num_of_SCs"
        %% Update of Larange Multipliers Eta Vector based on Adadelta Alogorithm
        grad_Eta_vector = (temp_Rate_vector - BS.rate_min)/BS.rate_min;
        % calculate the decaying average of gradient with respect to the objective function
        decay_ave_grad = lambda_grad*decay_ave_grad_pre +  (1-lambda_grad)*(grad_Eta_vector).^2;
        decay_ave_grad_pre =  decay_ave_grad;
        EMS_grad = sqrt(decay_ave_grad + epsilon);
        % calculate the step size Eta
        if inner_iteration_index <= max(window_size_grad, window_size_deacy)
            step_size = -0.01*(1/BS.rate_max)*grad_Eta_vector;
            step_size = max(step_size, -0.1/BS.rate_max);
            step_size = min(step_size,  0.2/BS.rate_max);
        else
            step_size = -(EMS_step_size_pre./EMS_grad).*grad_Eta_vector;
            step_size = max(step_size, -0.1/BS.rate_max);
            step_size = max(step_size, 0.2/BS.rate_max);
        end % end of if inner_iteration_index
        % calculate the RMS of stepsize
        decay_ave_step_size =  lambda_deacy*decay_ave_step_size_pre + ...
            (1-lambda_deacy)*(step_size).^2;
        decay_ave_step_size_pre =  decay_ave_step_size;
        EMS_step_size_pre =  sqrt(decay_ave_step_size + epsilon);
        % Renew Eta Vector
        Eta_vector = max(Eta_vector + step_size, 0);
        %% Termination condition of inner iteration
        if inner_iteration_index == max_inner_iter_num
            break; % end the programme
        end
        if outer_iteration_index == 1 && inner_iteration_index >= 60
            break;
        end
        if (inner_iteration_index > min_inner_iter_num) && (inner_iteration_index <= max_inner_iter_num/2)
            if  abs(J_fun_pre - J_fun) <= 0.05*abs(J_fun_pre) && sum(grad_Eta_vector >= -0.05) == BS.num_of_UE
                break; % end the programme
            end
        end
        if (inner_iteration_index > max_inner_iter_num/2) && (inner_iteration_index <= max_inner_iter_num*(2/3))
            if  abs(J_fun_pre - J_fun) <= 0.1*abs(J_fun_pre) && sum(grad_Eta_vector >= -0.1) == BS.num_of_UE
                break; % end the programme
            end
        end
        if (inner_iteration_index > max_inner_iter_num*(2/3))
            if  abs(J_fun_pre - J_fun) <= 0.1*abs(J_fun_pre) && sum(grad_Eta_vector >= -0.2) == BS.num_of_UE
                break; % end the programme
            end
        end
        inner_iteration_index = inner_iteration_index + 1; % Renew  inner_iteration_index
        fprintf("         Outer iteration index = %d, Inner iteration index = %d.\n", outer_iteration_index, inner_iteration_index);
        iteration_idx = iteration_idx + 1;
    end % end of inner while
    %% update Nu by bi-section search
    total_power = sum(temp_power_vector);
    Total_power_record(outer_iteration_index) = total_power;
    if  total_power <= BS.tx_power && outer_iteration_index == 1
        break;
    end
    if outer_iteration_index == 1
        Nu = (Nu_max + Nu_min)/2;
        outer_iteration_index = outer_iteration_index + 1; % the index of outer iteration
        continue;
    end
    Nu = (Nu_max + Nu_min)/2;
    if total_power > BS.tx_power
        Nu_min = Nu;
    end
    if total_power < BS.tx_power
        Nu_max = Nu;
    end
    if outer_iteration_index == max_outer_iter_num || ...
            abs(total_power - BS.tx_power) < 1e-1*BS.tx_power ||...
            abs(total_power - Total_power_record(outer_iteration_index-1))*(outer_iteration_index>15)...
                    < 1e-3*Total_power_record(outer_iteration_index-1)*(outer_iteration_index>15)
        break;
    end
    outer_iteration_index = outer_iteration_index + 1; % the index of outer iteration
end % end of outer while
%% OUTPUT
for Coalition_idx = 1 : BS.num_of_SCs % index of coalition
    UE_set = Coalitions(Coalition_idx).user_subset;
    UE_set_size = length(UE_set);
    index_ = find(Temp_SC_allocation_mat(Coalition_idx,:) == 1);
    SC_index = index_(1);
    Coalitions(Coalition_idx).SC_idx = SC_index;
    for loop_UE = 1 : UE_set_size
        UE_idx = UE_set(loop_UE);
        Users(UE_idx).SC_idx = SC_index;
        Users(UE_idx).SINR = temp_SINR_vector(UE_idx);
        Users(UE_idx).SINR_SIC = temp_SINR_SIC_vector(UE_idx);
        Users(UE_idx).SINR_dB = 10*log10(temp_SINR_vector(UE_idx));
        Users(UE_idx).rate =  temp_Rate_vector(UE_idx);
        Users(UE_idx).rate_norm =  temp_Rate_vector(UE_idx)/BS.rate_max;
        Users(UE_idx).pow = temp_power_vector(UE_idx);
        Users(UE_idx).u_fun = Users(UE_idx).rate/BS.rate_max - Users(UE_idx).cost_weight*g_cost_cal(UE_idx, SC_index, Users(UE_idx).SINR_SIC);
    end % end of loop_UE = 1 : UE_set_size
end % end of "Coalition_idx = 1 : BS.num_of_SCs"
clc;
end % end of function
