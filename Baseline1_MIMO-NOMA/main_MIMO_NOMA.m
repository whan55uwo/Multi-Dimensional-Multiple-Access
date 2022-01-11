clc;
clear -global; clearvars;
close all; dbstop if error; digits(64);  % 64-bit operation precision, the default

% MIMO-NOMA
global BS; % type: structure
global Users;
global Groups;

%% 1.Setting parameters
CASE_IDX = 3;
BS.cell_radius = 250;        % m
BS.sys_SC_bandwidth = 2*10^6;  % Hz
BS.tx_power_dBm = 33;         % dBm, TBD
BS.noise_density_dB = -140;       % dBm/Hz, TBD
BS.Nr = 1; % number of recieve antennas at UE
BS.Nt = 64; % number of transmit antennas
BS.num_of_UE   = 25; % the number of users in the cell, To Be Detemined (TBD)
BS.num_of_SCs = 12;  % the number of SCs
BS.L_max = 3; % the maximum number of UEs can share one SC
BS.N_max = 2; % the maximum number of UEs can share one SC
BS.Rician_fading_factor = 9; %  uncorrelated Rician fading
BS.num_beam = 6; % Number of beamspaces in the cell, (TBD) not being used
BS.sinr_max = 10^(20/10); % not in dB
BS.rate_max = 0.5*BS.sys_SC_bandwidth*log2(1+BS.sinr_max); % in bits
BS.sinr_th = 10^(6/10); % not in dB
BS.rate_min = 0.5*BS.sys_SC_bandwidth*log2(1+BS.sinr_th); % in bits
BS.a_k = BS.sinr_th/(1+BS.sinr_th); % formula (11a)
BS.b_k = log2(1+BS.sinr_th) - BS.a_k*log2(BS.sinr_th); % formula (11b)
BS.constant_G = 0.5; % the cofficient in equation (6)
BS.rho1 = 0.24; % the cofficient in equation (6) realted to power domain NOMA
BS.rho2 = 0.5; % the cofficient in equation (6) related to SD-NOMA
BS.ue_SIC_ratio = 0.9;
% Calculate some parameters
BS.tx_power      = 10^(BS.tx_power_dBm / 10);  % in mW
BS.noise_power   = BS.sys_SC_bandwidth*10^(BS.noise_density_dB/10);  % in mW
SNR_dB=10*log10(BS.tx_power/BS.noise_power); % TBD,
drop_num = 100; % run this simulation for "drop_num" times
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('%% **************    Simulation Start!!! 1~%d drops.   ******************* %%\n', drop_num);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for loop_drop = 1:drop_num
    fprintf('%**************   loop_drop = %d, Begin!  ******************* %\n', loop_drop);
    BS.transmit_pow = deal(0);
    
    [Users(1:BS.num_of_UE).ang] = deal(0);
    [Users(1:BS.num_of_UE).a_vector] = deal(0);
    [Users(1:BS.num_of_UE).H_matrix] = deal(0);
    [Users(1:BS.num_of_UE).beam_vector] = deal(0);
    [Users(1:BS.num_of_UE).dist]  = deal(0);
    [Users(1:BS.num_of_UE).coor] = deal(0);
    [Users(1:BS.num_of_UE).pathloss]   = deal(0);
    [Users(1:BS.num_of_UE).candidates] = deal(0);
    [Users(1:BS.num_of_UE).if_SIC_capability] = deal(0); % hardware constraints: 1 MEANS HAVE SIC capability; 0, otherwises
    [Users(1:BS.num_of_UE).partner]    = deal(0); % the set of users share one SC with this UE
    [Users(1:BS.num_of_UE).beam_idx]    = deal(0); % which beamspace this UE belongs to
    [Users(1:BS.num_of_UE).if_P_NOMA]    = deal(0); % if this UE is using PD_NOMA, pls read fun "g_cost_cal"
    [Users(1:BS.num_of_UE).if_S_NOMA]    = deal(0); % if this UE is using SD_NOMA, pls read fun "g_cost_cal"
    [Users(1:BS.num_of_UE).if_Strong_UE_P_NOMA]    = deal(0); % if this UE is the "STRONG" UE in the PD_NOMA pair, see fun "g_cost_cal"
    [Users(1:BS.num_of_UE).partner_P_NOMA]    = deal(0); % if this UE is using PD_NOMA, its partner UE index, see fun "g_cost_cal"
    [Users(1:BS.num_of_UE).partner_S_NOMA]    = deal(0); % if this UE is using SD_NOMA, its partner UE index, see fun "g_cost_cal"
    [Users(1:BS.num_of_UE).PD_ND]    = deal(0); % the power_domain (PD) ND of UE, see fun "g_cost_cal"
    [Users(1:BS.num_of_UE).SD_ND]    = deal(0); % the spatial_domain (SD) ND of UE, see fun "g_cost_cal"
    [Users(1:BS.num_of_UE).utilization_cost]    = deal(0); % the see fun "g_cost_cal"
    [Users(1:BS.num_of_UE).SC_idx]    = deal(0); % the assigned SC to this UE, see fun "g_cost_cal"
    [Users(1:BS.num_of_UE).SINR]    = deal(0); % the data rate of UE, see fun "g_cost_cal"
    [Users(1:BS.num_of_UE).SINR_SIC]    = deal(0); % The SINR experienced by the SIC receiver
    [Users(1:BS.num_of_UE).SINR_dB]    = deal(0);
    [Users(1:BS.num_of_UE).rate]    = deal(0);
    [Users(1:BS.num_of_UE).rate_norm]    = deal(0);
    [Users(1:BS.num_of_UE).u_fun]    = deal(0);
    [Users(1:BS.num_of_UE).pow]    = deal(0);
    
    [Users(1:BS.num_of_UE).cost_weight] = deal(0); % objective funtion
    
    [Beamspaces(1:BS.num_beam).users_set] = deal(0); % the user set assoicated with beamspace
    [Beamspaces(1:BS.num_beam).beam_theta] = deal(0); % the angle of each beamspace
    
    sum_rate = deal(0);
    
    %% 3. randomly generate UEs and channel modeling
    for u = 1 : BS.num_of_UE
        a = rand()*2*pi;  % randomly generate angle in(0,2*pi)
        if (a > pi/3 && a < 2*pi/3) || (a > 4*pi/3 && a < 5*pi/3)
            angle = a - pi/3;
        elseif (a > 2*pi/3 && a < pi) || (a > 5*pi/3 && a < 2*pi)
            angle = a - 2*pi/3;
        else
            angle = a;
        end
        border = sqrt(3)*BS.cell_radius / (2*sin(pi/3 + angle));  % Hexagonal border
        
        while 1
            dist = rand()*abs(border);  % distance m
            if dist > 35
                break
            end
        end
        pathloss = 128.1 + 37.6 * log10(dist / 1000);  % pathloss dB
        Users(u).pathloss = 10^(-pathloss/10);  % dB2linear
        H_matrix_temp = zeros(BS.Nt, BS.num_of_SCs);
        a_vector_temp = zeros(BS.Nt, 1);  % Nr x Nt
        for i = 1:BS.Nt
            a_vector_temp(i) = exp(1i*sin(a - 2*pi*(i-1)/BS.Nt));
        end
        for loop_SC = 1 : BS.num_of_SCs
            h = sqrt(BS.Rician_fading_factor/(BS.Rician_fading_factor+1))*a_vector_temp.' + ...
                sqrt(1/(BS.Rician_fading_factor+1))*(sqrt(1/2)*(randn(1,BS.Nt)+1i*randn(1,BS.Nt)));
            h = sqrt(Users(u).pathloss)*h.';
            H_matrix_temp(:, loop_SC) = h;
        end % end of loop_SC
        Users(u).ang = a;
        Users(u).a_vector = a_vector_temp/norm(a_vector_temp);
        Users(u).H_matrix = H_matrix_temp;  % Chnanel matrix
        Users(u).dist = dist;
        Users(u).coor = dist*exp(1i*a);  % coordinate
        Users(u).if_SIC_capability = (rand <= BS.ue_SIC_ratio);
        Users(u).cost_weight = 0.2;
        Users(u).pow = BS.tx_power/BS.num_of_UE;
    end
    %% user coalition: matching (calculate the non-orthogonality)
    % initalization: the user Groups
    [Groups(1:BS.num_of_SCs).user_subset] = deal(0);
    [Groups(1:BS.num_of_SCs).MA_mode] = deal(0);
    [Groups(1:BS.num_of_SCs).SC_idx] = deal(0);
    % 0: OMA mode, 1: P-NOMA mode; 2: S-NOMA mode; 3: Hybrid NOMA mode
    User_Group_Formation_Process(); % output:UE Groups
    fprintf('loop_drop = %d, user coaltion formation. \n', loop_drop);
    %% Power Allocation
    % Input:  p_vector: the transmit power for the each UE
    p0_vector = ones(BS.num_of_UE, 1)*(.9*BS.tx_power/BS.num_of_UE);
    for loop_UE_idx = 1 : BS.num_of_UE
        if Users(loop_UE_idx).if_P_NOMA == 1 && Users(loop_UE_idx).if_Strong_UE_P_NOMA == 0
            % this UE is a weak signal user
            p0_vector(loop_UE_idx) = p0_vector(loop_UE_idx)*1.6;
        end
        if Users(loop_UE_idx).if_P_NOMA == 1 && Users(loop_UE_idx).if_Strong_UE_P_NOMA == 1
            % this UE is a strong signal user
            p0_vector(loop_UE_idx) = p0_vector(loop_UE_idx)*0.4;
        end
    end
    %% linear inequalities: maximum power constraint
    a_vector = ones(1, BS.num_of_UE);
    b = BS.tx_power;
    %% Linear equality Constraint
    Aeq = [];
    beq = [];
    %% Lower and Upper bounds of variables
    lb = 1e-3*ones(BS.num_of_UE, 1);
    ub = BS.tx_power*ones(BS.num_of_UE, 1);
    %% Nonlinear Constraint
    nonlcon = @my_non_con;
    %% Optimization Process
    options = optimoptions('fmincon','Display','final-detailed','Algorithm','sqp', 'MaxFunctionEvaluations', 20000);
    [p_vector_opt, fval] = fmincon(@obj_fun_mimo_noma,p0_vector, a_vector, b, Aeq, beq, lb, ub, nonlcon, options);
    %% Output
    [f_tot, u_fun, rate_ue, rate_ue_norm, SINR, SINR_SIC] = data_rate_cal(p_vector_opt);
    for u = 1 : BS.num_of_UE
        SC_idx = Users(u).SC_idx;
        [cost_UE] = g_cost_cal(u, SC_idx, SINR_SIC(u), 1);
        Users(u).pow = p_vector_opt(u);
        Users(u).SINR = SINR(u);
        Users(u).SINR_dB = 10*log10(SINR(u));
        Users(u).SINR_SIC = SINR_SIC(u);
        Users(u).u_fun = u_fun(u);
        Users(u).rate = rate_ue(u);
        if Users(u).rate/BS.sys_SC_bandwidth < 1e-3
            fprintf('Check.\n');
        end
        Users(u).rate_norm = rate_ue(u)/BS.rate_max;
    end % end of "for i = 1 : BS.num_of_UE"
    %% Data Saving
    fprintf('loop_drop = %d, Finsih all and Save data!\n', loop_drop);
    savefile = sprintf('../Data/MIMO-NOMA/Case_%d_File_idx_%d_Result.mat', CASE_IDX, loop_drop);
    save(savefile, 'BS', 'Users', 'Groups');
end % end of loop_drop