clc;
clear -global; clearvars;
close all; dbstop if error; digits(64);  % 64-bit operation precision, the default
% IMDMA
global BS; % type: structure
global Users;
global Beamspaces;
global Coalitions;

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
BS.num_beam = 10; % Number of beamspaces in the cell, (TBD)
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
BS.tx_power = 10^(BS.tx_power_dBm / 10);  % in mW
BS.noise_power = BS.sys_SC_bandwidth*10^(BS.noise_density_dB/10);  % in mW
SNR_dB=10*log10(BS.tx_power/BS.noise_power); % TBD,
drop_num = 200; % run this simulation for "drop_num" times

% Computatinal Complexity of Algorithm 1 and 2
Num_iteration_alg1_record = zeros(drop_num, 1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('%% **************    Simulation Start!!! 1~%d drops.   ******************* %%\n', drop_num);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for loop_drop = 1:drop_num
    fprintf('%**************   loop_drop = %d, Begin!  ******************* %\n', loop_drop);
    BS.transmit_pow = deal(0);
    
    [Users(1:BS.num_of_UE).ang] = deal(0);
    [Users(1:BS.num_of_UE).a_vector] = deal(0);
    [Users(1:BS.num_of_UE).H_matrix] = deal(0);
    [Users(1:BS.num_of_UE).beam_vector_matrix] = deal(0);
    [Users(1:BS.num_of_UE).channel_gain] = deal(0);
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
    [Users(1:BS.num_of_UE).PD_ND_est]    = deal(0); % 
    [Users(1:BS.num_of_UE).SD_ND]    = deal(0); % the spatial_domain (SD) ND of UE, see fun "g_cost_cal"
    [Users(1:BS.num_of_UE).SD_ND_est]    = deal(0);
    [Users(1:BS.num_of_UE).utilization_cost]    = deal(0); % the see fun "g_cost_cal"
    [Users(1:BS.num_of_UE).utilization_cost_est]    = deal(0); % the see fun "g_cost_cal"
    [Users(1:BS.num_of_UE).SC_idx]    = deal(0); % the assigned SC to this UE, see fun "g_cost_cal"
    [Users(1:BS.num_of_UE).SINR_SIC]    = deal(0); 
    [Users(1:BS.num_of_UE).SINR]    = deal(0); % the data rate of UE, see fun "g_cost_cal"
    [Users(1:BS.num_of_UE).SINR_dB]    = deal(0); % the data rate of UE, see fun "obj_fun_P10"
    [Users(1:BS.num_of_UE).rate]    = deal(0); % the data rate of UE, see fun "obj_fun_P10"
    [Users(1:BS.num_of_UE).u_fun]    = deal(0); % the data rate of UE, see fun "obj_fun_P10"
    [Users(1:BS.num_of_UE).pow]    = deal(0); % the data rate of UE, see fun "obj_fun_P10"
    
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
        beam_vector_matrix_temp = zeros(BS.Nt, BS.num_of_SCs);
        a_vector_temp = zeros(BS.Nt, 1);  % Nr x Nt
        for i = 1:BS.Nt
            a_vector_temp(i) = exp(1i*sin(a - 2*pi*(i-1)/BS.Nt));
        end
        for loop_SC = 1 : BS.num_of_SCs
            h = sqrt(BS.Rician_fading_factor/(BS.Rician_fading_factor+1))*a_vector_temp.' + ...
                sqrt(1/(BS.Rician_fading_factor+1))*(sqrt(1/2)*(randn(1,BS.Nt)+1i*randn(1,BS.Nt)));
            h = sqrt(Users(u).pathloss)*h.';
            H_matrix_temp(:, loop_SC) = h;
            beam_vector_matrix_temp(:, loop_SC) = h/norm(h);
        end % end of loop_SC
        Users(u).ang = a;
        Users(u).a_vector = a_vector_temp/norm(a_vector_temp);
        Users(u).H_matrix = H_matrix_temp;  % Chnanel matrix
        Users(u).beam_vector_matrix = beam_vector_matrix_temp;
        Users(u).dist = dist;
        Users(u).coor = dist*exp(1i*a);  % coordinate
        Users(u).if_SIC_capability = (rand <= BS.ue_SIC_ratio);
        Users(u).cost_weight = 0.2;
    end
    
    %% 4.Clustering and Beamforming: [clustering by \Delta \theta ]and beamforming design
    Form_Beamspace()
    % beamforming vector
    fprintf('loop_drop = %d, forming beamspace. \n', loop_drop);
    %% user coalition: matching (calculate the non-orthogonality)
    % initalization: the user Coalitions
    [Coalitions(1:BS.num_of_SCs).user_subset] = deal(0);
    [Coalitions(1:BS.num_of_SCs).MA_mode] = deal(0);
    [Coalitions(1:BS.num_of_SCs).SC_idx] = deal(0);
    % 0: OMA mode, 1: P-NOMA mode; 2: S-NOMA mode; 3: Hybrid NOMA mode
    [sum_cost_est, utilization_cost_est_UE_array, iter_num_ag1]  = Many_to_one_matching_Yannan_TWC(); % output:UE coalitions
    Num_iteration_alg1_record(loop_drop) = iter_num_ag1;
    fprintf('loop_drop = %d, user coaltion formation. \n', loop_drop);
    %% 9. Data Saving
    fprintf('loop_drop = %d, Finsih all and Save data!\n', loop_drop);
    savefile = sprintf('../Data/IMDMA/Case_for_compare_matching_%d_File_idx_%d_Result.mat', CASE_IDX, loop_drop);
    save(savefile, 'BS', 'Users', 'Beamspaces', 'Coalitions','Num_iteration_alg1_record');
end % end of loop_drop