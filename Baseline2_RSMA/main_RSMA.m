clc;
clear -global; clearvars;
close all; dbstop if error; digits(64);  % 64-bit operation precision, the default 
global BS; % type: structure
global Users;
global Beamspaces;
global Groups;
%% 1.Setting parameters
CASE_IDX = 5;
BS.cell_radius = 250;        % m
BS.sys_SC_bandwidth = 2*10^6;  % Hz
BS.tx_power_dBm = 33;         % dBm, TBD
BS.noise_density_dB = -140;       % dBm/Hz, TBD
BS.Nr = 1; % number of recieve antennas at UE
BS.Nt = 64; % number of transmit antennas
BS.num_of_UE   = 35; % the number of users in the cell, To Be Detemined (TBD)
BS.num_of_SCs = 12;  % the number of SCs
BS.L_max = 3; % the maximum number of UEs can share one SC
BS.Rician_fading_factor = 9; %  uncorrelated Rician fading
BS.sinr_max = 10^(20/10); % not in dB
BS.rate_max = 0.5*BS.sys_SC_bandwidth*log2(1+BS.sinr_max); % in bits
BS.sinr_th = 10^(6/10); % not in dB
BS.rate_min = 0.5*BS.sys_SC_bandwidth*log2(1+BS.sinr_th); % in bits
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
for loop_drop = 67:drop_num
    fprintf('%**************   loop_drop = %d, Begin!  ******************* %\n', loop_drop);
    BS.transmit_pow = deal(0);
    
    [Users(1:BS.num_of_UE).ang]        = deal(0);
    [Users(1:BS.num_of_UE).a_vector] = deal(0);
    [Users(1:BS.num_of_UE).H_matrix]          = deal(0);
    [Users(1:BS.num_of_UE).channel_gain]          = deal(0);
    [Users(1:BS.num_of_UE).dist]       = deal(0);
    [Users(1:BS.num_of_UE).coor]       = deal(0);
    [Users(1:BS.num_of_UE).pathloss]   = deal(0);
    [Users(1:BS.num_of_UE).candidates] = deal(0);
    [Users(1:BS.num_of_UE).if_SIC_capability] = deal(0); % hardware constraints: 1 MEANS HAVE SIC capability; 0, otherwises
    [Users(1:BS.num_of_UE).partner]    = deal(0); % the set of users share one SC with this UE
    [Users(1:BS.num_of_UE).beam_idx]    = deal(0); % which beamspace this UE belongs to
    [Users(1:BS.num_of_UE).PD_ND]    = deal(0); % the power_domain (PD) ND of UE, see fun "g_cost_cal"
    [Users(1:BS.num_of_UE).SD_ND]    = deal(0); % the spatial_domain (SD) ND of UE, see fun "g_cost_cal"
    [Users(1:BS.num_of_UE).ND]    = deal(0); % the spatial_domain (SD) + power_domain (PD) ND of UE, see fun "g_cost_cal"
    [Users(1:BS.num_of_UE).utilization_cost]    = deal(0); % the see fun "g_cost_cal"
    [Users(1:BS.num_of_UE).SC_idx]    = deal(0); % the assigned SC to this UE, see fun "g_cost_cal"
    [Users(1:BS.num_of_UE).SINR_SIC]    = deal(0); % The SINR at SIC receiver
    [Users(1:BS.num_of_UE).SINR_SIC_dB]    = deal(0); % The SINR at SIC receiver
    [Users(1:BS.num_of_UE).SINR]    = deal(0); % the data rate of UE, see fun "g_cost_cal"
    [Users(1:BS.num_of_UE).SINR_dB]    = deal(0); % the data rate of UE, see fun "obj_fun_P10"
    [Users(1:BS.num_of_UE).rate_total]    = deal(0); % the data rate of PRIVATE message, see fun "obj_fun_P10"
    [Users(1:BS.num_of_UE).rate_total_norm]    = deal(0); % the data rate of PRIVATE message, see fun "obj_fun_P10"
    [Users(1:BS.num_of_UE).rate_private]    = deal(0); % the data rate of PRIVATE message, see fun "obj_fun_P10"
    [Users(1:BS.num_of_UE).rate_common]    = deal(0); % the data rate of 	COMMON message, see fun "obj_fun_P10"
    [Users(1:BS.num_of_UE).u_fun]    = deal(0); % the data rate of UE, see fun "obj_fun_P10"
    [Users(1:BS.num_of_UE).pow_private]    = deal(0); % the data rate of UE, see fun "obj_fun_P10"
    [Users(1:BS.num_of_UE).beam_vector_private] = deal(0); % private message for private message
    [Users(1:BS.num_of_UE).cost_weight] = deal(0); % objective funtion
    
    sum_rate = deal(0);
    
    %% 3. randomly generate UEs and channel modeling
    for u = 1:BS.num_of_UE
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
        for loop_SC = 1 : BS.num_of_SCs
            h = deal(0);  % Nr x Nt
            for i = 1:BS.Nt
                h(i) = exp(1i*sin(a - 2*pi*(i-1)/BS.Nt));
            end
            Users(u).a_vector = h/norm(h);
            h = sqrt(BS.Rician_fading_factor/(BS.Rician_fading_factor+1))*h + ...
                sqrt(1/(BS.Rician_fading_factor+1))*(sqrt(1/2)*(randn(1,BS.Nt)+1i*randn(1,BS.Nt)));
            h = sqrt(Users(u).pathloss)*h.';
            H_matrix_temp(:, loop_SC) = h;
        end % end of loop_SC
        Users(u).ang = a;
        Users(u).H_matrix = H_matrix_temp;  % Chnanel matrix
        Users(u).dist = dist;
        Users(u).coor = dist*exp(1i*a);  % coordinate
        Users(u).if_SIC_capability = (rand <= BS.ue_SIC_ratio);
        Users(u).cost_weight = 0.2;
    end
    %% 4.Clustering and Beamforming: [clustering by \Delta \theta ]and beamforming design
    Form_Beamspace()
    fprintf('loop_drop = %d, forming beamspace. \n', loop_drop);
    %% user grouping: matching (calculate the non-orthogonality) and SC allocation
    % initalization: the user groups
    [Groups(1:BS.num_of_SCs).user_group] = deal(0);
    SC_allocation_matrix = User_Grouping_SC_Allocation(); % output: a 0-1 matrix show the allocation resulits
    fprintf('loop_drop = %d, user grouping. \n', loop_drop);
    %% 6. power allocation
    [Groups(1:BS.num_of_SCs).beam_vector_common] = deal(0); % preocder for the common message
    [Groups(1:BS.num_of_SCs).pow_common] = deal(0); % transmit power for the common message
    [Groups(1:BS.num_of_SCs).rate_common_total] = deal(0); % transmit power for the common message
    [Groups(1:BS.num_of_SCs).rate_common_total_norm] = deal(0); % transmit power for the common message
    [Groups(1:BS.num_of_SCs).a_vector] = deal(0); % the precentage of common message allocated to the user in group
    [p_c_vector_opt, p_p_vector_opt, a_vector_opt] = precoding_power_allocation();
    %% 9. Data Saving
    fprintf('loop_drop = %d, Finsih all!\n', loop_drop);
    savefile = sprintf('../Data/RSMA/Case_%d_File_idx_%d_Result.mat', CASE_IDX, loop_drop);
    save(savefile, 'BS', 'Users', 'Beamspaces', 'Groups');
end % end of loop_drop