clc;
clear -global; clearvars;
%% Initalization
Utility_sum = cell(3, 1);
Utility_fun = cell(3, 1);
Utility_avg = cell(3,1);
SINR_dB = cell(3, 1);
SINR_dB_avg=cell(3,1);
Rate = cell(3, 1);
Rate_avg=cell(3,1);
pow_tot= cell(3,1);
ND_total = cell(3, 1);
ND_PD =  cell(3, 1);
ND_SD = cell(3, 1);
ND_total_avg = cell(3, 1);
ND_PD_avg =  cell(3, 1);
ND_SD_avg = cell(3, 1);
%% Data Collection: IMDMA for jounral paper
Case_IDX_proposed = 63;
loop_drop_num_proposed = 100;
Utility_sum_temp = []; % the utility sum of the proposed scheme
Utility_fun_temp = [];
Utility_avg_temp = [];
SINR_temp = [];
SINR_avg_temp=[];
Rate_avg_temp=[];
Rate_temp = [];
ND_total_temp = [];
ND_PD_temp =  [];
ND_SD_temp = [];
ND_total_avg_temp = [];
ND_PD_avg_temp =  [];
ND_SD_avg_temp = [];
pow_tot_temp = [];
for loop_drop =  1 :  loop_drop_num_proposed
    clear -global
    fprintf('loop_drop = %d, IMDMA Scheme!\n', loop_drop);
    openfile1 =sprintf('../Data/IMDMA/Case_%d_File_idx_%d_Result.mat', Case_IDX_proposed, loop_drop);
    load(openfile1,  'BS', 'Users');
    Utility_avg_temp=[[sum([Users(1:BS.num_of_UE).u_fun])/BS.num_of_UE];Utility_avg_temp];
    Utility_sum_temp=[[sum([Users(1:BS.num_of_UE).u_fun])];Utility_sum_temp];
    Utility_fun_temp = [[Users(1:BS.num_of_UE).u_fun]'; Utility_fun_temp];
    SINR_temp = [[Users(1:BS.num_of_UE).SINR_dB]'; SINR_temp];
    Rate_temp = [[Users(1:BS.num_of_UE).rate]'; Rate_temp];
    SINR_avg_temp=[[sum([Users(1:BS.num_of_UE).SINR_dB])/BS.num_of_UE];SINR_avg_temp];
    Rate_avg_temp=[[sum([Users(1:BS.num_of_UE).rate])/BS.num_of_UE];Rate_avg_temp];
    ND_total_temp =  [[Users(1:BS.num_of_UE).utilization_cost]'; ND_total_temp];
    ND_total_avg_temp = [sum([Users(1:BS.num_of_UE).utilization_cost])/BS.num_of_UE; ND_total_avg_temp];
    ND_PD_temp =  [[Users(1:BS.num_of_UE).PD_ND]'; ND_PD_temp];
    ND_SD_temp =  [[Users(1:BS.num_of_UE).SD_ND]'; ND_SD_temp];
    ND_PD_avg_temp =  [sum([Users(1:BS.num_of_UE).PD_ND])/BS.num_of_UE; ND_PD_avg_temp];
    ND_SD_avg_temp =  [sum([Users(1:BS.num_of_UE).SD_ND])/BS.num_of_UE; ND_SD_avg_temp];
    pow_tot_ = sum([Users(1:BS.num_of_UE).pow]);
    pow_tot_temp = [pow_tot_; pow_tot_temp];
end % end for loop_drop
Utility_sum{1,1} = Utility_sum_temp;
Utility_fun{1,1} = Utility_fun_temp;
Utility_avg{1,1} = Utility_avg_temp;
SINR_dB{1,1} = SINR_temp;
Rate{1,1} = Rate_temp;
SINR_dB_avg{1,1}=SINR_avg_temp;
Rate_avg{1,1}=Rate_avg_temp;
ND_total{1,1} = ND_total_temp;
ND_total_avg{1,1} = ND_total_avg_temp;
ND_PD{1,1} = ND_PD_temp;
ND_SD{1,1} = ND_SD_temp;
ND_PD_avg{1,1} = ND_PD_avg_temp;
ND_SD_avg{1,1} = ND_SD_avg_temp;
pow_tot{1,1} = pow_tot_temp;
%% Baseline 1 -RSMA
Case_IDX_Baseline_1 = 3;
loop_drop_num_Baseline_1 = 100;
Utility_sum_temp = [];
Utility_fun_temp = [];
Utility_avg_temp = [];
SINR_temp = [];
SINR_avg_temp=[];
Rate_avg_temp=[];
Rate_temp = [];
ND_total_temp = [];
ND_PD_temp =  [];
ND_SD_temp = [];
pow_tot_temp = [];
ND_total_avg_temp = [];
ND_PD_avg_temp =  [];
ND_SD_avg_temp = [];
for loop_drop =  1 : loop_drop_num_Baseline_1
    clear -global
    fprintf('loop_drop = %d, Baseline 1 RSMA!\n', loop_drop);
    openfile1 =sprintf('../Data/RSMA/Case_%d_File_idx_%d_Result.mat', Case_IDX_Baseline_1, loop_drop);
    load(openfile1,  'BS', 'Users', 'Groups');
    Utility_avg_temp=[[sum([Users(1:BS.num_of_UE).u_fun])/BS.num_of_UE];Utility_avg_temp];
    Utility_sum_temp=[[sum([Users(1:BS.num_of_UE).u_fun])];Utility_sum_temp];
    Utility_fun_temp = [[Users(1:BS.num_of_UE).u_fun]'; Utility_fun_temp];
    SINR_temp = [[Users(1:BS.num_of_UE).SINR_dB]'; SINR_temp];
    Rate_temp = [[Users(1:BS.num_of_UE).rate_total]'; Rate_temp];
    SINR_avg_temp=[[sum([Users(1:BS.num_of_UE).SINR_dB])/BS.num_of_UE];SINR_avg_temp];
    Rate_avg_temp=[[sum([Users(1:BS.num_of_UE).rate_total])/BS.num_of_UE];Rate_avg_temp];
    ND_total_temp =  [[Users(1:BS.num_of_UE).utilization_cost]'; ND_total_temp];
    ND_total_avg_temp = [sum([Users(1:BS.num_of_UE).utilization_cost])/BS.num_of_UE; ND_total_avg_temp];
    ND_PD_temp =  [[Users(1:BS.num_of_UE).PD_ND]'; ND_PD_temp];
    ND_SD_temp =  [[Users(1:BS.num_of_UE).SD_ND]'; ND_SD_temp];
    ND_PD_avg_temp =  [sum([Users(1:BS.num_of_UE).PD_ND])/BS.num_of_UE; ND_PD_avg_temp];
    ND_SD_avg_temp =  [sum([Users(1:BS.num_of_UE).SD_ND])/BS.num_of_UE; ND_SD_avg_temp];
    pow_tot_ = sum([Users(1:BS.num_of_UE).pow_private]) + ...
        sum([Groups(1:BS.num_of_SCs).pow_common]);
    pow_tot_temp = [pow_tot_; pow_tot_temp];
end % end for loop_drop
Utility_sum{2,1} = Utility_sum_temp;
Utility_fun{2,1} = Utility_fun_temp;
Utility_avg{2,1} = Utility_avg_temp;
SINR_dB{2,1} = SINR_temp;
Rate{2,1} = Rate_temp;
SINR_dB_avg{2,1}=SINR_avg_temp;
Rate_avg{2,1}=Rate_avg_temp;
ND_total_avg{2,1} = ND_total_avg_temp;
ND_total{2,1} = ND_total_temp;
ND_PD{2,1} = ND_PD_temp;
ND_SD{2,1} = ND_SD_temp;
ND_PD_avg{2,1} = ND_PD_avg_temp;
ND_SD_avg{2,1} = ND_SD_avg_temp;
pow_tot{2,1} = pow_tot_temp;
%% Baseline 2 - MIMO-NOMA 
Case_IDX_Baseline_2 = 3;
loop_drop_num_Baseline_2 = 100;
Utility_sum_temp = [];
Utility_fun_temp = [];
Utility_avg_temp = [];
SINR_temp = [];
SINR_avg_temp=[];
Rate_avg_temp=[];
Rate_temp = [];
ND_total_temp = [];
ND_PD_temp =  [];
ND_SD_temp = [];
ND_total_avg_temp = [];
ND_PD_avg_temp =  [];
ND_SD_avg_temp = [];
for loop_drop =  1 : loop_drop_num_Baseline_2
    clear -global
    fprintf('loop_drop = %d, Baseline 2 MIMO-NOMA!\n', loop_drop);
    openfile1 =sprintf('../Data/MIMO-NOMA/Case_%d_File_idx_%d_Result.mat', Case_IDX_Baseline_2, loop_drop);
    load(openfile1,  'BS', 'Users');
    Utility_avg_temp=[[sum([Users(1:BS.num_of_UE).u_fun])/BS.num_of_UE];Utility_avg_temp];
    Utility_sum_temp=[[sum([Users(1:BS.num_of_UE).u_fun])];Utility_sum_temp];
    Utility_fun_temp = [[Users(1:BS.num_of_UE).u_fun]'; Utility_fun_temp];
    SINR_temp = [[Users(1:BS.num_of_UE).SINR_dB]'; SINR_temp];
    Rate_temp = [[Users(1:BS.num_of_UE).rate]'; Rate_temp];
    SINR_avg_temp=[[sum([Users(1:BS.num_of_UE).SINR_dB])/BS.num_of_UE];SINR_avg_temp];
    Rate_avg_temp=[[sum([Users(1:BS.num_of_UE).rate])/BS.num_of_UE];Rate_avg_temp];
    ND_total_temp =  [[Users(1:BS.num_of_UE).utilization_cost]'; ND_total_temp];
    ND_total_avg_temp = [sum([Users(1:BS.num_of_UE).utilization_cost])/BS.num_of_UE; ND_total_avg_temp];
    ND_PD_temp =  [[Users(1:BS.num_of_UE).PD_ND]'; ND_PD_temp];
    ND_SD_temp =  [[Users(1:BS.num_of_UE).SD_ND]'; ND_SD_temp];
    ND_PD_avg_temp =  [sum([Users(1:BS.num_of_UE).PD_ND])/BS.num_of_UE; ND_PD_avg_temp];
    ND_SD_avg_temp =  [sum([Users(1:BS.num_of_UE).SD_ND])/BS.num_of_UE; ND_SD_avg_temp];
end % end for loop_drop
Utility_sum{3,1} = Utility_sum_temp;
Utility_fun{3,1} = Utility_fun_temp;
Utility_avg{3,1} = Utility_avg_temp;
SINR_dB{3,1} = SINR_temp;
Rate{3,1} = Rate_temp;
SINR_dB_avg{3,1}=SINR_avg_temp;
Rate_avg{3,1}=Rate_avg_temp;
ND_total{3,1} = ND_total_temp;
ND_total_avg{3,1} = ND_total_avg_temp;
ND_PD{3,1} = ND_PD_temp;
ND_SD{3,1} = ND_SD_temp; 
ND_PD_avg{3,1} = ND_PD_avg_temp;
ND_SD_avg{3,1} = ND_SD_avg_temp;
%% Data Analysis and ploting
Seg_Num = 100;
figure(1)
[ALL_Seg_Utility_sum_Array1, ALL_Utility_sum_CDF1]=cdf_calculate(Utility_sum{1,1},Seg_Num);
[ALL_Seg_Utility_sum_Array2, ALL_Utility_sum_CDF2]=cdf_calculate(Utility_sum{2,1},Seg_Num);
[ALL_Seg_Utility_sum_Array3, ALL_Utility_sum_CDF3]=cdf_calculate(Utility_sum{3,1},Seg_Num);
plot(ALL_Seg_Utility_sum_Array1, ALL_Utility_sum_CDF1, '--or', ALL_Seg_Utility_sum_Array2, ALL_Utility_sum_CDF2,'--pg',......
    ALL_Seg_Utility_sum_Array3, ALL_Utility_sum_CDF3,'--hb','LineWidth', 1);
xlabel('Sum of Utility functions ','FontSize',12);
ylabel('CDF(cumulative distribution function)','FontSize',12);
grid on
grid minor
legend('proposed-IMDMA','RSMA', 'MIMO-NOMA');
title([num2str(BS.num_of_UE),' UEs']);

figure(2)
[ALL_Seg_Utility_Array1, ALL_Utility_CDF1]=cdf_calculate(Utility_fun{1,1},Seg_Num);
[ALL_Seg_Utility_Array2, ALL_Utility_CDF2]=cdf_calculate(Utility_fun{2,1},Seg_Num);
[ALL_Seg_Utility_Array3, ALL_Utility_CDF3]=cdf_calculate(Utility_fun{3,1},Seg_Num);
plot(ALL_Seg_Utility_Array1, ALL_Utility_CDF1, '--or', ALL_Seg_Utility_Array2, ALL_Utility_CDF2,'--pg',......
    ALL_Seg_Utility_Array3, ALL_Utility_CDF3,'--hb','LineWidth', 1);
xlabel('Utility functions ','FontSize',12);
ylabel('CDF(cumulative distribution function)','FontSize',12);
grid on
grid minor
legend('proposed-IMDMA','RSMA', 'MIMO-NOMA');
title([num2str(BS.num_of_UE),' UEs']);

figure(3)
[ALL_Seg_Utility_avg_Array1, ALL_Utility_avg_CDF1]=cdf_calculate(Utility_avg{1,1},Seg_Num);
[ALL_Seg_Utility_avg_Array2, ALL_Utility_avg_CDF2]=cdf_calculate(Utility_avg{2,1},Seg_Num);
[ALL_Seg_Utility_avg_Array3, ALL_Utility_avg_CDF3]=cdf_calculate(Utility_avg{3,1},Seg_Num);
plot(ALL_Seg_Utility_avg_Array1, ALL_Utility_avg_CDF1, '--or', ALL_Seg_Utility_avg_Array2, ALL_Utility_avg_CDF2,'--pg',......
    ALL_Seg_Utility_avg_Array3, ALL_Utility_avg_CDF3,'--hb','LineWidth', 1);
xlabel('Averaged Utility function value ','FontSize',12);
ylabel('CDF(cumulative distribution function)','FontSize',12);
grid on
grid minor
legend('proposed-IMDMA','RSMA', 'MIMO-NOMA');
title([num2str(BS.num_of_UE),' UEs']);

figure(4)
[ALL_Seg_Rate_Array1, ALL_Rate_CDF1]=cdf_calculate(Rate{1,1}/BS.sys_SC_bandwidth,Seg_Num);
[ALL_Seg_Rate_Array2, ALL_Rate_CDF2]=cdf_calculate(Rate{2,1}/BS.sys_SC_bandwidth,Seg_Num);
[ALL_Seg_Rate_Array3, ALL_Rate_CDF3]=cdf_calculate(Rate{3,1}/BS.sys_SC_bandwidth,Seg_Num);
plot(ALL_Seg_Rate_Array1, ALL_Rate_CDF1, '--or', ALL_Seg_Rate_Array2, ALL_Rate_CDF2,'--pg',...
    ALL_Seg_Rate_Array3, ALL_Rate_CDF3,'--hb','LineWidth',1);
xlabel('Rate (bps/Hz)','FontSize',12);
ylabel('CDF(cumulative distribution function)','FontSize',12);
grid on
grid minor
legend('proposed-IMDMA','RSMA', 'MIMO-NOMA');
title([num2str(BS.num_of_UE),' UEs']);

figure(6)
[ALL_Seg_ND_Array1, ALL_ND_CDF1]=cdf_calculate(ND_total_avg{1,1},Seg_Num);
[ALL_Seg_ND_Array2, ALL_ND_CDF2]=cdf_calculate(ND_total_avg{2,1},Seg_Num);
[ALL_Seg_ND_Array3, ALL_ND_CDF3]=cdf_calculate(ND_total_avg{3,1},Seg_Num);
plot(ALL_Seg_ND_Array1, ALL_ND_CDF1, '--or', ALL_Seg_ND_Array2, ALL_ND_CDF2,'--pg',...
    ALL_Seg_ND_Array3, ALL_ND_CDF3,'--hb','LineWidth',1);
xlabel('Resource utilization cost','FontSize',12);
ylabel('CDF(cumulative distribution function)','FontSize',12);
grid on
grid minor
legend('proposed-IMDMA','RSMA', 'MIMO-NOMA');
title([num2str(BS.num_of_UE),' UEs']);
