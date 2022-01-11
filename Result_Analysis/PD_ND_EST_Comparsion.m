%% Add simulation results:  PD_ND_est VS PD_ND
Case_IDX_proposed = 63;
loop_drop_num_proposed = 100;
ND_PD_temp =  [];
ND_PD_EST_temp =  [];
u_cost_temp =  [];
u_cost_EST_temp =  [];
for loop_drop =  1 :  loop_drop_num_proposed
    clear -global
    fprintf('loop_drop = %d, IMDMA Scheme!\n', loop_drop);
    openfile1 =sprintf('../Data/IMDMA/Case_%d_File_idx_%d_Result.mat', Case_IDX_proposed, loop_drop);
    load(openfile1,  'BS', 'Users');
    ND_PD_temp =  [[Users(1:BS.num_of_UE).PD_ND]'; ND_PD_temp];
    ND_PD_EST_temp =  [[Users(1:BS.num_of_UE).PD_ND_est]'; ND_PD_EST_temp];
    u_cost_temp =  [[Users(1:BS.num_of_UE).utilization_cost]'; u_cost_temp];
    u_cost_EST_temp =  [[Users(1:BS.num_of_UE).utilization_cost_est]'; u_cost_EST_temp];
end
ND_PD_EST_error_temp = ND_PD_temp - ND_PD_EST_temp;
ND_PD_temp = ND_PD_temp(ND_PD_temp > 0.02);
ND_PD_EST_temp = ND_PD_EST_temp(ND_PD_EST_temp > 0.02);
ND_PD_EST_error_temp = ND_PD_EST_error_temp(ND_PD_EST_error_temp ~= 0);
u_cost_error_temp = u_cost_temp - u_cost_EST_temp;
%% Figure
figure(2)
Seg_Num = 100;
[ALL_Seg_ND_PD_Array1, ALL_ND_PD_CDF1]=cdf_calculate(ND_PD_temp,Seg_Num);
[ALL_Seg_ND_PD_EST_Array2, ALL_ND_PD_EST_CDF2]=cdf_calculate(ND_PD_EST_temp,Seg_Num);
plot(ALL_Seg_ND_PD_Array1, ALL_ND_PD_CDF1, '-', ALL_Seg_ND_PD_EST_Array2, ALL_ND_PD_EST_CDF2,'--');
xlabel('value ','FontSize',12);
ylabel('CDF(cumulative distribution function)','FontSize',12);
grid on
grid minor
legend('ND PD','ND PD EST');
figure(3)
[ALL_Seg_ND_PD_Array1, ALL_ND_PD_CDF1]=cdf_calculate(ND_PD_EST_error_temp,Seg_Num);
plot(ALL_Seg_ND_PD_Array1, ALL_ND_PD_CDF1, '-');
xlabel('value ','FontSize',12);
ylabel('CDF(cumulative distribution function)','FontSize',12);
grid on
grid minor
legend('Estimation Error of ND PD');
figure(4)
[f_error_ND_PD,xi_error_ND_PD]=ksdensity(ND_PD_EST_error_temp);
plot(xi_error_ND_PD, f_error_ND_PD, '-');
xlabel('value ','FontSize',12);
ylabel('CDF(cumulative distribution function)','FontSize',12);
grid on
grid minor
legend('Estimation Error of ND PD');
%%
figure(6)
Seg_Num = 100;
[ALL_Seg_u_cost_Array1, ALL_u_cost_CDF1]=cdf_calculate(u_cost_temp,Seg_Num);
[ALL_Seg_u_cost_est_Array2, ALL_u_cost_est_CDF2]=cdf_calculate(u_cost_EST_temp,Seg_Num);
plot(ALL_Seg_u_cost_Array1, ALL_u_cost_CDF1, '-', ALL_Seg_u_cost_est_Array2, ALL_u_cost_est_CDF2,'--');
xlabel('value ','FontSize',12);
ylabel('CDF(cumulative distribution function)','FontSize',12);
grid on
grid minor
legend('utiliztion cost','utiliztion cost estimation');
figure(7)
[ALL_Seg_u_cost_error_Array1, ALL_u_cost_error_CDF1]=cdf_calculate(u_cost_error_temp,Seg_Num);
plot(ALL_Seg_u_cost_error_Array1, ALL_u_cost_error_CDF1, '-');
xlabel('value ','FontSize',12);
ylabel('CDF(cumulative distribution function)','FontSize',12);
grid on
grid minor
legend('Estimation Error of utilization cost');
