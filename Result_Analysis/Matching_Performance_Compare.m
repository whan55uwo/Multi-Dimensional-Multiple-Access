%% Add simulation results:  vs 
Case_IDX_proposed = 63;
loop_drop_num_proposed = 100;
u_cost_EST_temp =  [];
for loop_drop =  1 :  loop_drop_num_proposed
    clear -global
    fprintf('loop_drop = %d, IMDMA Scheme!\n', loop_drop);
    openfile1 =sprintf('../Data/IMDMA/Case_%d_File_idx_%d_Result.mat', Case_IDX_proposed, loop_drop);
    load(openfile1,  'BS', 'Users');
    u_cost_EST_temp =  [[Users(1:BS.num_of_UE).utilization_cost_est]'; u_cost_EST_temp];
end
%% TWC by Yanan
Case_IDX_TWC = 3;
loop_drop_num_TWC = 200;
u_cost_EST_temp_TWC =  [];
for loop_drop =  1 :  loop_drop_num_TWC
    clear -global
    fprintf('loop_drop = %d, IMDMA Scheme!\n', loop_drop);
    openfile1 =sprintf('../Data/IMDMA/Case_for_compare_matching_%d_File_idx_%d_Result.mat', Case_IDX_TWC, loop_drop);
    load(openfile1,  'BS', 'Users');
    u_cost_EST_temp_TWC =  [[Users(1:BS.num_of_UE).utilization_cost_est]'; u_cost_EST_temp_TWC];
end
%% Figure
figure(1)
Seg_Num = 100;
[ALL_Seg_u_cost_Array1, ALL_u_cost_CDF1]=cdf_calculate(u_cost_EST_temp,Seg_Num);
[ALL_Seg_u_cost_est_Array2, ALL_u_cost_est_CDF2]=cdf_calculate(u_cost_EST_temp_TWC,Seg_Num);
plot(ALL_Seg_u_cost_Array1, ALL_u_cost_CDF1, '-', ALL_Seg_u_cost_est_Array2, ALL_u_cost_est_CDF2,'--');
xlabel('intensity of radio resource sharing conflicts','FontSize',12);
ylabel('CDF(cumulative distribution function)','FontSize',12);
grid on
grid minor
legend('Proposed Scheme','TWC');
