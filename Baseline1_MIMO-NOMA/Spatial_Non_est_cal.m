function ND_S_est = Spatial_Non_est_cal(ue_n_idx, ue_i_idx)
%% Input
% mode = 0, only calculate based on angel info.
% mode = 1, calculate based on full CSI
%% Initializtion
global Users;
%% Processing
ND_S_est = norm((Users(ue_n_idx).a_vector)'*Users(ue_i_idx).a_vector);
end