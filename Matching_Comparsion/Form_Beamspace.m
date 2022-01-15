function [ ] = Form_Beamspace()
global BS; % type: structure
global Users;
global Beamspaces;
%% Calculation
K = BS.num_beam;
X = [Users(1:BS.num_of_UE).ang].'; % vector: num_of_BS X 1

[Beam_idx_vector, Theta_beam] = kmeans(X, K);

%% Output
for loop_UE_idx = 1 : BS.num_of_UE
    Users(loop_UE_idx).beam_idx = Beam_idx_vector(loop_UE_idx); 
end % end of loop_UE_idx
for loop_beam_idx = 1 : BS.num_beam
    % determine the user set of each beamspace
    Beamspaces(loop_beam_idx).users_set = find(Beam_idx_vector == loop_beam_idx);
    % determine the angle of each beamspace
    Beamspaces(loop_beam_idx).beam_theta =  Theta_beam(loop_beam_idx);
end % end of loop_beam_idx
end