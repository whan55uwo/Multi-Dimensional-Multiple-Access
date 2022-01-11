function [ ] = Form_Beamspace()
global BS; % type: structure
global Users;
global Beamspaces;
%% Inialization
BS.num_beam = BS.num_of_UE; % Number of beamspaces in the cell;
[Beamspaces(1:BS.num_beam).users_set] = deal(0); % the user set assoicated with beamspace
%% Output
for loop_UE_idx = 1 : BS.num_of_UE
    Users(loop_UE_idx).beam_idx = loop_UE_idx;
end % end of loop_UE_idx
for loop_beam_idx = 1 : BS.num_beam
    % determine the user set of each beamspace
    Beamspaces(loop_beam_idx).users_set = loop_beam_idx;
end % end of loop_beam_idx
end