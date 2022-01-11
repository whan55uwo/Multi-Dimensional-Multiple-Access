function [flag] = is_forbidden_coalition(UE_coalition_new, UE_idx)
%% INPUT:
% 1. UE coalition new Set, it contains the idx of each UE
% 2. UE index, the new UE add inthe coalition
%% Import the user data
global Users;
global Beamspaces;
global BS;
%% OUTPUT
flag = 0;
if length(UE_coalition_new) > BS.L_max
    flag = 1;
    return
end
% 0) OMA set
if length(UE_coalition_new) == 1 
    return
end
beam_space_idx = Users(UE_idx).beam_idx;
beam_user_set = Beamspaces(beam_space_idx).users_set;
set_ = intersect(UE_coalition_new, beam_user_set);
%1) only max two users in the same beamsapce can use PD-NOMA
if length(set_) > BS.N_max
    flag = 1;
    return;
end
% 2) hardware constraints
if length(set_) == BS.N_max && Users(UE_idx).if_SIC_capability == 0
    % if this UE does not have SIC capability
    distance = Users(UE_idx).dist;
    pair_UE_idx = setdiff(set_, [UE_idx]);
    distance_ = Users(pair_UE_idx).dist;
    if distance <= distance_
        flag = 1;
    end
end
end