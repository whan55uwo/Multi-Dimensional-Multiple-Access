function [ ] = ZF_precoder_at_SC()
global BS; % type: structure
global Users;
global Coalitions;
%% Initialization
for Coalition_idx = 1 : BS.num_of_SCs % index of coalition
    UE_set = Coalitions(Coalition_idx).user_subset;
    UE_set_size = length(UE_set);
    if UE_set_size == 1
        continue;
    end
    if Coalitions(Coalition_idx).MA_mode == 1
        continue;
    end
    %% ZF PRECODING
    for SC_idx = 1 : BS.num_of_SCs % index of SC
        %% Coalitions(Coalition_idx).MA_mode == 2, Spatial-domain NOMA
        if Coalitions(Coalition_idx).MA_mode == 2
            H_mat = zeros(BS.Nt, UE_set_size);
            for j = 1 : UE_set_size
                UE_idx = UE_set(j);
                H_mat(:, j) = Users(UE_idx).H_matrix(:, SC_idx);
            end
            W_ZF = H_mat*pinv(H_mat'*H_mat);
            W_ZF=W_ZF./repmat(sqrt(sum(W_ZF.*conj(W_ZF),1)),size(W_ZF,1),1); % normalization
            for j = 1 : UE_set_size
                UE_idx = UE_set(j);
                Users(UE_idx).beam_vector_matrix(:,SC_idx) = W_ZF(:, j);
            end
        end
        %% Coalitions(Coalition_idx).MA_mode == 3, Hybrid NOMA
        if Coalitions(Coalition_idx).MA_mode == 3
            H_mat = zeros(BS.Nt, UE_set_size);
            for j = 1 : UE_set_size
                UE_idx = UE_set(j);
                H_mat(:, j) = Users(UE_idx).H_matrix(:, SC_idx);
            end
            W_ZF = H_mat*pinv(H_mat'*H_mat);
            W_ZF=W_ZF./repmat(sqrt(sum(W_ZF.*conj(W_ZF),1)),size(W_ZF,1),1); % normalization
            for j = 1 : UE_set_size
                UE_idx = UE_set(j);
                if isempty(Users(UE_idx).partner_P_NOMA)
                    Users(UE_idx).beam_vector_matrix(:,SC_idx) = W_ZF(:, j);
                end
            end
            for j = 1 : UE_set_size
                UE_idx = UE_set(j);
                if Users(UE_idx).if_Strong_UE_P_NOMA == 1
                    UE_weak_idx = Users(UE_idx).partner_P_NOMA;
                    angle_set_strong_ue = Users(UE_idx).ang;
                    angle_set_weak_ue = Users(UE_weak_idx).ang;
                    ang_ = (angle_set_strong_ue + angle_set_weak_ue)/2;
                        a_vector_temp = zeros(BS.Nt, 1);
                        for i = 1:BS.Nt
                            a_vector_temp(i) = exp(1i*sin(ang_ - 2*pi*(i-1)/BS.Nt));
                        end
                    Users(UE_idx).beam_vector_matrix(:,SC_idx) = a_vector_temp/norm(a_vector_temp);%W_ZF(:, j);
                    Users(UE_weak_idx).beam_vector_matrix(:,SC_idx) = a_vector_temp/norm(a_vector_temp);%W_ZF(:, j);
                end
            end
        end
    end % end of SC_idx
end % end of Coalition_idx
end