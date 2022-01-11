function [f] = obj_fun_mimo_noma(p_vector)
[f_tot, u_fun, rate_ue, rate_ue_norm, SINR, SINR_SIC] = data_rate_cal(p_vector);
f = -f_tot;
if ~isreal(f)
    fprintf('error in function "obj_fun_mimo_noma.m"!\n');
end
end % end of function