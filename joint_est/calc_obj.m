function [f] = calc_obj(R, S_R, T)
logdetR = 2*sum(log(diag(chol(R))));
f = T*logdetR + trace(R*inv(S_R));
end

