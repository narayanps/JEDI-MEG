function [A_reconst] = re_arrange_A(loc_true,loc_est, A_est, num_sources, P)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
[comb_id, ~, ~] = check_loc_error(loc_true, loc_est, num_sources);
A_reconst=zeros(size(A_est));
for p=1:P
    for i=1:num_sources
        for j=1:num_sources
            A_reconst(comb_id(i),comb_id(j),p) = A_est(i,j,p);
        end
        
    end
end

