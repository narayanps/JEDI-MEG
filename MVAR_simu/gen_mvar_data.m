function [data, A] = gen_mvar_data(A, T, T0, P, M)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
x=randn(M,T+T0);
y=x;
for i=P+1:T+T0
    yloc=reshape(fliplr(y(:,i-P:i-1)),[],1);
    y(:,i)=A*yloc+x(:,i);
end
data=y(:,T0+1:end);
end

