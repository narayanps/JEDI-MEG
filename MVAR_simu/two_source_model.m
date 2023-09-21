function [data, A] = two_source_model(T)
M = 2; %number of sources;
P=2; %order of the model
T0=1000; %length of ignored start 

%Generate stable AR matrix
lambdamax=10;
while lambdamax > 1 || lambdamax < 0.9
  A=[];
  for k=1:P
    aloc = zeros(M);
    aloc([1 4 ]) = -0.9; %diagonal elements
    aloc([2]) = randn(1, 1);
    A=[A,aloc];
  end
  E=eye(M*P);AA=[A;E(1:end-M,:)];lambda=eig(AA);lambdamax=max(abs(lambda));
end

%Generate MVAR data
[data, A] = gen_mvar_data(A, T, T0, P, M);









