function [data, A] = three_source_model(T, type)
%if type = 0, one non-interacting source else source 1 -> source 2 ->
%source 3
M = 3; %number of sources;
P=2; %order of the model
T0=1000; %length of ignored start 

%Generate stable AR matrix
lambdamax=10;
while lambdamax > 1 || lambdamax < 0.9
  A=[];
  for k=1:P
    aloc = zeros(M);
    aloc([1 5 9]) = -0.9; %diagonal elements
    if type==0
        aloc([2]) = randn(1, 1);
    else
        aloc([2 6]) = randn(2, 1);
    end
    A=[A,aloc];
  end
  E=eye(M*P);AA=[A;E(1:end-M,:)];lambda=eig(AA);lambdamax=max(abs(lambda));
end


%Generate MVAR data
[data, A] = gen_mvar_data(A, T, T0, P, M);









