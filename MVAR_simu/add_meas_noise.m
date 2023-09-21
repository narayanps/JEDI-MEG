function [Y, sigma, est_snr] = add_meas_noise(Y, snr)
%UNTITLED2 Summary of this function goes here
%author : Narayan Subramaniyam
%Aalto/NBE
% resonable to use this if you have just one type of sensors

num_chan = size(Y,1);

sig_var = zeros(num_chan,1);

for i = 1:1:num_chan
    sig_var(i,1) = var(Y(i,:));
end

avg_var = mean(sig_var);
C=(avg_var / snr) * eye(num_chan);
noise = chol(C) * randn(size(Y));
est_snr = 0;
for i = 1:1:size(Y,1)
    est_snr = est_snr + (var(Y(i,:)) / var(noise(i,:)));
end

est_snr = est_snr/num_chan ;
sigma=chol(C(1,1));
Y = Y + noise;
end

