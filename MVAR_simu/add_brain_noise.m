function [Y, sigma_b] = add_brain_noise(Y, G, snr, n_noise_sources)
%UNTITLED2 Summary of this function goes here
T=size(Y,2);
ind_noise_rand=randperm(size(G,2));
noise_ids=ind_noise_rand(1:n_noise_sources)';   
pn = mkpinknoise(T, n_noise_sources)';
Y_brain_noise = G(:, noise_ids)*pn;
bio_snr_=trace(Y*Y')./trace(Y_brain_noise*Y_brain_noise');
Y_brain_noise = Y_brain_noise*sqrt(bio_snr_/snr);
Y = Y + Y_brain_noise;
sigma_b = sqrt(bio_snr_/snr);
end

