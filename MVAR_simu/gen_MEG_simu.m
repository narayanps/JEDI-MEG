function [Y, sigma_b, sigma_m, est_snr] = gen_MEG_simu(q,G,source_ids, meas_snr, bio_snr, num_bio_ns)
%GEN_MEG_SIMU Summary of this function goes here
%   Detailed explanation goes here


Y_brain_sig = G(:, source_ids) * q;

% add brain noise
[Y, sigma_b] = add_brain_noise(Y_brain_sig, G, bio_snr, num_bio_ns);

%add measurement noise
[Y, sigma_m, est_snr] = add_meas_noise(Y, meas_snr);


end

