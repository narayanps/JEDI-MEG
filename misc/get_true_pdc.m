function PDC = get_true_pdc(A, nfft, fc, Ns)
[~,~,~,gpdc,~,~,~,~,~,~,~] = fdMVAR(A,eye(Ns),nfft,fc);
PDC=abs(gpdc).^2; % partial directed coherence
PDC=mean(PDC,3);
end