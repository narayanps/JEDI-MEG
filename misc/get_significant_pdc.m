function [PDC_sig, PDC_bin, f] = get_significant_pdc(x, Am, S, numsurr, nfft, fc, Ns, p)
S1=S .* eye(Ns);
for is=1:numsurr

    for ii=1:Ns
        for jj=1:Ns
            if ii~=jj
                xs=surrVCFTd(x,Am,S1,ii,jj);
                [Ams,Sus]=idMVAR(xs,p,2);
                Sus = Sus .* eye (Ns);
                [~,~,~,gpdcs,~,~,~,~,~,~,~] = fdMVAR(Ams,Sus,nfft,fc);
                PDCs(ii,jj,:,is)=abs(gpdcs(ii,jj,:)).^2;
            end
        end
    end
end

PDCsth=prctile(PDCs,95,4);

[~,~,~,gpdc,~,~,~,~,~,~,f] = fdMVAR(Am,S1,nfft,fc);
PDC=abs(gpdc).^2;
PDC(PDC < PDCsth) = 0;
PDC_sig=squeeze(mean(PDC(:,:,:),3));
pdc_surr = mean(PDCsth(:,:,:),3);
PDC_sig(PDC_sig < pdc_surr) = 0;
PDC_bin = PDC_sig - diag(diag(PDC_sig));
PDC_bin(PDC_bin>0) = 1;