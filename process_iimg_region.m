function iimg = process_iimg_region(iimg,cc1,cc2,rr1,rr2,reversezeroedFT,zeroedFTvec,NFFT)


for nn = cc1:cc2
    for mm = rr1:rr2
        if iimg(mm,nn) == 127
            a = reversezeroedFT(:,mm,nn).*zeroedFTvec;
            crosscorr = ifft(a);        
            [~, b] = max(crosscorr);
            if b > (NFFT)/2
                iimg(mm,nn) = b - NFFT;
            else
                iimg(mm,nn) = b;
            end    
        end
    end
end