function [smask, ring, iimg] = dilatefill_mat_zerotest(iimg,smask,reversezeroedFT,zeroedFTvec,NFFT,rindex,cindex,step,ignore_areas_larger_than)


step = int8(step);
[r, c] = size(smask);
dilateimg = false(r,c);
ring = false(r,c);
startsmask = false(r,c);
while ~isequal(smask,startsmask)
    startsmask = smask;
    rr1 = max(1,(min(rindex)-1));
    rr2 = min(r,(max(rindex)+1));
    cc1 = max(1,(min(cindex)-1));
    cc2 = min(c,(max(cindex)+1));
    subsmask = smask(rr1:rr2,cc1:cc2);
    try
        DD = imdilate(subsmask,ones(3,3));
    catch
        DD = morphmex('dilate_binary_ones33',logical(subsmask),ones(3),zeros(3),-1);
    end
    dilateimg(rr1:rr2,cc1:cc2) = DD;
    [rindex, cindex] = find(dilateimg);
    indsize = size(rindex,1);
    if indsize >= ignore_areas_larger_than
        smask = [];
        break
    end
    for ind = 1:indsize
        if iimg(rindex(ind),cindex(ind)) == 127
            a = reversezeroedFT(:,rindex(ind),cindex(ind)).*zeroedFTvec;
            xcorrcoeff = ifft(a);        
            [~, b] = max(xcorrcoeff);
            if b > (NFFT)/2
                iimg(rindex(ind),cindex(ind)) = b - NFFT;
            else
                iimg(rindex(ind),cindex(ind)) = b;
            end    
        end
    end
    try
        smask(rr1:rr2,cc1:cc2) = findstep_mex(dilateimg(rr1:rr2,cc1:cc2), iimg(rr1:rr2,cc1:cc2), step);
    catch
        smask(rr1:rr2,cc1:cc2) = dilateimg(rr1:rr2,cc1:cc2) & (iimg(rr1:rr2,cc1:cc2) == step(rr1:rr2,cc1:cc2));
    end
end

[rindex, cindex] = find(smask);
rr1 = max(1,(min(rindex)-1));
rr2 = min(r,(max(rindex)+1));
cc1 = max(1,(min(cindex)-1));
cc2 = min(c,(max(cindex)+1));
subsmask = smask(rr1:rr2,cc1:cc2);
try
    ring(rr1:rr2,cc1:cc2) = ~subsmask & morphmex('dilate_binary_ones33',logical(subsmask),ones(3),zeros(3),-1);
catch
    ring(rr1:rr2,cc1:cc2) = ~subsmask & imdilate(subsmask,ones(3,3));
end