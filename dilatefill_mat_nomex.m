function [smask ring iimg] = dilatefill_mat(iimg,smask,reversezeroedFT,zeroedFTvec,stdzeroedstack,stdF,NFFT,rindex,cindex,step)


step = int8(step);
[r c] = size(smask);
dilateimg = false(r,c);
ring = false(r,c);
startsmask = false(r,c);
% totarea = 1;
while ~isequal(smask,startsmask)% && totarea < 100
    startsmask = smask;
    rr1 = max(1,(min(rindex)-1));
    rr2 = min(r,(max(rindex)+1));
    cc1 = max(1,(min(cindex)-1));
    cc2 = min(c,(max(cindex)+1));
    subsmask = smask(rr1:rr2,cc1:cc2);
    dilateimg(rr1:rr2,cc1:cc2) = imdilate(subsmask,ones(3,3));
%     DD = morphmex('dilate_binary_ones33',logical(subsmask),ones(3),zeros(3),-1);
    dilateimg(rr1:rr2,cc1:cc2) = DD;
    
    [rindex cindex] = find(dilateimg);
    indsize = size(rindex,1);
%     totarea = totarea + indsize;
    for ind = 1:indsize
        if iimg(rindex(ind),cindex(ind)) == 127
            a = reversezeroedFT(:,rindex(ind),cindex(ind)).*zeroedFTvec;
            xcorrcoeff = ifft(a);
%             crosscorr = ifft(a);
%             xcorrcoeff = crosscorr/(stdzeroedstack(rindex(ind),cindex(ind))*stdF);         
            [~, b] = max(xcorrcoeff);
            if b > (NFFT)/2
                iimg(rindex(ind),cindex(ind)) = b - NFFT;
            else
                iimg(rindex(ind),cindex(ind)) = b;
            end    
        end
    end
    smask(rr1:rr2,cc1:cc2) = findstep_mex(dilateimg(rr1:rr2,cc1:cc2), iimg(rr1:rr2,cc1:cc2), step);
end


if step ~= 0
    subsmask = false(rr2-rr1+1,cc2-cc1+1);
    smaskCC = bwconncomp(smask(rr1:rr2,cc1:cc2));
    numPixels = cellfun(@numel,smaskCC.PixelIdxList);
    [~,idx] = max(numPixels);
    subsmask(smaskCC.PixelIdxList{idx}) = 1;
    smask = false(r,c);
    smask(rr1:rr2,cc1:cc2) = subsmask;
end




[rindex cindex] = find(smask);
rr1 = max(1,(min(rindex)-1));
rr2 = min(r,(max(rindex)+1));
cc1 = max(1,(min(cindex)-1));
cc2 = min(c,(max(cindex)+1));
subsmask = smask(rr1:rr2,cc1:cc2);
% ring(rr1:rr2,cc1:cc2) = ~subsmask & morphmex('dilate_binary_ones33',logical(subsmask),ones(3),zeros(3),-1);
ring(rr1:rr2,cc1:cc2) = ~subsmask & imdilate(subsmask,ones(3,3));