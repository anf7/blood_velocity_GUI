function outimg = diff2lags(inimg1,inimg2,rnlags,rplags,cnlags,cplags)

[r1 c1] = size(inimg1);
[r2 c2] = size(inimg2);

outimg = zeros(rplags - rnlags + 1,cplags - cnlags + 1);

if r1 ~= r2 || c1 ~= c2
    error('Images must be the same size')
else
    for n = cnlags:cplags
        for m = rnlags:rplags
            a = double(inimg1(max(1,1+m):min(r1,r1+m),max(1,1+n):min(c1,c1+n)));
            b = double(inimg2(max(1,1-m):min(r1,r1-m),max(1,1-n):min(c1,c1-n)));
            diff = (a - b).^2;
            outimg(m-rnlags+1,n-cnlags+1) = nansum(diff(:))/((size(a,1)*size(a,2))-sum(isnan(diff(:))));
        end
    end
end