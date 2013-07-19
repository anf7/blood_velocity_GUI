function imstack = imreg(imstack,interpfact,lags,ignoremask)

[z,r,c] = size(imstack);
imstack = double(imstack);
if nargin < 4
    ignoremask = zeros(r,c);
else
    ignoremask = double(ignoremask);
end

imstack_nans = imstack;

ignoremask(ignoremask == 1) = nan;
ignoremask(ignoremask == 0) = 1;
for n = 1:z
    imstack_nans(n,:,:) = squeeze(imstack_nans(n,:,:)).*ignoremask;
end
    

a = squeeze(imstack_nans(1,:,:));
    
if interpfact > 0
    a_i = interp2(a,interpfact);
else
    a_i = a;
end

[r_i c_i] = size(a_i);
    
rnlags = -lags;
rplags = lags;
cnlags = -lags;
cplags = lags;
rshift = 0;
cshift = 0;
rwindowsize = rplags - rnlags + 1;
cwindowsize = cplags - cnlags + 1;
    
for n = 1:z 
    clc
    disp(['Registering frame ' num2str(n) ' of ' num2str(z)]);
    
    b = squeeze(imstack_nans(n,:,:));
    if interpfact > 0
        b_i = interp2(b,interpfact);
        b_i_nonans = interp2(squeeze(imstack(n,:,:)),interpfact);
    else
        b_i = b;
        b_i_nonans = squeeze(imstack(n,:,:));
    end
    d_i = nan(size(b_i));
    
    repeat = true;
    while repeat
        repeat = false;
        cc = diff2lags(a_i,b_i,rnlags,rplags,cnlags,cplags);

        minv = min(cc(:));
        d = (cc == minv);
        [i,j] = find(d,1);

        rshift = ((i-(lags+1)) + rshift);
        cshift = ((j-(lags+1)) + cshift);

        if i == 1 || i == rwindowsize || j == 1 || j == cwindowsize
            rnlags = -lags + rshift;
            rplags = lags + rshift;
            cnlags = -lags + cshift;
            cplags = lags + cshift;
            repeat = true;
        end
    end

    d_i(max(1,1+rshift):min(r_i,r_i+rshift),max(1,1+cshift):min(c_i,c_i+cshift)) = ...
        b_i_nonans(max(1,1-rshift):min(r_i,r_i-rshift),max(1,1-cshift):min(c_i,c_i-cshift));

    if interpfact > 0
        for ir = 1:r-1
            for ic = 1:c-1
                imstack(n,ir,ic) = sum(sum(d_i((ir-1)*2^interpfact+1:ir*2^interpfact,(ic-1)*2^interpfact+1:ic*2^interpfact)))/(2^interpfact)^2;
            end
        end
    else
        imstack(n,:,:) = d_i;
    end
    
    if n ~= z
        rnlags = -lags + rshift;
        rplags = lags + rshift;
        cnlags = -lags + cshift;
        cplags = lags + cshift;
    end
end

[~, r, c] = size(imstack);
if interpfact ~= 0
    imstack = imstack(:,1:r-1,1:c-1);
end
[~, r, c] = size(imstack);   
nanmask = squeeze(isnan(sum(imstack,1)));
firstrow = r;
lastrow = 1;
firstcol = c;
lastcol = 1;
for n = 1:r
    if sum(nanmask(n,:)) ~= c
        if n < firstrow
            firstrow = n;
        end
        if n > lastrow
            lastrow = n;
        end
    end
end

for n = 1:c
    if sum(nanmask(:,n)) ~= r
        if n < firstcol
            firstcol = n;
        end
        if n > lastcol
            lastcol = n;
        end
    end
end

cropimg = imstack(:,firstrow:lastrow,firstcol:lastcol);
imstack = cropimg;


