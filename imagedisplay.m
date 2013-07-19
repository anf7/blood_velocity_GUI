function imagedisplay(movement_matrix,scalefactor,togetherimg)

noisethresh = 1;


if nargin < 3
    togetherimg = ones(size(meanimg));
end

if nargin < 2
    scalefactor = 1;
end





r = size(movement_matrix,1);
c = size(movement_matrix,2);

udimg = -1*movement_matrix(:,:,1);
lrimg = movement_matrix(:,:,2);

noisemaskimg = logical(udimg) + logical(lrimg);
noisemaskimg = bwareaopen(noisemaskimg,noisethresh);
udimg = udimg.*noisemaskimg;
lrimg = lrimg.*noisemaskimg;


wheelrad = 50;
[wheelimg rangle gangle bangle] = generate_colorwheel(wheelrad);

dirimg = zeros(r,c,3);
for n = 1:c
    for m = 1:r
        angle = round(atand(udimg(m,n)/lrimg(m,n)));
        if ~isnan(angle)
            if lrimg(m,n) < 0
                angle = 180 + angle; 
            elseif udimg(m,n) < 0 && lrimg(m,n) >= 0
                angle = 360 + angle;
            end
            if angle == 0
                angle = 360;
            end
            dirimg(m,n,1) = rangle(1,angle);
            dirimg(m,n,2) = gangle(1,angle);
            dirimg(m,n,3) = bangle(1,angle);
        end
    end
end


velocity = noisemaskimg.*squeeze(movement_matrix(:,:,3));

maxvel = max(velocity(:));

% velocity = uint8(255*velocity/35);
% velrgb = ind2rgb(velocity,jet);
% for m = 1:r
%     for n = 1:c
%         if isequal(dirimg(m,n,:),zeros(1,1,3))
%             velrgb(m,n,:) = zeros(1,1,3);
%         end
%     end
% end

% figure(10),image(dirimg),axis off
% figure(11),imagesc(velocity),axis off,caxis([0 30])
% figure(1),imagesc(meanimg),axis off
% figure(1),image(wheelimg),axis off,axis image





velocity(velocity > maxvel) = maxvel;

uint8vel = uint8(255*velocity/maxvel);

rgbvel = zeros(r,c,3);
jetmat = jet(256);

for n = 1:c
    for m = 1:r
        if ~isequal(dirimg(m,n,:),zeros(1,1,3))
            indx = uint8vel(m,n)+1;
            rr = jetmat(indx,1);
            gg = jetmat(indx,2);
            bb = jetmat(indx,3);
            rgbvel(m,n,1) = rr;
            rgbvel(m,n,2) = gg;
            rgbvel(m,n,3) = bb;
        end
    end
end

% filteredvelimg = filterimg(velocity,5);
% filtereddirimg = filterimg(dirimg,5);

togetherimg(togetherimg == 0) = NaN;
togetherimg = togetherimg - min(togetherimg(:));
togetherimg(isnan(togetherimg)) = 0;
togetherimg = togetherimg/max(togetherimg(:));


dirimg(:,:,1) = dirimg(:,:,1).*togetherimg;
dirimg(:,:,2) = dirimg(:,:,2).*togetherimg;
dirimg(:,:,3) = dirimg(:,:,3).*togetherimg;

rgbvel(:,:,1) = rgbvel(:,:,1).*togetherimg;
rgbvel(:,:,2) = rgbvel(:,:,2).*togetherimg;
rgbvel(:,:,3) = rgbvel(:,:,3).*togetherimg;



maxmmpersec = maxvel*scalefactor;


figure(2),image(dirimg),title('Direction'),axis off,axis image

figure(3),image(rgbvel),title('Velocity Magnitude'), axis off,axis image,colormap(jet(67)),
colorbar('YTickLabel',...
    {'0',num2str((1/6)*maxmmpersec),num2str((2/6)*maxmmpersec),num2str((3/6)*maxmmpersec),...
     num2str((4/6)*maxmmpersec),num2str((5/6)*maxmmpersec),num2str(maxmmpersec)},'YTick',[1 12 23 34 45 56 67])





% maskimg = (squeeze(movement_matrix(:,:,3)) ~= 0);
% newdiam_img = findgeometery(maskimg);
% save images dirimg filtereddirimg filteredvelimg movement_matrix wheelimg
% save('gfpmcherry_post.mat')