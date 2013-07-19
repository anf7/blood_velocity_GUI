function [wheelimg rgbvel dirimg] = drawmaps(movement_matrix,minvel,maxvel,minspeed_threshold,maxspeed_threshold)

if nargin < 3
    minspeed_threshold = 0;
    maxspeed_threshold = inf;
end
if isempty(minvel)
    minvel = 0;
end

r = size(movement_matrix,1);
c = size(movement_matrix,2);

velocity = squeeze(movement_matrix(:,:,3));

udimg = -1*movement_matrix(:,:,1);
lrimg = movement_matrix(:,:,2);

wheelrad = 50;
[wheelimg rangle gangle bangle] = generate_colorwheel(wheelrad);

dirimg = zeros(r,c,3);
for n = 1:c
    for m = 1:r
        if velocity(m,n) > minspeed_threshold && velocity(m,n) < maxspeed_threshold
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
end

velocity(isinf(velocity)) = nan;
velocity(velocity > maxvel) = maxvel;
velocity(velocity < minvel) = minvel;
scalevelocity = velocity - minvel;
uint8vel = uint8(255*scalevelocity/max(scalevelocity(:)));


rgbvel = zeros(r,c,3);
jetmat = jet(256);

for n = 1:c
    for m = 1:r
        if ~isequal(dirimg(m,n,:),zeros(1,1,3)) && velocity(m,n) > minspeed_threshold && velocity(m,n) < maxspeed_threshold
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
