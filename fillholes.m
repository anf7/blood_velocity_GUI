function outputimg = fillholes(inputimg, hbtot)


ringlut = makelut('x(2,2) == 0 && sum(x(:)) ~= 0',3');
outputimg = zeros(size(inputimg));

invimg = ~(inputimg);

[label_img, num] = bwlabeln(invimg, 4);

for n = 1:num
    holeimg = (label_img == n);
    ringimg = applylut(holeimg, ringlut);
    sizehole = size(find(holeimg),1);
    sizering = size(find(ringimg),1);
      
    holehbtot = holeimg.*hbtot;
    ringhbtot = ringimg.*hbtot;
    ringavghb = sum(ringhbtot(:))/sizering;
    
    highhole = size(find(holehbtot > (0.9*ringavghb)),1);
    
    
    if highhole/sizehole > 0.9 || sizehole < 9
        outputimg = outputimg + holeimg;
    end
end

outputimg = inputimg + outputimg;