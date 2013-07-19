function rgbimg = setbackgroundtoblack(img,maxval,cmap)


cropimg = img;
cropimg(cropimg > maxval) = maxval;
cropimg = 255*cropimg/max(cropimg(:));
rgbimg = ind2rgb(uint8(cropimg),cmap);

for m = 1:size(rgbimg,1)
    for n = 1:size(rgbimg,2)
        if img(m,n) == 0
            rgbimg(m,n,:) = zeros(1,1,3); 
        end
    end
end


