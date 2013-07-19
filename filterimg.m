function filteredimg = filterimg(input,gblur)

h = fspecial('gaussian',[gblur gblur]);
[r c z] = size(input);

filteredimg = zeros(r,c,z);

for n = 1:z
    slice = input(:,:,n);
    done = false;
    holemask = (slice == 0);
    smallholes = holemask - bwareaopen(holemask,20);

    while ~done
        startimg = smallholes;
        [i j] = find(smallholes(2:r-1,2:c-1));
        slicep = zeros(r,c);
        for m = 1:length(i)
            boundpix = smallholes(i(m):i(m)+2,j(m):j(m)+2);
            denom = 9 - sum(boundpix(:));
            if denom ~= 0
                avgbound = sum(sum(slice(i(m):i(m)+2,j(m):j(m)+2)))/denom;
                slicep(i(m)+1,j(m)+1) = avgbound;
            end
        end
        slice = slice + slicep;
        holemask = (slice == 0);
        smallholes = holemask - bwareaopen(holemask,9);
        if isequal(smallholes,startimg);
            done = true;
        end
    end
    
    mask = (slice ~= 0);
    
    slicef = medfilt2(slice,[3 3]);   
    if gblur ~= 0
        slicef = imfilter(slice,h);
    end
    slicef = slicef.*mask;
    filteredimg(:,:,n) = slicef;
end
