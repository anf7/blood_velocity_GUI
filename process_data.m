function [movement_matrix, togetherimg] = process_data(instruct)

% instruct.downsample
% instruct.registervid 
% instruct.dynamic_processing
% instruct.subpixel 
% instruct.createvid 
% instruct.framerate 
% instruct.pixpermm 
% instruct.fshift 
% instruct.seqlength 
% instruct.seqspacing 
% instruct.filename
% instruct.first 
% instruct.last 
% instruct.premask
mask_exclusion = false;

% (pix/frame) * Scale Factor = (mm/sec)
% Scale Factor = (frame/sec) * (mm/pix)
scale_factor = instruct.framerate*(instruct.pixpermm)^-1;

if ispc % Use correct syntax for current OS
    divider = '\';
else
    divider = '/';
end
divind = strfind(instruct.filename,divider);
divind = divind(end); % Find last directory divider symbol
filepath = instruct.filename(1:divind); % File path is all the text up to last divider symbol

imginfo = imfinfo(instruct.filename); % Get size info on TIFF file
r = imginfo(1).Height;
c = imginfo(1).Width;
fulllength = length(imginfo);

%% Load specified frames from multi-page TIFF file into the variable 'imstack'
if isequal(instruct.last,'end')
    instruct.last = fulllength;
elseif instruct.last > fulllength
    instruct.last = fulllength;
end

imstack = zeros((instruct.last - instruct.first + 1),r,c);
for n = instruct.first:instruct.last
    imstack(n-instruct.first+1,:,:) = double(imread(instruct.filename,n));
end


%%  If option is selected, registered video frames
if instruct.registervid 
    if instruct.subpixel % Use image interpolation to correct for sub-pixel shifts
        interp_fact = 1;
    else
        interp_fact = 0;
    end
    
    imstack = imreg(imstack,interp_fact,instruct.fshift); % Custom image registration function
    z = size(imstack,1);
    
    %% If option in selected, create registered video in AVI format
    if instruct.createvid
        maxv = max(imstack(:));
        tempimg = uint8(squeeze(64*(imstack(1,:,:))/maxv)); % Represent first frame as 8-bit unsigned int. values
        
        regfilename = strcat(instruct.filename(1:end-4),'_registered.avi'); % Name registered movie
        if exist(regfilename,'file') % If registered movie with same filename exists, delete it
            delete(regfilename)
        end
        
        %% Movie creation
        vidObj = VideoWriter(regfilename);
        open(vidObj);       
        image(tempimg)
        colormap(gray)
        axis tight
        axis off
        set(gca,'nextplot','replacechildren');
        
        currFrame = getframe;
        writeVideo(vidObj,currFrame);
                
        for n = 2:z
            tempimg = round(squeeze(64*(imstack(n,:,:))/maxv));
            image(tempimg)
            currFrame = getframe;
            writeVideo(vidObj,currFrame);
        end
        close(vidObj);
    end
end  
[z,r,c] = size(imstack); 

%% If option is selected, downsample image by binning
if instruct.downsample
    dsimgstack = zeros(z,floor(r/2),floor(c/2));
    for n = 1:z
        for p = 1:floor(r/2)
            for q = 1:floor(c/2)
                dsimgstack(n,p,q) = sum(sum(imstack(n,2*p-1:2*p,2*q-1:2*q)))/4;
            end
        end
    end     
    [z,r,c] = size(dsimgstack);
    imstack = dsimgstack;
    clear dsimstack
    scale_factor = 2*scale_factor;
end   

%% Correct for fluctuations in overall image intensity 
mean_img = squeeze(mean(imstack,1)); 
zimg = zeros(z,r,c);
for n = 1:z
    zimg(n,:,:) = (squeeze(imstack(n,:,:)) - mean_img)./mean_img;
end

%% Calculate 'togetherimg' weighting factor (degree pixels fluctuate coherently with surrounding pixels) 
togetherimg = zeros(r,c);
for m = 1:r
    clc
    disp(['Processing Flow Data: ' num2str(round(100*m/r)) '%'])
    for n = 1:c
        nvec = zeros(1,8);
        if (m > 1)
            nvec(1,2) = sum((zimg(:,m,n).*zimg(:,m-1,n)) > 0);
            if (n > 1)
                nvec(1,1) = sum((zimg(:,m,n).*zimg(:,m-1,n-1)) > 0);
            end
            if (n < c)
                nvec(1,3) = sum((zimg(:,m,n).*zimg(:,m-1,n+1)) > 0);
            end
        end

        if (m < r)
            nvec(1,7) = sum((zimg(:,m,n).*zimg(:,m+1,n)) > 0);
            if (n > 1)
                nvec(1,6) = sum((zimg(:,m,n).*zimg(:,m+1,n-1) > 0));
            end
            if (n < c)
                nvec(1,8) = sum((zimg(:,m,n).*zimg(:,m+1,n+1)) > 0); 
            end
        end

        if (n > 1)
            nvec(1,4) = sum((zimg(:,m,n).*zimg(:,m,n-1)) > 0);
        end
        if (n < c)
            nvec(1,5) = sum((zimg(:,m,n).*zimg(:,m,n+1)) > 0);
        end                    

        togetherimg(m,n) = sum(nvec);
    end
end

togetherimginner = togetherimg(2:r-1,2:c-1);
togetherimg = togetherimg - min(togetherimginner(:));
togetherimg(togetherimg < 0) = 0;

%% If option is selected, use togetherimg to pre-mask in order to increase speed
if instruct.premask
    
    globalvarmax = 2.0;
    togetherimg(togetherimg > globalvarmax*std(togetherimg(:)) + mean(togetherimg(:))) =...
        globalvarmax*std(togetherimg(:)) + mean(togetherimg(:));

    mask_img = (togetherimg > (1/3)*max(togetherimg(:)));
else
    mask_img = ones(r,c);
end


%% Divides processing over blocks in order to conserve memory usage
movement_matrix = zeros(r,c,3);
meanpixval = squeeze(sum(imstack,1)/z);
multipartfactor = ceil(sqrt(r*c)/200);
overlap = 75;

for n = 1:multipartfactor
    for m = 1:multipartfactor

        r1 = max([1 round(r*(m-1)/multipartfactor)+1-(overlap-1)]);
        r2 = min([r round(r*m/multipartfactor)+(overlap-1)]);
        c1 = max([1 round(c*(n-1)/multipartfactor)+1-(overlap-1)]);
        c2 = min([c round(c*n/multipartfactor)+(overlap-1)]);

        r1w = round(r*(m-1)/multipartfactor)+1;
        r2w = round(r*m/multipartfactor);
        c1w = round(c*(n-1)/multipartfactor)+1;
        c2w = round(c*n/multipartfactor);

        startstop(1) = ((round(r*(m-1)/multipartfactor)+1) - r1) + 1;
        startstop(2) = (startstop(1) - 1) + (round(r*m/multipartfactor) - round(r*(m-1)/multipartfactor));
        startstop(3) = ((round(c*(n-1)/multipartfactor)+1) - c1) + 1;
        startstop(4) = (startstop(3) - 1) + (round(c*n/multipartfactor) - round(c*(n-1)/multipartfactor));       

        subimg = imstack(:,r1:r2,c1:c2);
        submeanpixval = meanpixval(r1:r2,c1:c2);
        submask = mask_img(r1:r2,c1:c2);     
             
        submatrix = flow_analysis(subimg,submeanpixval,startstop,instruct.downsample,submask,mask_exclusion);
        %           ^^^^^^^^^^^^^
        % This is the function that does the mapping operations
        
        movement_matrix(r1w:r2w,c1w:c2w,:) = movement_matrix(r1w:r2w,c1w:c2w,:) + submatrix;

    end
end 

movement_matrix = scale_factor*movement_matrix;

    
    
    
    
    
