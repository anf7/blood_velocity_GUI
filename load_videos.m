function [movement_matrix, togetherimg_p] = process_data(filename,filepath,scale_factor,registration_lags,interpolation_factor,labeled_cells,create_reg_video,downsampleimg)


%% Load multipage TIFFs
home = pwd;

if registration_lags > 0
    register_stack = true;
end

if nargin < 3
    [filename, filepath] = uigetfile('*.tif', 'Select multipage .TIF files...', 'MultiSelect', 'On');
end

if ~iscell(filename)
    tempfname = filename;
    filename = cell(1);
    filename{1} = tempfname;
end




for i = 1:length(filename);
    
    if isempty(matfile)
        clear imstack
        fullfilename = strcat(filepath,filename{i});

        info = imfinfo(fullfilename);
        r = info(1).Height;
        c = info(1).Width;
        fulllength = length(info);

        if ~(fulllength > 1)
            clc
            disp('Please select only multipage TIFFs')
            pause(5)
            clc
            continue
        end

        if usepow2
            z = 2^(floor(log2(fulllength)));
            startframe = 1 + round((fulllength - z)/2);
        else
            z = fulllength - 1;
            startframe = 2;
        end


        imstack = zeros(z,r,c);
        for n = 1:z
            clc
            disp(['Loading frame ' num2str(n) ' of ' num2str(z)])
            imstack(n,:,:) = imread(fullfilename,startframe+n-1);
        end
    else
        imstack = matfile;
        clear matfile
    end
    
  
    if register_stack
        
       imstack = imreg(imstack,interpolation_factor,registration_lags); 
        
        [z,~,~] = size(imstack);
        maxv = max(imstack(:));
        cd(filepath)
        regfilename = strcat(filename{i}(1:end-4),'_registered.tif');
        if exist(regfilename,'file')
            delete(regfilename)
        end
        if ~isempty(imstack) && create_reg_video
            for n = 1:z
                tempimg = uint8(squeeze(255*(imstack(n,:,:))/maxv));
                imwrite(tempimg, regfilename, 'WriteMode', 'append', 'Compression', 'none')
            end
        end
        cd(home)
    end
    
    [z,r,c] = size(imstack);
    if downsampleimg
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
    end   
    
    togetherimg = zeros(r,c);
    movement_matrix = [];
    
    if isempty(imstack)
        continue
    end
    
    
    if ~isa(imstack,'double')
        imstack = double(imstack);
    end
    
 
        
    mean_img = squeeze(mean(imstack,1)); 
    zimg = zeros(z,r,c);
    for n = 1:z
        zimg(n,:,:) = (squeeze(imstack(n,:,:)) - mean_img)./mean_img;
    end



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

    togetherimg_p = togetherimg;

    if ~labeled_cells

        togetherimginner = togetherimg(2:r-1,2:c-1);
        togetherimg = togetherimg - min(togetherimginner(:));
        togetherimg(togetherimg < 0) = 0;


        globalvarmax = 2.0;
        togetherimg(togetherimg > globalvarmax*std(togetherimg(:)) + mean(togetherimg(:))) =...
            globalvarmax*std(togetherimg(:)) + mean(togetherimg(:));


        mask_img = (togetherimg > (1/3)*max(togetherimg(:)));
    else
        mask_img = ones(r,c);

    end

    
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
            

            submatrix = flow_analysis(subimg,submeanpixval,startstop,downsampleimg,submask);
            
            movement_matrix(r1w:r2w,c1w:c2w,:) = movement_matrix(r1w:r2w,c1w:c2w,:) + submatrix;

        end
    end   
    movement_matrix = scale_factor*movement_matrix;
      
end
    
    
    
    
    
