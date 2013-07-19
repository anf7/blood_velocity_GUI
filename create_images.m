function avgvelocity = create_images

framerate = 24; %(frames/sec)
framesize_pix = 1392; %(pixels)
framesize_um = 5300; %(microns)


framerate = 21.4; %(frames/sec)
framesize_pix = 512; %(pixels)
framesize_um = 830; %(microns)


% rescale_velocity_value = framerate*((framesize_um/1000)/framesize_pix);
rescale_velocity_value = 5;

% rescale_velocity_value = 0.0914/0.0944;
% rescale_velocity_value = 1;
% rescale_velocity_value = 1;
% rescale_velocity_value = 2*0.0472;
usecolorbars = true;
register_imagepairs = false;
combine_images = false;

% rescale_velocity_value = 4.272;
% rescale_velocity_value = 0.1375; % Rat Lung - 296 Scope

homedir = pwd;
fluor_channel_base = 0.25;

H = fspecial('gaussian',10,1);

indcolormap = gray(256);
jmap = jet(325);
RGBcolormap = jmap(56:311,:);
% redmap = gray(256);
% redmap(:,2) = zeros(256,1);
% redmap(:,3) = zeros(256,1);

usesubfolders = false;

%% Hb Parameters

create_hbsatcolorimg = 0;
create_hbsatwithflourcolorimg = 0;
create_hbsatimg = 0;
create_hbtotimg = 0;
create_hbsatmaskedimg = 0;
create_hbtotmaskedimg = 0;
create_parameterimg = 0;
create_diameterimg = 0;
create_rfpimg = 0;

minhbtot = 0;
maxhbsat = 75;
minhbsat = 0;

%% Flow Parameters
create_flowdirectionmap = 1;
create_flowdirectionmap_overlay = 0;
create_flowspeedmap = 1;
create_shearstressmap = 0;

use_log_scaling = 1;
log_range = 2; %% Set to 2
together_weighting = 1;
together_scale_factor = 1;
mintogethervalue = 0;
eliminate_noise = false;

autoscale_speed = false;
manual_maxspeed = 10;
minspeed_threshold = 0;
maxspeed_threshold = inf;

maxss = 7.5;

%%


nodata = true;

if register_imagepairs
    [basefilename,pathname] = uigetfile('*.mat','Select .mat file for base images...','MultiSelect','off');
    cd(pathname)
    [tformfilename,pathname] = uigetfile('*.mat','Select .mat files for transformed images...','MultiSelect','off');
    cd(homedir)
    
    filename = cell(1,2);
    filename{1} = basefilename;
    filename{2} = tformfilename;
else
    [filename,pathname] = uigetfile('*.mat','Select .mat files from which to generate images...','MultiSelect','on');
end

if ~iscell(filename)
    tempfname = filename;
    filename = cell(1);
    filename{1} = tempfname;
end

mkdir(pathname,'Images')


if combine_images
    for i = 1:length(filename)
        load(strcat(pathname,filename{i}));
        endpt = length(flowdata);
        for j = 1:endpt
            if ~isempty(strfind(filename{i},flowdata(j).filename(1:end-4)))
                mm_p = flowdata(j).movement_matrix;
                t_p = flowdata(j).togetherimg;
                
                if i == 1
                    [rrr,ccc,~] = size(mm_p);
                elseif i > 1 && (rrr < size(mm_p,1) || ccc < size(mm_p,2))
                    mm_p2 = mm_p(1:rrr,1:ccc,3);
                    t_p2 = t_p(1:rrr,1:ccc);
                    mm_p = mmp2;
                    t_p = t_p2;
                elseif i > 1 && (rrr > size(mm_p,1) || ccc > size(mm_p,2))
                    mm_p2 = zeros(rrr,ccc,3);
                    t_p2 = zeros(rrr,ccc);
                    mm_p2(1:size(mm_p,1),1:size(mm_p,2),3) = mm_p;
                    t_p2(1:size(mm_p,1),1:size(mm_p,2)) = t_p;
                    mm_p = mmp2;
                    t_p = t_p2;
                end
                
                mm(:,:,:,i) = mm_p;
                together(:,:,i) = t_p;
                continue
            end
        end
    end
    
    combmm = zeros(size(mm,1),size(mm,2),size(mm,3));
    combt = zeros(size(together,1),size(together,2));
    for i = 1:size(mm,4)
        for m = 1:size(mm,1)
            for n = 1:size(mm,2)
                if isequal(combmm(m,n,:),zeros(1,1,size(mm,3)))
                    combmm(m,n,:) = squeeze(mm(m,n,:,i));
                end
            end
        end
    end
    combt = mean(together,3);
    
    flowdata(1).movement_matrix = combmm;
    flowdata(1).togetherimg = combt;
    
    x = flowdata(1);
    flowdata = x;
    
    newname = strcat(filename{1}(1:end-4),'_combined.mat');
    cd(pathname)
    save(newname,'flowdata')
    cd(homedir)
    clear filename
    filename{1} = newname;
end



for n = 1:length(filename)   
    
    clc
    disp(['Working on ' num2str(n) ' of ' num2str(length(filename))])
    
    load(strcat(pathname,filename{n}));
    if usesubfolders
        mkdir(strcat(pathname,'Images'),filename{n}(1:end-4))
    end
    
    
    
       
    if exist('hbdata','var') && ~isempty(hbdata)
        nodata = false;
        
        if create_hbsatcolorimg 
            for m = 1:length(hbdata)
                imagefilename = strcat(filename{n}(1:end-4),'_',hbdata(m).filename(1:end-4),'_hbsatimg.tif'); 
                HbSat_img_slice = (255*(squeeze(hbdata(m).hbsatimg) - minhbsat)/(maxhbsat-minhbsat));
                HbSat_img_slice = imfilter(HbSat_img_slice,H);
                HbSat_img_slice = uint8(HbSat_img_slice); 
                HbTotfactor = 255*(squeeze(hbdata(m).hbtotimg) - minhbtot)/(max(max(squeeze(hbdata(m).hbtotimg))) - minhbtot);
                HbTotfactor3d(:,:,1) = HbTotfactor;
                HbTotfactor3d(:,:,2) = HbTotfactor;
                HbTotfactor3d(:,:,3) = HbTotfactor;
                RGBhbsat = ind2rgb(HbSat_img_slice,RGBcolormap);
                HbSatRGB = uint8(double(RGBhbsat).*double(HbTotfactor3d));
                if register_imagepairs && n == 1
                    baseimg(:,:,:,m) = HbSatRGB;
                    bfilename{m} = imagefilename;
                end     
                if register_imagepairs && n == 2
                    load('L:\All Staff\Andrew Fontanella\Hemodynamics\registration_pairs.mat')
%                     [xyinput_out, xybase_out] = cpselect(HbSatRGB,squeeze(baseimg(:,:,:,m)),'Wait',true); 
                    imgtform = cp2tform(xyinput_out,xybase_out,'lwm');
                    HbSatRGB = imtransform(HbSatRGB,imgtform,...
                        'XData',[1 1392], 'YData',[1 1040]);
                end                
                
            
                f = figure('visible','off');
                colormap(RGBcolormap)
                image(HbSatRGB)
                a = gca;
                axis(a,'off')
                axis(a,'image')
                set(gca,'position',[0 0 1 1],'units','normalized')
                if usecolorbars
                    colorbar('peer',a);
                    colorbar('YTickLabel',{num2str(round(minhbsat)),num2str(round(0.25*(maxhbsat-minhbsat)+minhbsat)),...
                        num2str(round(0.5*(maxhbsat-minhbsat)+minhbsat)),num2str(round(0.75*(maxhbsat-minhbsat)+minhbsat)),...
                        num2str(round(maxhbsat))},'YTick',[1 64 128 192 255])
                end                
                if usesubfolders
                    cd(strcat(pathname,'Images\',filename{n}(1:end-4)))
                    saveimgpath = strcat(pathname,'Images\',filename{n}(1:end-4));
                else
                    cd(strcat(pathname,'Images'))
                    saveimgpath = strcat(pathname,'Images');
                end
                saveas(f,imagefilename)
                cd(homedir)
                close(f) 
            end
        end
         
        if create_hbsatwithflourcolorimg 
            for m = 1:1%length(hbdata)
%                 if exist('hbdata(m).hbsatimg')
                    imagefilename = strcat(filename{n}(1:end-4),'_',hbdata(m).filename(1:end-4),'_hbsatgfpimg.tif');
                    HbSat_img_slice = uint8(255*(squeeze(hbdata(m).hbsatimg) - minhbsat)/(maxhbsat-minhbsat));
                    HbTotfactor = 255*(squeeze(hbdata(m).hbtotimg) - minhbtot)/(max(max(squeeze(hbdata(m).hbtotimg))) - minhbtot);
%                     HbTotfactor = 255*squeeze(hbdata(m).hbtotimg)/max(max(squeeze(hbdata(m).hbtotimg)));
                    HbTotfactor3d(:,:,1) = HbTotfactor;
                    HbTotfactor3d(:,:,2) = HbTotfactor;
                    HbTotfactor3d(:,:,3) = HbTotfactor;
                    RGBhbsat = ind2rgb(HbSat_img_slice,RGBcolormap);
                    HbSatRGB = uint8(double(RGBhbsat).*double(HbTotfactor3d));
                    
                    q = hbdata(m).associatedgfpimg;
                    
%                     fluor_img_slice = uint8(255*(1-fluor_channel_base)*squeeze(fluordata1(q).signalimg)/fluordata1(q).maxsignal);
                    fluor_img_slice = uint8(255*(1-fluor_channel_base)*squeeze(fluordata1(q).signalimg)/10);
                    Combined_img = HbSatRGB;
                    Combined_img(:,:,1) = Combined_img(:,:,1) + fluor_img_slice;
                    Combined_img(:,:,2) = Combined_img(:,:,2);% + fluor_img_slice;
                    Combined_img(:,:,3) = Combined_img(:,:,3) + fluor_img_slice;
                    
                    if register_imagepairs && n == 2
                        Combined_img = imtransform(Combined_img,imgtform,...
                            'XData',[1 1392], 'YData',[1 1040]);
                    end       
                    
                    
                    
                    f = figure('visible','off');
                    colormap(RGBcolormap)
                    image(Combined_img)
                    a = gca;
                    axis(a,'off')
                    axis(a,'image')
                    set(gca,'position',[0 0 1 1],'units','normalized')
                    if usecolorbars
                        colorbar('peer',a);
                        colorbar('YTickLabel',{num2str(round(minhbsat)),num2str(round(0.25*(maxhbsat-minhbsat)+minhbsat)),...
                            num2str(round(0.5*(maxhbsat-minhbsat)+minhbsat)),num2str(round(0.75*(maxhbsat-minhbsat)+minhbsat)),...
                            num2str(round(maxhbsat))},'YTick',[1 64 128 192 255])
                    end
                    if usesubfolders
                        cd(strcat(pathname,'Images\',filename{n}(1:end-4)))
                    else
                        cd(strcat(pathname,'Images'))
                    end
                    saveas(f,imagefilename)
                    cd(homedir)
                    close(f) 
%                 end
            end
        end
               
        if create_hbsatimg
            for m = 1:length(hbdata)
                imagefilename = strcat(filename{n}(1:end-4),'_',hbdata(m).filename(1:end-4),'_hbsatraw.tif');
                
                hbsatimg = hbdata(m).hbsatimg;
                
                if register_imagepairs && n == 2
                    hbsatimg = imtransform(hbsatimg,imgtform,...
                        'XData',[1 1392], 'YData',[1 1040]);
                end     
                
                f = figure('visible','off');
                colormap(indcolormap)
                imagesc(hbsatimg),caxis([minhbsat maxhbsat])
                a = gca;
                axis(a,'off')
                axis(a,'image')
                set(gca,'position',[0 0 1 1],'units','normalized')
                if usecolorbars
                    colorbar
                end
                if usesubfolders
                    cd(strcat(pathname,'Images\',filename{n}(1:end-4)))
                else
                    cd(strcat(pathname,'Images'))
                end
                saveas(f,imagefilename)
                cd(homedir)
                close(f) 
            end
        end
                
        if create_rfpimg
            for m = 1:length(fluordata2)
                imagefilename = strcat(filename{n}(1:end-4),'_',fluordata2(m).filename(1:end-4),'_rfp.tif');
                rfpimg = fluordata2(m).signalimg;
                f = figure('visible','off');
                colormap(gray)
                imagesc(rfpimg)
                a = gca;
                axis(a,'off')
                axis(a,'image')
                set(gca,'position',[0 0 1 1],'units','normalized')
                if usesubfolders
                    cd(strcat(pathname,'Images\',filename{n}(1:end-4)))
                else
                    cd(strcat(pathname,'Images'))
                end
                saveas(f,imagefilename)
                cd(homedir)
                close(f) 
            end
        end
        
        if create_hbtotimg
            for m = 1:length(hbdata)
                imagefilename = strcat(filename{n}(1:end-4),'_',hbdata(m).filename(1:end-4),'_hbtotraw.tif');
                hbtotimg = hbdata(m).hbtotimg;
                if register_imagepairs && n == 2
                    hbtotimg = imtransform(hbtotimg,imgtform,...
                        'XData',[1 1392], 'YData',[1 1040]);
                end  
                f = figure('visible','off');
                colormap(indcolormap)
                imagesc(hbtotimg)
                a = gca;
                axis(a,'off')
                axis(a,'image')
                set(gca,'position',[0 0 1 1],'units','normalized')
                if usecolorbars
                    colorbar
                end
                if usesubfolders
                    cd(strcat(pathname,'Images\',filename{n}(1:end-4)))
                else
                    cd(strcat(pathname,'Images'))
                end
                saveas(f,imagefilename)
                cd(homedir)
                close(f) 
            end
        end
        if create_hbsatmaskedimg
            for m = 1:length(hbdata)
                imagefilename = strcat(filename{n}(1:end-4),'_',hbdata(m).filename(1:end-4),'_hbsatmasked.tif');
                hbsatmasked = hbdata(m).hbsatmaskedimg;
                if register_imagepairs && n == 2
                    hbsatmasked = imtransform(hbsatmasked,imgtform,...
                        'XData',[1 1392], 'YData',[1 1040]);
                end  
                f = figure('visible','off');
                colormap(indcolormap)
                imagesc(hbsatmasked)
                a = gca;
                axis(a,'off')
                axis(a,'image')
                set(gca,'position',[0 0 1 1],'units','normalized')
                if usecolorbars
                    colorbar
                end
                if usesubfolders
                    cd(strcat(pathname,'Images\',filename{n}(1:end-4)))
                else
                    cd(strcat(pathname,'Images'))
                end
                saveas(f,imagefilename)
                cd(homedir)
                close(f) 
            end
        end
        if create_hbtotmaskedimg
            for m = 1:length(hbdata)
                imagefilename = strcat(filename{n}(1:end-4),'_',hbdata(m).filename(1:end-4),'_hbtotmasked.tif');
                hbtotmasked = hbdata(m).hbtotmaskedimg;
                if register_imagepairs && n == 2
                    hbtotmasked = imtransform(hbtotmasked,imgtform,...
                        'XData',[1 1392], 'YData',[1 1040]);
                end  
                
                f = figure('visible','off');
                colormap(indcolormap)
                imagesc(hbtotmasked)
                a = gca;
                axis(a,'off')
                axis(a,'image')
                set(gca,'position',[0 0 1 1],'units','normalized')
                if usecolorbars
                    colorbar
                end
                if usesubfolders
                    cd(strcat(pathname,'Images\',filename{n}(1:end-4)))
                else
                    cd(strcat(pathname,'Images'))
                end
                saveas(f,imagefilename)
                cd(homedir)
                close(f) 
            end
        end
        if create_parameterimg
             for m = 1:length(hbdata)
                imagefilename = strcat(filename{n}(1:end-4),'_',hbdata(m).filename(1:end-4),'_parameterimg.tif');
                paramimg = hbdata(m).parameterimg_newroi;
                if register_imagepairs && n == 2
                    paramimg = imtransform(paramimg,imgtform,...
                        'XData',[1 1392], 'YData',[1 1040]);
                end  
                
                f = figure('visible','off');
                image(paramimg)
                a = gca;
                axis(a,'off')
                axis(a,'image')
                set(gca,'position',[0 0 1 1],'units','normalized')
                if usesubfolders
                    cd(strcat(pathname,'Images\',filename{n}(1:end-4)))
                else
                    cd(strcat(pathname,'Images'))
                end
                saveas(f,imagefilename)
                cd(homedir)
                close(f) 
            end
        end
        if create_diameterimg
            for m = 1:length(hbdata)
                imagefilename = strcat(filename{n}(1:end-4),'_',hbdata(m).filename(1:end-4),'_diameterimg.tif');
                daimimg = hbdata(m).diameterimg;
                if register_imagepairs && n == 2
                    daimimg = imtransform(daimimg,imgtform,...
                        'XData',[1 1392], 'YData',[1 1040]);
                end  
                f = figure('visible','off');
                colormap(indcolormap)
                imagesc(daimimg)
                a = gca;
                axis(a,'off')
                axis(a,'image')
                set(gca,'position',[0 0 1 1],'units','normalized')
                
                if usecolorbars
                    colorbar
                end
                if usesubfolders
                    cd(strcat(pathname,'Images\',filename{n}(1:end-4)))
                else
                    cd(strcat(pathname,'Images'))
                end
                saveas(f,imagefilename)
                cd(homedir)
                close(f) 
            end
        end
    end
    
    if exist('ssimg','var') && create_shearstressmap
        nodata = false;
        imagefilename = strcat(filename{n}(1:end-4),'_shearstressimg.tif');
        
        SE = strel('disk',4);
        mmax = max(ssimg(:));
        ssimg = imdilate(ssimg,SE);
        ssimg = mmax*ssimg/max(ssimg(:));
        
        gmap = gray(276);
        gmap256 = gmap(21:276,:);
        
        ssrgb = setbackgroundtoblack(ssimg,maxss,gmap256);
        
               
        if register_imagepairs && n == 2
            ssrgb = imtransform(ssrgb,imgtform,...
                'XData',[1 1392], 'YData',[1 1040]);
        end  
        
        f = figure('visible','off');
        colormap(gmap256)
        image(ssrgb)
        a = gca;
        axis(a,'off')
        axis(a,'image')
        set(gca,'position',[0 0 1 1],'units','normalized')
        maxss
        if usecolorbars
            colorbar('peer',a);
            colorbar('YTickLabel',{num2str(0),num2str(round(0.25*maxss)),...
                num2str(round(0.5*maxss)),num2str(round(0.75*maxss)),...
                num2str(round(maxss))},'YTick',[1 64 128 192 255])
        end

        if usesubfolders
            cd(strcat(pathname,'Images\',filename{n}(1:end-4)))
        else
            cd(strcat(pathname,'Images'))
        end
        saveas(f,imagefilename)
        cd(homedir)
        close(f) 
    end   
    
    if (exist('flowdata','var') && ~isempty(flowdata)) || exist('movement_matrix','var')
        nodata = false;
        if create_flowspeedmap 
            if ~exist('flowdata','var')
                
                velocityimg = rescale_velocity_value*squeeze(movement_matrix(:,:,3));
                
                avgvelocity(n) = sum(velocityimg(:))/nnz(velocityimg);
                                                   
                imagefilename = strcat('Speed_',filename{n},'.tif');
                if autoscale_speed
                    max_velocity = max(velocityimg(:));
                else
                    max_velocity = manual_maxspeed;
                end       
                
                if use_log_scaling
                    p_maxvelocity = max_velocity;
                    max_velocity = log10(max_velocity);
                    min_velocity = max_velocity-log_range;
                    p_maxspeed_threshold = maxspeed_threshold;
                    p_minspeed_threshold = minspeed_threshold;
                    maxspeed_threshold = log10(maxspeed_threshold);
                    minspeed_threshold = log10(minspeed_threshold);
                    velocityimg = log10(velocityimg);
                    velocityimg((velocityimg < max_velocity-log_range) & (~isinf(velocityimg))) = max_velocity-log_range;
                    
                    movement_matrix(:,:,3) = velocityimg;
                else
                    p_maxvelocity = max_velocity;
                    p_maxspeed_threshold = maxspeed_threshold;
                    p_minspeed_threshold = minspeed_threshold;
                    min_velocity = 0;
                    velocityimg(velocityimg > maxspeed_threshold) = 0;
                    velocityimg(velocityimg < minspeed_threshold) = 0;
                    velocityimg(velocityimg > max_velocity) = max_velocity;
                    movement_matrix(:,:,3) = velocityimg;
                end
                
                                                 
                [~, rgbvel, ~] = drawmaps(movement_matrix,min_velocity,max_velocity,minspeed_threshold,maxspeed_threshold);
                  
                
                
                if eliminate_noise
                    
                    rgbvel(:,:,1) = medfilt2(squeeze(rgbvel(:,:,1)));
                    rgbvel(:,:,2) = medfilt2(squeeze(rgbvel(:,:,2)));
                    rgbvel(:,:,3) = medfilt2(squeeze(rgbvel(:,:,3)));
                    
                    
%                     rgbmask = ~isequal(rgbimg,zeros(1,1,3));
                        
                    

%                     rgbvel(:,:,1) = imfilter(squeeze(rgbvel(:,:,1)),H).*rgbmask;
%                     rgbvel(:,:,2) = imfilter(squeeze(rgbvel(:,:,2)),H).*rgbmask;
%                     rgbvel(:,:,3) = imfilter(squeeze(rgbvel(:,:,3)),H).*rgbmask;
                    
                    rgbvel(:,:,1) = imfilter(squeeze(rgbvel(:,:,1)),H);
                    rgbvel(:,:,2) = imfilter(squeeze(rgbvel(:,:,2)),H);
                    rgbvel(:,:,3) = imfilter(squeeze(rgbvel(:,:,3)),H);

                end
                
                if together_weighting && exist('togetherimg','var')
                    
                     togetherimg(togetherimg < mintogethervalue) = 0;
                     togetherimg(togetherimg <= 0) = NaN;
                     togetherimg = togetherimg - min(min(togetherimg(2:end-1,2:end-1)));
                     togetherimg(isnan(togetherimg)) = 0;
                     togetherimg(togetherimg < 0) = 0;
                     togetherimg = together_scale_factor*(togetherimg/max(togetherimg(:)));
                     togetherimg(togetherimg > 1) = 1;
                     
                     
                     rgbvel(:,:,1) = (squeeze(rgbvel(:,:,1)).*togetherimg);
                     rgbvel(:,:,2) = (squeeze(rgbvel(:,:,2)).*togetherimg);
                     rgbvel(:,:,3) = (squeeze(rgbvel(:,:,3)).*togetherimg);
                     
                end
                
                rgbvel = rgbvel/max(rgbvel(:));

                if register_imagepairs && n == 2
                    rgbvel = imtransform(rgbvel,imgtform,...
                        'XData',[1 1392], 'YData',[1 1040]);
                end  
                
                rgbvel = rgbvel/max(rgbvel(:));
                
                figure('visible','off')
                image(rgbvel), axis off, axis image
                colormap(jet(128))
                set(gca,'position',[0 0 1 1],'units','normalized')
                if usecolorbars
                    if use_log_scaling
                        firsttic = floor(max_velocity);

                        tic_pos = zeros(1,26);
                        tic_pos(18) = 10^firsttic;
                        tic_pos(9) = 10^(firsttic-1);
                        for tt = 19:26
                            tic_pos(tt) = 10^firsttic + 10^firsttic*(tt-18);
                        end
                        for tt = 10:17
                            tic_pos(tt) = 10^(firsttic-1) + 10^(firsttic-1)*(tt-9);
                        end
                        for tt = 1:8
                            tic_pos(tt) = 10^(firsttic-2) + 10^(firsttic-2)*(tt);
                        end
                        tic_pos = log10(tic_pos);

                        tic_pos = (127*((tic_pos - (max_velocity-log_range))/((max_velocity)-((max_velocity-log_range)))));

                        tic_label = {'','','','','','','','',num2str(10^(firsttic-1),'%0.3g'),...
                            '','','','','','','','',num2str(10^(firsttic),'%0.3g'),'','','','','','','',''};
                        colorbar('YTickLabel',tic_label,'YTick',tic_pos)

                    else 
                        colorbar('YTickLabel',...
                            {'0',num2str((1/6)*max_velocity),num2str((2/6)*max_velocity),num2str((3/6)*max_velocity),...
                             num2str((4/6)*max_velocity),num2str((5/6)*max_velocity),num2str(max_velocity)},'YTick',...
                             [1 round(128/6) round(2*128/6) round(3*128/6) round(4*128/6) round(5*128/6) round(6*128/6)])
                    end
                end
                max_velocity = p_maxvelocity;
                maxspeed_threshold = p_maxspeed_threshold;
                minspeed_threshold = p_minspeed_threshold;
                    
                figh = gcf;
                if usesubfolders
                    cd(strcat(pathname,'Images\',filename{n}(1:end-4)))
                else
                    cd(strcat(pathname,'Images'))
                end
                saveas(figh,imagefilename)
                cd(homedir)
                close(figh) 
                
                
%                 bandwlinearimg = rescale_velocity_value*squeeze(movement_matrix_p(:,:,3));
%                 bandwimagefilename = strcat(filename{n},'_flowspeed_BandW.tif');
%                 f = figure('visible','off');
%                 imagesc(bandwlinearimg),colormap(gray(2^16)),colorbar,caxis([0 max_velocity])
%                 axis off
%                 axis image
%                 set(gca,'position',[0 0 1 1],'units','normalized')
%                 if usesubfolders
%                     cd(strcat(pathname,'Images\',filename{n}(1:end-4)))
%                 else
%                     cd(strcat(pathname,'Images'))
%                 end
%                 saveas(f,bandwimagefilename)
%                 cd(homedir)
%                 close(figh)
                
            else
                
                for m = 1:length(flowdata)
                    imagefilename = strcat('Speed_',filename{n}(1:end-4),'_',flowdata(m).filename(1:end-4),'.tif');
                    movement_matrix = flowdata(m).movement_matrix;
                    movement_matrix_p = movement_matrix;
                    velocityimg = rescale_velocity_value*squeeze(movement_matrix(:,:,3));
                    
                    avgvelocity(n) = sum(velocityimg(:))/nnz(velocityimg);
                    
                    if autoscale_speed
                        max_velocity = max(velocityimg(:));
                    else
                        max_velocity = manual_maxspeed;
                    end     
                    
                    if use_log_scaling
                        p_maxvelocity = max_velocity;
                        max_velocity = log10(max_velocity);
                        min_velocity = max_velocity-log_range;
                        p_maxspeed_threshold = maxspeed_threshold;
                        p_minspeed_threshold = minspeed_threshold;
                        maxspeed_threshold = log10(maxspeed_threshold);
                        minspeed_threshold = log10(minspeed_threshold);
                        velocityimg = log10(velocityimg);
                        velocityimg((velocityimg < max_velocity-log_range) & (~isinf(velocityimg))) = max_velocity-log_range;

                        movement_matrix(:,:,3) = velocityimg;
                    else
                        p_maxvelocity = max_velocity;
                        p_maxspeed_threshold = maxspeed_threshold;
                        p_minspeed_threshold = minspeed_threshold;
                        min_velocity = 0;
                        velocityimg(velocityimg > maxspeed_threshold) = 0;
                        velocityimg(velocityimg < minspeed_threshold) = 0;
                        velocityimg(velocityimg > max_velocity) = max_velocity;
                        movement_matrix(:,:,3) = velocityimg;
                    end              

                    [~, rgbvel, ~] = drawmaps(movement_matrix,min_velocity,max_velocity,minspeed_threshold,maxspeed_threshold);
                    
                                       
                    if eliminate_noise
                        
%                         rgbmask = ~isequal(rgbvel,zeros(1,1,3));
                        
                        rgbvel(:,:,1) = medfilt2(squeeze(rgbvel(:,:,1)));
                        rgbvel(:,:,2) = medfilt2(squeeze(rgbvel(:,:,2)));
                        rgbvel(:,:,3) = medfilt2(squeeze(rgbvel(:,:,3)));
                        
%                         rgbvel(:,:,1) = imfilter(squeeze(rgbvel(:,:,1)),H).*rgbmask;
%                         rgbvel(:,:,2) = imfilter(squeeze(rgbvel(:,:,2)),H).*rgbmask;
%                         rgbvel(:,:,3) = imfilter(squeeze(rgbvel(:,:,3)),H).*rgbmask;
                        
                        rgbvel(:,:,1) = imfilter(squeeze(rgbvel(:,:,1)),H);
                        rgbvel(:,:,2) = imfilter(squeeze(rgbvel(:,:,2)),H);
                        rgbvel(:,:,3) = imfilter(squeeze(rgbvel(:,:,3)),H);
                       
                    end
                    
                    if together_weighting
                        
                        togetherimg = flowdata(m).togetherimg;
                        
                        togetherimg(togetherimg < mintogethervalue) = 0;
                        togetherimg(togetherimg <= 0) = NaN;
                        togetherimg = togetherimg - min(min(togetherimg(2:end-1,2:end-1)));
                        togetherimg(isnan(togetherimg)) = 0;
                        togetherimg(togetherimg < 0) = 0;
                        togetherimg = together_scale_factor*(togetherimg/max(togetherimg(:)));
                        togetherimg(togetherimg > 1) = 1;
                        
                        rgbvel(:,:,1) = (squeeze(rgbvel(:,:,1)).*togetherimg);
                        rgbvel(:,:,2) = (squeeze(rgbvel(:,:,2)).*togetherimg);
                        rgbvel(:,:,3) = (squeeze(rgbvel(:,:,3)).*togetherimg);

                    end
                    
                    rgbvel = rgbvel/max(rgbvel(:));
                    
                    if register_imagepairs && n == 2
                        rgbvel = imtransform(rgbvel,imgtform,...
                            'XData',[1 1392], 'YData',[1 1040]);
                    end  
                    rgbvel = rgbvel/max(rgbvel(:));

                    figure('visible','off')
                    image(rgbvel), axis off, axis image
                    colormap(jet(128))
                    set(gca,'position',[0 0 1 1],'units','normalized')

                    if usecolorbars
                        if use_log_scaling
                            firsttic = floor(max_velocity);

                            tic_pos = zeros(1,26);
                            tic_pos(18) = 10^firsttic;
                            tic_pos(9) = 10^(firsttic-1);
                            for tt = 19:26
                                tic_pos(tt) = 10^firsttic + 10^firsttic*(tt-18);
                            end
                            for tt = 10:17
                                tic_pos(tt) = 10^(firsttic-1) + 10^(firsttic-1)*(tt-9);
                            end
                            for tt = 1:8
                                tic_pos(tt) = 10^(firsttic-2) + 10^(firsttic-2)*(tt);
                            end
                            tic_pos = log10(tic_pos);

                            tic_pos = (127*((tic_pos - (max_velocity-log_range))/((max_velocity)-((max_velocity-log_range)))));

                            tic_label = {'','','','','','','','',num2str(10^(firsttic-1),'%0.3g'),...
                                '','','','','','','','',num2str(10^(firsttic),'%0.3g'),'','','','','','','',''};
                            colorbar('YTickLabel',tic_label,'YTick',tic_pos)

                        else 
                            colorbar('YTickLabel',...
                                {'0',num2str((1/6)*max_velocity),num2str((2/6)*max_velocity),num2str((3/6)*max_velocity),...
                                 num2str((4/6)*max_velocity),num2str((5/6)*max_velocity),num2str(max_velocity)},'YTick',...
                                 [1 round(128/6) round(2*128/6) round(3*128/6) round(4*128/6) round(5*128/6) round(6*128/6)])
                        end
                    end
                    max_velocity = p_maxvelocity;
                    maxspeed_threshold = p_maxspeed_threshold;
                    minspeed_threshold = p_minspeed_threshold;

                    figh = gcf;
                    if usesubfolders
                        cd(strcat(pathname,'Images\',filename{n}(1:end-4)))
                    else
                        cd(strcat(pathname,'Images'))
                    end
                    saveas(figh,imagefilename)
                    cd(homedir)
                    close(figh) 
                    
%                     bandwlinearimg = rescale_velocity_value*squeeze(movement_matrix_p(:,:,3));
%                     bandwimagefilename = strcat(filename{n}(1:end-4),'_',flowdata(m).filename(1:end-4),'_flowspeed_BandW.tif');
%                     f = figure('visible','off');
%                     imagesc(bandwlinearimg),colormap(gray(2^16)),colorbar,caxis([0 max_velocity])
%                     axis off 
%                     axis image
%                     set(gca,'position',[0 0 1 1],'units','normalized')
%                     if usesubfolders
%                         cd(strcat(pathname,'Images\',filename{n}(1:end-4)))
%                     else
%                         cd(strcat(pathname,'Images'))
%                     end
%                     saveas(f,bandwimagefilename)
%                     cd(homedir)
%                     close(figh)
                    
                end
            end
        end
        
        if create_flowdirectionmap 
            
            if ~exist('flowdata','var')
                
                imagefilename = strcat('Direction_',filename{n},'.tif');
                [~, ~, dirimg] = drawmaps(movement_matrix,[],max_velocity,minspeed_threshold,maxspeed_threshold);
                
                               
                if eliminate_noise
                    
%                     dirmask = ~isequal(dirimg,zeros(1,1,3));
                    dirimg(:,:,1) = medfilt2(squeeze(dirimg(:,:,1)));
                    dirimg(:,:,2) = medfilt2(squeeze(dirimg(:,:,2)));
                    dirimg(:,:,3) = medfilt2(squeeze(dirimg(:,:,3)));
                    

%                     dirimg(:,:,1) = imfilter(squeeze(dirimg(:,:,1)),H).*dirmask;
%                     dirimg(:,:,2) = imfilter(squeeze(dirimg(:,:,2)),H).*dirmask;
%                     dirimg(:,:,3) = imfilter(squeeze(dirimg(:,:,3)),H).*dirmask;
                    
                    dirimg(:,:,1) = imfilter(squeeze(dirimg(:,:,1)),H);
                    dirimg(:,:,2) = imfilter(squeeze(dirimg(:,:,2)),H);
                    dirimg(:,:,3) = imfilter(squeeze(dirimg(:,:,3)),H);

                end
                
                if together_weighting && exist('togetherimg','var')
                    
                     togetherimg(togetherimg < mintogethervalue) = 0;
                     togetherimg(togetherimg <= 0) = NaN;
                     togetherimg = togetherimg - min(min(togetherimg(2:end-1,2:end-1)));
                     togetherimg(isnan(togetherimg)) = 0;
                     togetherimg(togetherimg < 0) = 0;
                     togetherimg = together_scale_factor*(togetherimg/max(togetherimg(:)));
                     togetherimg(togetherimg > 1) = 1;
                
                     dirimg(:,:,1) = (squeeze(dirimg(:,:,1)).*togetherimg);
                     dirimg(:,:,2) = (squeeze(dirimg(:,:,2)).*togetherimg);
                     dirimg(:,:,3) = (squeeze(dirimg(:,:,3)).*togetherimg);
                end
                
                dirimg = dirimg/max(dirimg(:));
                
                if create_flowdirectionmap_overlay
                    rawfilename = strcat(flowdata(m).filepath,flowdata(m).filename);
                    bright = double(imread(fullfilename,2));
                    bright = bright(1:size(dirimg,1),1:size(dirimg,2));
                    bright = max(bright(:))-bright;
                    bright = bright/max(bright(:));
                    dirimg(:,:,1) = dirimg(:,:,1).*bright;
                    dirimg(:,:,2) = dirimg(:,:,2).*bright;
                    dirimg(:,:,3) = dirimg(:,:,3).*bright;
                end

                
                
                if register_imagepairs && n == 2
                    dirimg = imtransform(dirimg,imgtform,...
                        'XData',[1 1392], 'YData',[1 1040]);
                end  
                
                dirimg = dirimg/max(dirimg(:));
                
                figure('visible','off')
                image(dirimg), axis off, axis image
                set(gca,'position',[0 0 1 1],'units','normalized')
                figh = gcf;
                if usesubfolders
                    cd(strcat(pathname,'Images\',filename{n}(1:end-4)))
                else
                    cd(strcat(pathname,'Images'))
                end
                saveas(figh,imagefilename)
                cd(homedir)
                close(figh) 
                
            else
            
                for m = 1:length(flowdata)
                    imagefilename = strcat('Direction_',filename{n}(1:end-4),'_',flowdata(m).filename(1:end-4),'.tif');
                    movement_matrix = flowdata(m).movement_matrix;

                    [~, ~, dirimg] = drawmaps(movement_matrix,[],max_velocity,minspeed_threshold,maxspeed_threshold);
                    
                                       
                    if eliminate_noise
                        
%                         dirmask = ~isequal(dirimg,zeros(1,1,3));
                    
                        dirimg(:,:,1) = medfilt2(squeeze(dirimg(:,:,1)));
                        dirimg(:,:,2) = medfilt2(squeeze(dirimg(:,:,2)));
                        dirimg(:,:,3) = medfilt2(squeeze(dirimg(:,:,3)));

%                         dirimg(:,:,1) = imfilter(squeeze(dirimg(:,:,1)),H).*dirmask;
%                         dirimg(:,:,2) = imfilter(squeeze(dirimg(:,:,2)),H).*dirmask;
%                         dirimg(:,:,3) = imfilter(squeeze(dirimg(:,:,3)),H).*dirmask;
                        
                        dirimg(:,:,1) = imfilter(squeeze(dirimg(:,:,1)),H);
                        dirimg(:,:,2) = imfilter(squeeze(dirimg(:,:,2)),H);
                        dirimg(:,:,3) = imfilter(squeeze(dirimg(:,:,3)),H);

                    end
                    
                    if together_weighting
                         togetherimg = flowdata(m).togetherimg;                  
                         
                         togetherimg(togetherimg < mintogethervalue) = 0;
                         togetherimg(togetherimg <= 0) = NaN;
                         togetherimg = togetherimg - min(min(togetherimg(2:end-1,2:end-1)));
                         togetherimg(isnan(togetherimg)) = 0;
                         togetherimg(togetherimg < 0) = 0;
                         togetherimg = together_scale_factor*(togetherimg/max(togetherimg(:)));
                         togetherimg(togetherimg > 1) = 1;
                         
                         
                         
                         dirimg(:,:,1) = (squeeze(dirimg(:,:,1)).*togetherimg);
                         dirimg(:,:,2) = (squeeze(dirimg(:,:,2)).*togetherimg);
                         dirimg(:,:,3) = (squeeze(dirimg(:,:,3)).*togetherimg);
                    end
                    
                    
                    dirimg = dirimg/max(dirimg(:));
                    
                    if create_flowdirectionmap_overlay
                        rawfilename = strcat(flowdata(m).filepath,flowdata(m).filename);
                        bright = double(imread(rawfilename,2));
                                               
                        bright = max(bright(:))-bright;
                        bright = bright - min(bright(:));
                        
                        bright = bright/max(bright(:));
                        
%                         bright = downsample(bright,2);
%                         bright = downsample(bright',2);
%                         bright = bright';
                            
                        if size(bright,1) > size(dirimg,1)
                            bright = bright(1:size(dirimg,1),:);
                        else
                            bright(size(bright,1):size(dirimg,1),:) = zeros(size(dirimg,1)-size(bright,1)+1,size(bright,2));
                        end
                        
                        if size(bright,2) > size(dirimg,2)
                            bright = bright(:,2:size(dirimg,2));
                        else
                            bright(:,size(bright,2):size(dirimg,2)) = zeros(size(bright,1),size(dirimg,2)-size(bright,2)+1);
                        end
                        
                        
                        for ii = 1:size(dirimg,1)
                            for jj = 1:size(dirimg,2)
                                if isequal(dirimg(ii,jj,:),zeros(1,1,3))
                                    dirimg(ii,jj,:) = ones(1,1,3);
                                end
                            end
                        end
                        
                        dirimg(:,:,1) = dirimg(:,:,1).*bright;
                        dirimg(:,:,2) = dirimg(:,:,2).*bright;
                        dirimg(:,:,3) = dirimg(:,:,3).*bright;
                        
                        
                    end

                    if register_imagepairs && n == 2
                        dirimg = imtransform(dirimg,imgtform,...
                            'XData',[1 1392], 'YData',[1 1040]);
                    end  
                    
                    dirimg = dirimg/max(dirimg(:));
                    
                    figure('visible','off')
                    image(dirimg), axis off, axis image
                    set(gca,'position',[0 0 1 1],'units','normalized')
                    figh = gcf;
                    if usesubfolders
                        cd(strcat(pathname,'Images\',filename{n}(1:end-4)))
                    else
                        cd(strcat(pathname,'Images'))
                    end
                    saveas(figh,imagefilename)
                    cd(homedir)
                    close(figh) 
                end
            end
        end
    end
end

if nodata
    clc
    disp('No data to process')
else
    clc
    disp('Done')
end



