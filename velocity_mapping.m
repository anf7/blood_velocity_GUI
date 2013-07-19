function velocity_mapping(instruct)

framerate = 21.4; %(frames/sec)
framesize_pix = 512; %(pixels)
framesize_um = 830; %(microns)


create_reg_video = false;
labeled_cells = false;
downsampleimg = false;



% scale_factor = framerate*((framesize_um/1000)/framesize_pix);
% scale_factor = 1; % Uncalibrated Data
% scale_factor = 0.322; % Primate Brain
% scale_factor = 0.1375; % Rat Lung - 296 Scope - 4x4 binning
scale_factor = 0.0472; % Thies's Scope - 5X
scale_factor = 2*0.0472; % Thies's Scope - 2.5X

start_dir = 'L:\All Staff\Andrew Fontanella\Flow Phantom';


if downsampleimg
    scale_factor = scale_factor*2;
end

registration_lags = 0;
interpolation_factor = 0;


homedir = pwd;

buttonval = 'Yes';
startind = 1;
while strcmp(buttonval,'Yes')
    
    [tdirlist, ttopdir] = find_subdirectories(start_dir);
    if strcmp(ttopdir(end),'\')
        ttopdir = ttopdir(1:end-1);
    end

    for n = startind:startind+length(tdirlist)-1
        dirlist(n).topdir = ttopdir;
        dirlist(n).dir = tdirlist(n-startind+1).dir;
        dirlist(n).path = tdirlist(n-startind+1).path;
    end   
        
    startind = startind + length(tdirlist);
    
    buttonval = questdlg('Add another directory?', ...
        'Input Selection...','Yes', 'No', 'No');

end



for i = 1:length(dirlist)

    currentpath = (strcat(dirlist(i).path,'\',dirlist(i).dir));
    if ~strcmp(currentpath(end),'\')
        currentpath = strcat(currentpath,'\');
    end
    cd(currentpath)
    tifflist = dir('*.tif');
    if isempty(tifflist)
        continue
    end
    
    
    savename = strcat(currentpath(length(dirlist(i).topdir)+1:end),'_compiled_data.mat');
    if strcmp(savename(1),'\')
        savename = savename(2:end);
    end
    if strcmp(savename(1),'_')
        savename = savename(2:end);
    end
    
    savename = savename + 0;
    savename(savename == 92) = 45;
    savename = char(savename);
        
    cd(homedir)
    
    
%% Process Raw Data
   
           
    for j = 1:length(tifflist)
        velocitydata = [];
        fulltiffname = strcat(currentpath,tifflist(j).name);
        savename2 = strcat(savename(1:end-18),tifflist(j).name(1:end-4),'_processed.mat');
        cd(currentpath)
        if ~isempty(dir('savename2'))
            continue
        end
        cd(homedir)

        tiffinfo = imfinfo(fulltiffname);
        if ~length(tiffinfo) > 1 || (strcmp(fulltiffname(end-14:end),'_registered.tif'))
            continue
        end


        [movement_matrix, togetherimg] = load_videos(tifflist(j).name,...
            currentpath,scale_factor,registration_lags,interpolation_factor,labeled_cells,create_reg_video,downsampleimg);
        vel_img = squeeze(movement_matrix(:,:,3));                   
        mask_img = logical(vel_img);


        velocitydata.movement_matrix = movement_matrix;
        velocitydata.vesselmask = mask_img;
        velocitydata.maxvelocity = max(vel_img(:));
        velocitydata.filename = tifflist(j).name;
        velocitydata.filepath = currentpath;
        velocitydata.scalefactor = scale_factor;
        velocitydata.togetherimg = togetherimg;


        cd(currentpath)
        save(savename2,'velocity_data')
        cd(homedir)

    end      
end  

cd(homedir)
    


     

