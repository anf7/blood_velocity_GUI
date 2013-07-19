function showflowimages
homedir = pwd;


start_dir = 'C:\Users\anf7\Desktop\Radiation Study Flow Images';


ctform = makecform('cmyk2srgb');
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
    
    buttonval = questdlg('Add another top level directory?', ...
        'Input Selection...','Yes', 'No', 'No');

end



for i = 1:length(dirlist)    
    
    clc
    currentpath = (strcat(dirlist(i).path,'\',dirlist(i).dir))
    if ~strcmp(currentpath(end),'\')
        currentpath = strcat(currentpath,'\');
    end
    cd(currentpath)
    tifflist = dir('Direction*flow.tif');
    if isempty(tifflist) || length(tifflist) < 6
        continue
    end
    bfimgname = dir('bright.tif');
    if isempty(bfimgname)
        continue
    end
    
    
    savename = strcat(currentpath(length(dirlist(i).topdir)+1:end),'_compiled_flow.tif');
    if strcmp(savename(1),'\')
        savename = savename(2:end);
    end
    if strcmp(savename(1),'_')
        savename = savename(2:end);
    end
    
    savename = savename + 0;
    savename(savename == 92) = 45;
    savename = char(savename);
    

    cimg = cell(6,1);
    rr = cell(6,1);
    cc = cell(6,1);
    
    
    for j = 1:length(tifflist)
        fulltiffname = strcat(currentpath,tifflist(j).name);
        img = imread(fulltiffname);
        if size(img,3) == 4
%             img = double(img);
%             img = uint16(2^16*img/max(img(:)));
            img = applycform(img,ctform);
            img = double(img);
            img = img/max(img(:));
            img = 1-img;
            cutimg = img(30:end-30,30:end-30);
            img = img - min(cutimg(:));
            img = img/max(img(:));
        else
            img = double(img);
            img = img/max(img(:)); 
        end
        
        
        [r, c, ~] = size(img);
        
        if ~isempty(strfind(fulltiffname,'Pre'))
            cimg{1} = img;
            rr{1} = r;
            cc{1} = c;
        elseif ~isempty(strfind(fulltiffname,'Post'))
            cimg{2} = img;
            rr{2} = r;
            cc{2} = c;
        elseif ~isempty(strfind(fulltiffname,'Day1'))
            cimg{3} = img;
            rr{3} = r;
            cc{3} = c;
        elseif ~isempty(strfind(fulltiffname,'Day2'))
            cimg{4} = img;
            rr{4} = r;
            cc{4} = c;
        elseif ~isempty(strfind(fulltiffname,'Day3'))
            cimg{5} = img;
            rr{5} = r;
            cc{5} = c;
        elseif ~isempty(strfind(fulltiffname,'Day4'))
            cimg{6} = img;
            rr{6} = r;
            cc{6} = c;
        else
            continue
        end
    end
    bfimg = double(imread(strcat(currentpath,bfimgname(1).name)));
    bfimg = bfimg/max(bfimg(:));
    [bfr, bfc] = size(bfimg);
    bfimgc = cat(3,bfimg,bfimg);
    bfimgc = cat(3,bfimgc,bfimg);
    
    
    if isempty(cimg(1))
        cimg{1} = zeros(size(bfimgc));
    end
    if isempty(cimg(2))
        cimg{2} = zeros(size(cimg{1}));
    end
    if isempty(cimg(3))
        cimg{3} = zeros(size(cimg{2}));
    end
    if isempty(cimg(4))
        cimg{4} = zeros(size(cimg{3}));
    end
    if isempty(cimg(5))
        cimg{5} = zeros(size(cimg{4}));
    end
    if isempty(cimg(6))
        cimg{6} = zeros(size(cimg{5}));
    end
    
        
    combimg = bfimgc;
    combimg(bfr+1:bfr+rr{1},1:cc{1},:) = cimg{1};
    combimg(bfr+1:bfr+rr{2},cc{1}+1:cc{1}+cc{2},:) = cimg{2};
    combimg(bfr+1:bfr+rr{3},cc{1}+cc{2}+1:cc{1}+cc{2}+cc{3},:) = cimg{3};
    
    [rnew, cnew, ~] = size(combimg);
    
    combimg(rnew+1:rnew+rr{4},1:cc{4},:) = cimg{4};
    combimg(rnew+1:rnew+rr{5},cc{4}+1:cc{4}+cc{5},:) = cimg{5};
    combimg(rnew+1:rnew+rr{6},cc{4}+cc{5}+1:cc{4}+cc{5}+cc{6},:) = cimg{6};
    
    [rnew2, cnew2, ~] = size(combimg);
    
    combimg(rnew-2:rnew+2,:,:) = ones(5,cnew2,3);
    combimg(bfr+1:bfr+rr{1},cc{1}-2:cc{1}+2,:) = ones(rr{1},5,3);
    combimg(bfr+1:bfr+rr{2},cc{1}+cc{2}-2:cc{1}+cc{2}+2,:) = ones(rr{2},5,3);
    
    combimg(rnew+1:rnew+rr{4},cc{4}-2:cc{4}+2,:) = ones(rr{4},5,3);
    combimg(rnew+1:rnew+rr{5},cc{4}+cc{5}-2:cc{4}+cc{5}+2,:) = ones(rr{5},5,3);
    

    imwrite(combimg,savename)
    
end
    

cd(homedir)


    
    
    
    
    
 
    
    