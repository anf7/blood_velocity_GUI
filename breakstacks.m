function breakstacks

breaksize = 40;

[filename,pathname] = uigetfile('*.tif','Select .tif files to break...','MultiSelect','on');

if ~iscell(filename)
    tempfname = filename;
    filename = cell(1);
    filename{1} = tempfname;
end


for i = 1:length(filename)
    fullfilename = strcat(pathname,filename{i});
    
    
    info = imfinfo(fullfilename);
    r = info(1).Height;
    c = info(1).Width;
    z = length(info);

    for j = 1:2
        
    
        imstack = zeros(breaksize,r,c);
        if j == 1
            for n = 1:breaksize
                imstack(n,:,:) = imread(fullfilename,n);
            end
            newname = strcat(fullfilename(1:end-4),num2str(1),'to',num2str(breaksize),'.tif');
        else
            for n = 1:breaksize
                imstack(n,:,:) = imread(fullfilename,z-breaksize+n);
            end
            newname = strcat(fullfilename(1:end-4),num2str(z-breaksize+1),'to',num2str(z),'.tif');
        end
            
            
        imstack = double(imstack)/double(max(imstack(:)));
        
        for n = 1:breaksize
            if j == 1 && n == 1
                continue
            end
            frame = squeeze(imstack(n,:,:));
            imwrite(frame, newname, 'WriteMode', 'append', 'Compression', 'none')
        end
    end
end

