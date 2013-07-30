function filetree = find_subdirectories(startdir) 
homedir = pwd;

%% Change directory naming conventions if Mac OS.
if ismac
    divider = '/';
else
    divider = '\';
end

%% User selects top level directory
topdir = uigetdir(startdir,...
    'Select the highest level data directory');

currentdir = topdir;

%% Find number of directory divider symbols in order to identify level in filetree
topdivider = length(strfind(topdir,divider));

%% Define first directory on list as top level directory
dirlist(1).path = currentdir;
dirlist(1).dir = [];
dirlist(1).level = 0;


%% Find subdirectories of top level directory
outercount = 0;
innercount = 1;

while outercount ~= innercount
 
outercount = outercount + 1;

    if outercount == 1
        cd(currentdir)
    else
        cd(strcat(dirlist(outercount).path, divider, dirlist(outercount).dir));
        currentdir = pwd;
    end

    td = dir('*.tif');
    counter = 1;
    tiffs = [];
    for n = 1:length(td)
        if length(imfinfo(td(n).name)) > 1
             tiffs{counter} = td(n).name;
             counter = counter + 1;
        end
    end
    dirlist(outercount).tiffs = tiffs; 

    
    d = dir;
    for n = 1:size(d,1)
        if d(n).isdir == 1 && ~(strcmp(d(n).name,'.') || strcmp(d(n).name,'..'))
            innercount = innercount + 1;
            directory = d(n).name;
            level = length(strfind(currentdir,divider)) - topdivider + ~isempty(directory);
            dirlist(innercount).path = currentdir;
            dirlist(innercount).dir = directory;
            dirlist(innercount).level = level;
            dirlist(innercount).preserve = 0;
        end
    end
end

cd(homedir)

for n = 1:length(dirlist)
    if ~isempty(dirlist(n).tiffs) 
        dirlist(n).preserve = 1;
        for p = 1:length(dirlist)
            if strfind(dirlist(n).path,strcat(dirlist(p).path, divider, dirlist(p).dir))
                dirlist(p).preserve = 1;
            end
        end
    end
end

counter = 1;
for n = 1:length(dirlist)
    if dirlist(n).preserve == 1
        filetree_p(counter).path = strcat(dirlist(n).path, divider, dirlist(n).dir);
        filetree_p(counter).tiffs = dirlist(n).tiffs;
        filetree_p(counter).level = dirlist(n).level;
        counter = counter + 1;
    end
end


[tempstr, i] = sort({filetree_p.path});
for n = 1:length(tempstr)
    if ~isequal(tempstr{n}(end),divider)
        tempstr{n}(end+1) = divider;
    end
    filetree(n).path = tempstr(n);
end 
for n = 1:length(filetree)
    filetree(n).tiffs = filetree_p(i(n)).tiffs;
    filetree(n).level = filetree_p(i(n)).level;
end


% if length(filetree_p) > 1   
%     [tempstr, i] = sort({filetree_p.path});
%     for n = 1:length(tempstr)
%         if ~isequal(tempstr{n}(end),divider)
%             tempstr{n}(end+1) = divider;
%         end
%         filetree(n).path = tempstr(n);
%     end 
%     for n = 1:length(filetree)
%         filetree(n).tiffs = filetree_p(i(n)).tiffs;
%         filetree(n).level = filetree_p(i(n)).level;
%     end
% else
%     for n = 1:length(filetree_p)
%         filetree(n).path = filetree_p(n).path;
%         filetree(n).tiffs = filetree_p(n).tiffs;
%         filetree(n).level = filetree_p(n).level;
%     end
% end







