function [filestr, indexcell] = filetree_string(instruct)


if ispc
    divsymb = '\';
else
    divsymb = '/';
end

filestr = cell(0);
indexcell = cell(0);
count = 1;
for n = 1:length(instruct)
    level = instruct(n).level;
    pathstr = char(instruct(n).path);
    if isequal(pathstr(end),divsymb)
        ps = pathstr;
    else
        ps = strcat(pathstr,divsymb);
    end
    if level == 0
        filestr{count} = ['<html><font color="white">',...
            pathstr(1:end-1),'</font></html>'];
    else
        divsymbind = strfind(pathstr,divsymb);
        a = divsymbind(end-1);
        filestr{count} = ['<html><font color="white">',repmat('.....',1,level),...
            '[',pathstr(a+1:end-1),']','</font></html>'];
    end
    indexcell{count} = ps;
    count = count + 1;
    for p = 1:length(instruct(n).tiffs)
        for q = 1:length(instruct(n).tiffs)
            if (q < p) && isequal(instruct(n).tiffs(p),instruct(n).tiffs(q))
                break
            elseif (q == p) && isequal(instruct(n).tiffs(p),instruct(n).tiffs(q))
                filestr{count} = [repmat('.....',1,level+1),char(instruct(n).tiffs(p))];
                indexcell{count} = strcat(ps,char(instruct(n).tiffs(p)));
                count = count + 1;
                break
            end
        end
    end
end
    



