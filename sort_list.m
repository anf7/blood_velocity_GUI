function outs = sort_list(ins)

if ispc
    divider = '\';
else
    divider = '/';
end

chararray = cell(0);
outs = struct([]);
if length(ins) > 1
    for n = 1:length(ins)
        chararray{n} = ins(n).path{1};
    end
    [tempstr, i] = sort(chararray);
    for n = 1:length(tempstr)
        if ~isequal(tempstr{n}(end),divider)
            tempstr{n}(end+1) = divider;
        end
        outs(n).path = tempstr(n);
    end 
    for n = 1:length(outs)
        outs(n).tiffs = ins(i(n)).tiffs;
        outs(n).level = ins(i(n)).level;
    end
else
    outs = ins;
end
