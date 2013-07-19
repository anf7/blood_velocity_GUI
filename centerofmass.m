function [rpos cpos] = centerofmass(input)

[r c] = size(input);
rlinsp = linspace(1,r,r)';
clinsp = linspace(1,c,c);

sumacrossrows = sum(input,2);
sumdowncolumns = sum(input,1);

rpos = sum(sumacrossrows.*rlinsp)/sum(sumacrossrows);
cpos = sum(sumdowncolumns.*clinsp)/sum(sumdowncolumns);