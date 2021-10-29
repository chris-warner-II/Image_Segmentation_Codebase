function [imgName, Method, Params] = parse_Fname(filename)

% syntax: [imgName, Method, Params] = parse_Fname(filename)
%
%


x = strfind(filename,'-'); 
space = find(filename==' ');
imgName = filename(1:space(1));
Method = filename(space(1)+1:x-1);
Params = filename(x+1:end-4);
