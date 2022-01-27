function [ peaklist ] = batchmzxmlread( prefix,rows )
% import a set of mzXML format MS data files
%   Detailed explanation goes here
peaklist={};
for i=1:9
    file_name=[prefix,'0',num2str(i),'.mzXML']
    temp=mzxmlread(file_name);
    [peaklist{i,1},~]=mzxml2peaks(temp);
end    
for i=10:rows
    temp=mzxmlread([prefix,num2str(i),'.mzXML']);
    [peaklist{i,1},~]=mzxml2peaks(temp);
end
clear temp i
end

