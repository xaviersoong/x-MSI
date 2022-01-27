function [peaklist] = batchcdfread(prefix,rows)
%UNTITLED import a set of cdf format MS data files
%   Detailed explanation goes here
peaklist={};
for i=1:9
    file_name=[prefix,'0',num2str(i),'.cdf']
    temp=mzcdfread(file_name);
    [peaklist{i,1},~]=mzcdf2peaks(temp);
end    
for i=10:rows
    temp=mzcdfread([prefix,num2str(i),'.cdf']);
    [peaklist{i,1},~]=mzcdf2peaks(temp);
end
clear temp i
end

