function [image] = indximg(vector,row,column )
% construct an MS image
%   vector contain all elements that need to be filled according to their 
% row and column index
for i=1:length(vector)
    image(row(i,1),column(i,1))=vector(i,1);
end
clear i
end

