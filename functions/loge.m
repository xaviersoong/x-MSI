function [ y ] = loge( x )
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here
y=log(x);
for i=1:size(y,1)
    for j=1:size(y,2)
        if y(i,j)==inf|y(i,j)==-inf
            y(i,j)=0;
        end
    end
end
            

end

