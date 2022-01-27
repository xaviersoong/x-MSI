function cmz = msAverage(mz_matrix)
%generate an average mass spectrum from a set of mass spectra
%   Detailed explanation goes here
cache=[];
for i=1:size(mz_matrix,2)
   cache=[cache;mz_matrix(:,i)];
end
clear i delta
cache=cache(find(cache>0));
cache=sort(cache);
delta=[0;diff(cache)];
indx=find(delta>0.01);
for i=1:length(indx)+1
    if i==1
       cmz_mean(i,1)=mean(cache(1:indx(1)-1));
    elseif i==length(indx)+1
       cmz_mean(i,1)=mean(cache(indx(i-1):end));
    else
       cmz_mean(i,1)=mean(cache(indx(i-1):indx(i)-1));
    end
end
cmz=cmz_mean;
clear i indx
clc
end

