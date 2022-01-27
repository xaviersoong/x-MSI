function [datacube] = batchmassimage(peaklist,cmz,ppm_value)
% construct a set of ion images from peaklist data
% cmz: a vector consisted of ions exact m/z values
% ppm_value: mass tolerance (default 0.005)
iter=0;
 for m=1:length(cmz)
   target=cmz(m,1);
   datacube{m,1}=target;
   datacube{m,2}=massimage(peaklist,target,ppm_value);
   target=[];
   iter=iter+1
 end
 clear target m iter
 clc
end

