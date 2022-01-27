function [image] = massimage(peaklist,mz_target,mass_tolerance)
% construct certain specified ion image out of the peaklist
% batch searching the intensity values of the target mz within given 
% mass_tolerance in each scan pixels
% the mass_tolerance was defined as the relative error
for i=1:length(peaklist)
     for j=1:length(peaklist{i,1})
        distance=abs(mz_target-peaklist{i,1}{j,1}(:,1));
        if min(distance)<=mass_tolerance
          [r,~]=find(distance<=mass_tolerance);
          image(i,j)=sum(peaklist{i,1}{j,1}(r,2));
        else
          image(i,j)=0;
        end
     end
end
clear i j r distance mz_target
clear mass_tolerance
end