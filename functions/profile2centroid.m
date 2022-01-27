function peak = profile2centroid(ms)
%convert the profile MS to the centroid MS
%the input variable is a n row-by- 2 columns profile mode mass spectrum
  delta1=diff(ms(:,2));
  delta2=flipud(diff(flipud(ms(:,2))));
  a=sign(delta1);
  b=sign(delta2);
  indx=find([0;a].*[b;0]>0);
  ms_centroid=ms(indx,:);
  ms=[];a=[];b=[];delta1=[];delta2=[];indx=[];
  peak=ms_centroid(find(ms_centroid(:,2)>200),:);
  ms_centroid=[];
  clc
end

