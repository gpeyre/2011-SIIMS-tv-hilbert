function [f1,f2,maxi] = instantaneous_frequency(F)
  [maxi,freq] = max(F(:));
  q = size(F,1);
  [f1,f2] = ind2sub(q,freq);
  f1 = q/2+1-f1;
  f2 = q/2+1-f2;