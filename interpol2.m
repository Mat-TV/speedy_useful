%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Example: To interpolate NCEP w500 mb to T42
%                nc=ncload('ORO_pogaml.nc','lon','lat');
%                nc=ncload('w500mb_ncep_1949-2002.cdf');
%                [omega]=interpol2(X,Y,squeeze(vvel),lon,lat)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [reg]=interpol2(xaxis1,yaxis1,X,xaxis0,yaxis0)

[lt,ly,lx]=size(X);

for i=1:lt 
  reg(i,:,:)=interp2(xaxis1',yaxis1,squeeze(X(i,:,:)),xaxis0',yaxis0);
end
