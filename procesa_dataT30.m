%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 700s


clear all

mods=[700:703 705:732 734:750];

for nn=1:length(mods)
    n=mods(nn)

lat=ncread(['attm' num2str(n) '_18702017.nc'],'lat');
lon=ncread(['attm' num2str(n) '_18702017.nc'],'lon');

c=0;
for k=1870:2016
    for i=1:12
        c=c+1;
        sst(:,:,c)=squeeze(ncread(['attm' num2str(n) '_18702017.nc'],'sst',[1 1 c],[96 48 1],[1 1 1]));
        pr1=squeeze(ncread(['attm' num2str(n) '_18702017.nc'],'precls',[1 1 c],[96 48 1],[1 1 1]));
        pr2=squeeze(ncread(['attm' num2str(n) '_18702017.nc'],'precnv',[1 1 c],[96 48 1],[1 1 1]));
        aux=mean(mean(squeeze(sst(52:65,23:26,c))));
        auy=mean(mean(squeeze(sst(52:65,24:25,c))));
        n34(c)=(aux+auy)/2; 
        pp(:,:,c)=pr1+pr2;
        slp(:,:,c)=ncread(['attm' num2str(n) '_18702017.nc'],'mslp',[1 1 c],[96 48 1],[1 1 1]);
        h5(:,:,c)=ncread(['attm' num2str(n) '_18702017.nc'],'gh',[1 1 6 c],[96 48 1 1],[1 1 1 1]);    
    end    
end
clear aux auy
dat=slp(:,1:19,:);
for i=1:12
%    mdat(:,:,i)=mean(dat(:,:,1332+i:12:1692),3); % clim 1981-2010
    mdat(:,:,i)=mean(dat(:,:,i:12:end),3); % clim 1981-2010
end
c=0;
for k=1870:2016
    for i=1:12
        c=c+1;
        aux(:,:,c)=dat(:,:,c)-mdat(:,:,i);
    end
end

% Weighting of geophysical data in principal component analysis, Chung and
% Nigan 1999, J. Geos. Res.

for i=1:length(19)
    factor=sqrt(cosd(lat(i)));
    aux(:,i,:)=aux(:,i,:)*factor;
end
N=length(n34);
for i=1:N
    F(:,i)=reshape(squeeze(aux(:,:,i)),96*19,1);
end
M=length(F(:,1));            % espacio
N=length(F(1,:));             % tiempo
[L,A,E,error]=EOF(F,100); % F enter with (N,M)!!!!!!!
L1=L(1)/sum(L)
auy=reshape(E(:,1),96,19)';
if mean(auy(13,:)-auy(6,:))>0
    fac=1;
else
    fac=-1;
end
sam=fac*A(:,1)/std(A(:,1));

eval(['save tropical' num2str(n) ' n34 sam sst slp h5 pp lat lon L1'])

end



%% 600s


clear all

mods=[600:605 607:619 622:650];

for nn=1:length(mods)
    n=mods(nn)

lat=ncread(['attm' num2str(n) '_18702017.nc'],'lat');
lon=ncread(['attm' num2str(n) '_18702017.nc'],'lon');

c=0;
for k=1870:2016
    for i=1:12
        c=c+1;
        sst(:,:,c)=squeeze(ncread(['attm' num2str(n) '_18702017.nc'],'sst',[1 1 c],[96 48 1],[1 1 1]));
        pr1=squeeze(ncread(['attm' num2str(n) '_18702017.nc'],'precls',[1 1 c],[96 48 1],[1 1 1]));
        pr2=squeeze(ncread(['attm' num2str(n) '_18702017.nc'],'precnv',[1 1 c],[96 48 1],[1 1 1]));
        aux=mean(mean(squeeze(sst(52:65,23:26,c))));
        auy=mean(mean(squeeze(sst(52:65,24:25,c))));
        n34(c)=(aux+auy)/2; 
        pp(:,:,c)=pr1+pr2;
        slp(:,:,c)=ncread(['attm' num2str(n) '_18702017.nc'],'mslp',[1 1 c],[96 48 1],[1 1 1]);
        h5(:,:,c)=ncread(['attm' num2str(n) '_18702017.nc'],'gh',[1 1 6 c],[96 48 1 1],[1 1 1 1]);    
    end    
end
clear aux auy
dat=slp(:,1:19,:);
for i=1:12
%    mdat(:,:,i)=mean(dat(:,:,1332+i:12:1692),3); % clim 1981-2010
    mdat(:,:,i)=mean(dat(:,:,i:12:end),3); % clim 1981-2010
end
c=0;
for k=1870:2016
    for i=1:12
        c=c+1;
        aux(:,:,c)=dat(:,:,c)-mdat(:,:,i);
    end
end

% Weighting of geophysical data in principal component analysis, Chung and
% Nigan 1999, J. Geos. Res.

for i=1:length(19)
    factor=sqrt(cosd(lat(i)));
    aux(:,i,:)=aux(:,i,:)*factor;
end
N=length(n34);
for i=1:N
    F(:,i)=reshape(squeeze(aux(:,:,i)),96*19,1);
end
M=length(F(:,1));            % espacio
N=length(F(1,:));             % tiempo
[L,A,E,error]=EOF(F,100); % F enter with (N,M)!!!!!!!
L1=L(1)/sum(L)
auy=reshape(E(:,1),96,19)';
if mean(auy(13,:)-auy(6,:))>0
    fac=1;
else
    fac=-1;
end
sam=fac*A(:,1)/std(A(:,1));

eval(['save extratrop' num2str(n) ' n34 sam sst slp h5 pp lat lon L1'])

end

%% 500s

clear all

mods=[500:534 536:546 548:549];

for nn=1:length(mods)
    n=mods(nn)

lat=ncread(['attm' num2str(n) '_18702017.nc'],'lat');
lon=ncread(['attm' num2str(n) '_18702017.nc'],'lon');

c=0;
for k=1870:2016
    for i=1:12
        c=c+1;
        sst(:,:,c)=squeeze(ncread(['attm' num2str(n) '_18702017.nc'],'sst',[1 1 c],[96 48 1],[1 1 1]));
        pr1=squeeze(ncread(['attm' num2str(n) '_18702017.nc'],'precls',[1 1 c],[96 48 1],[1 1 1]));
        pr2=squeeze(ncread(['attm' num2str(n) '_18702017.nc'],'precnv',[1 1 c],[96 48 1],[1 1 1]));
        aux=mean(mean(squeeze(sst(52:65,23:26,c))));
        auy=mean(mean(squeeze(sst(52:65,24:25,c))));
        n34(c)=(aux+auy)/2; 
        pp(:,:,c)=pr1+pr2;
        slp(:,:,c)=ncread(['attm' num2str(n) '_18702017.nc'],'mslp',[1 1 c],[96 48 1],[1 1 1]);
        h5(:,:,c)=ncread(['attm' num2str(n) '_18702017.nc'],'gh',[1 1 6 c],[96 48 1 1],[1 1 1 1]);    
    end    
end
clear aux auy
dat=slp(:,1:19,:);
for i=1:12
%    mdat(:,:,i)=mean(dat(:,:,1332+i:12:1692),3); % clim 1981-2010
    mdat(:,:,i)=mean(dat(:,:,i:12:end),3); % clim 1981-2010
end
c=0;
for k=1870:2016
    for i=1:12
        c=c+1;
        aux(:,:,c)=dat(:,:,c)-mdat(:,:,i);
    end
end

% Weighting of geophysical data in principal component analysis, Chung and
% Nigan 1999, J. Geos. Res.

for i=1:length(19)
    factor=sqrt(cosd(lat(i)));
    aux(:,i,:)=aux(:,i,:)*factor;
end
N=length(n34);
for i=1:N
    F(:,i)=reshape(squeeze(aux(:,:,i)),96*19,1);
end
M=length(F(:,1));            % espacio
N=length(F(1,:));             % tiempo
[L,A,E,error]=EOF(F,100); % F enter with (N,M)!!!!!!!
L1=L(1)/sum(L)
auy=reshape(E(:,1),96,19)';
if mean(auy(13,:)-auy(6,:))>0
    fac=1;
else
    fac=-1;
end
sam=fac*A(:,1)/std(A(:,1));

eval(['save control' num2str(n) ' n34 sam sst slp h5 pp lat lon L1'])

end

