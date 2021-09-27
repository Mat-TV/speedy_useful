%% Programa adaptado para calcular Cond. Borde para SPEEDY, de SST y SeaICE HadISST
% ESTA VEZ PARA RESOLUCION T47, por lo que cambia definiciones lon (X) y lat (Y)
% 2018-04-24

%ruta a los archivos
path1 = '/home/danveloso/Documentos/Fondecyt_Garreaud/SPEEDY_SST-SIC/';

%----Climatology that comes with Speedy package

fid = fopen([path1,'t47/clim/sst_7908clim.t47.sea.grd'], 'r')
variable = fread(fid, 'float','b'); % archivo de dim 55320x1

% sstc=reshape(variable,4610,12); %4610*12=55320
% sstc=sstc(2:48*96+1,:); %48*96=4608 --> valor 1 y end es 2.5829e-41 (distinto a -9.999e+19 = NaN)
% sstc=reshape(sstc,96,48,12);
sstc=reshape(variable,10370,12); %10370*12=124440
sstc=sstc(2:72*144+1,:); %48*96=4608 --> valor 1 y end es 5.8115e-41 (distinto a -9.999e+19 = NaN)
sstc=reshape(sstc,144,72,12);

sstc=flipdim(sstc,2); % --> averiguar si parte desde +90° o desde -90°--> Parece que desde -90°
sstc=permute(sstc,[3 2 1 ]); %%%--------> IMPORTANTE! Dimensiones = (tiempo,latitud,longitud)

mm=min(min(min(sstc))); %-9.999e+19 = NaN (el valor + pequeño)
ii=find(sstc==mm); % ii contiene posiciones de los valores NAN,
SSTC = sstc; % para climatologia sin NaN
sstc(ii)=NaN;


%%% IMPORTANTE: Resolución y puntos de grilla usados en SPEEDY
X=0:2.5:360; X=X(1:144);
Y=[-88.100 -85.638  -83.161  -80.681  -78.200  -75.719  -73.237  -70.755  -68.272  -65.790  -63.308  -60.825  -58.343  -55.860  -53.377  -50.895  -48.412  -45.930  -43.447 ...
-40.964  -38.482  -35.999  -33.516  -31.034  -28.551  -26.068  -23.586  -21.103  -18.620  -16.138  -13.655  -11.172  -8.689   -6.207   -3.724   -1.241    1.241    3.724 ...
   6.207    8.689    11.172   13.655   16.138   18.620   21.103   23.586   26.068   28.551   31.034   33.516   35.999   38.482   40.964   43.447   45.930   48.412   50.895 ...
   53.377   55.860   58.343   60.825   63.308   65.790   68.272   70.755   73.237   75.719   78.200   80.681   83.161   85.638   88.100];


%%%% CREATE NEW CLIMATOLOGY AND ANOMALIES USING HadISST  %%%%%%%%%%%

%----ERSSTv.2 SST Jan1854-Jan2007 --> largo 1837
%----HadISST SST 1870/01-2017/04, dim=360x180x1768

SST=ncread([path1,'HadISST_sst.nc'],'sst');
%SST2=ncread('/home/marcelo/data/ERSSTv.2_Jan2007.cdf','SST');
lon=ncread([path1,'HadISST_sst.nc'],'longitude'); % -180:+180
lat=ncread([path1,'HadISST_sst.nc'],'latitude'); % +89:-89

%cambiamos dimension de longitud de -180:180 a 0:360
lon = lon([181:end 1:180]); lon(181:end) = lon(181:end)+360;
SST = SST([181:end 1:180],:,:);

SST=permute(SST,[3 2 1]); %---> (tiempo,lat,lon) --> así se necesita?
%SST2=reshape(SST2,1,180,89); SST2=permute(SST2,[1 3 2]);
%SST=[SST; SST2];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Se calcula primero la climatologia. Para ello se eliminarán los -1000°C
id1000 = find(SST == -1000);
SST(id1000) = NaN;

%---Se desea obtener climatologia 1961-2010 (50 años)
tiempo = datenum(1870,1:1768,1);
% 01/1961 = posicion 1093 ; 12/2010 = posicion 1692
SSTclim = SST(1093:1692,:,:);
SSTclim = reshape(SSTclim,[12 50 180 360]);
clim1 = squeeze(nanmean(SSTclim,2));

%Fills land with values interpolated from SST.
for j=1:12 %--> largo vector fecha
    clim1(j,:,:)=inpaint_nans(squeeze(clim1(j,:,:)));
    j
end
% Para corroborar que vamos bien
% clim2 = permute(clim1,[3 2 1]);
% figure
% contourf(Lon,lat,clim2(:,:,11)')
% colorbar

% Rellenamos la longitud 0 con valores para evitar NaN al interpolar
clim1(:,:,2:361) = clim1(:,:,1:360);
clim1(:,:,1) = (clim1(:,:,2) + clim1(:,:,end))./2;

%---Interpolation
Lon = single([0 0.5:1:359.5])';
[clim]=interpol2(Lon,flipdim(lat,1),flipdim(clim1,2),X,Y); %--> hacer flipdim SST y lat

%Set to -9.9990e+19 climatology over continents
clim(ii) = mm;

% Para ver mapa de climatologia
% clim2 = permute(clim,[3 2 1]);
% figure
% contourf(X,Y,clim2(:,:,11)')
% colorbar

% %Set values for sea ice as in the climatology that comes with Speedy.
% jj=find(sstc-273.15<0);
% clim(jj)=sstc(jj)-273.15;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Luego se calculan las anomalias
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Fills land with values interpolated from SST.
%Otherwise in the interpolation procedure the ERSST doesn't have
%good values along boundaries with continents.
for j=1:1768 %--> largo vector fecha
    SST2(j,:,:)=inpaint_nans(squeeze(SST(j,:,:)));
    j
end
% Rellenamos la longitud 0 con valores para evitar NaN al interpolar
SST2(:,:,2:361) = SST2(:,:,1:360);
SST2(:,:,1) = (SST2(:,:,2) + SST2(:,:,end))./2;

% % Rellenamos la longitud 0 con valores para evitar NaN al interpolar
% SST(:,:,2:361) = SST(:,:,1:360);
% SST(:,:,1) = (SST(:,:,2) + SST(:,:,end))./2;

%---Interpolation
Lon = single([0 0.5:1:359.5])';
[SST2]=interpol2(Lon,flipdim(lat,1),flipdim(SST2,2),X,Y); %--> hacer flipdim SST y lat

%Compute HadISST anomalies based on the climatology 1961-2010
for j=1:12
  m=j;
  while (m<=1768)
    SSTa(m,:,:)=SST2(m,:,:)-clim(j,:,:);
    m=m+12;
  end
end

%Set to 0 anomalies over continents
for t=1:1768
    ll = find(SSTC(1,:,:) == mm);
    SSTa(t,ll)=0;
end

clim=clim+273.15;

% break

% Para ver mapa final
% climaux = permute(SSTa,[3 2 1]);
% figure
% contourf(X,Y,climaux(:,:,500)')
% colorbar


%------ Save New Clim 

clim2=flipdim(clim,2);   %flips Y dim
clim2=permute(clim2,[3 2 1]);
clim2=reshape(clim2,[10368,12]);
dummy=ones(10370,12)*variable(1);   %creates a matrix of 5.8115e-41
dummy(2:72*144+1,:)=clim2;          %puts good values within matrix
dummy=dummy(:);

fidw = fopen([path1,'Cond_borde/sst_clim6110Hadisst.t47.sea.grd'],'w')
count = fwrite(fidw,dummy,'float','b');

clear dummy

%------ Save New Anomalies

SSTa2=flipdim(SSTa,2);   %flips Y dim
% SSTa2=permute(SSTa2,[1 3 2]);
SSTa2=permute(SSTa2,[3 2 1]);
% SSTa2=SSTa(:,:)';
SSTa2=reshape(SSTa2,[10368,1768]);
dummy=ones(10370,1768)*variable(1);   %creates a matrix of 5.8115e-41
dummy(2:72*144+1,:)=SSTa2;            %puts good values within matrix
dummy=dummy(:);

fidw =fopen([path1,'Cond_borde/sst_anom6110Hadisst.t47.sea.grd'],'w')
count = fwrite(fidw,dummy,'float','b');

%% Climatologia para el hielo marino

%ruta a los archivos
path1 = '/home/danveloso/Escritorio/GEOFISICA/Fondecyt_Garreaud/SPEEDY_SST-SIC/';

fid = fopen([path1,'seaice_7908clim.t30.sea.grd'], 'r')
variable = fread(fid, 'float','b'); % archivo de dim 55320x1

icec=reshape(variable,4610,12); %4610*12=55320
icec=icec(2:48*96+1,:); %48*96=4608 --> valor 1 y end es 2.5829e-41 (distinto a -9.999e+19 = NaN)
icec=reshape(icec,96,48,12);

icec=flipdim(icec,2); % --> averiguar si parte desde +90° o desde -90°--> Parece que desde -90°
icec=permute(icec,[3 2 1 ]); %%%--------> IMPORTANTE! Dimensiones = (tiempo,latitud,longitud)

% mm=min(min(min(icec))); %-9.999e+19 = NaN (el valor + pequeño)
% ii=find(icec==mm); % ii contiene posiciones de los valores NAN,
% ICEC = icec; % para climatologia sin NaN
% icec(ii)=NaN;

%%% IMPORTANTE: Resolución y puntos de grilla usados en SPEEDY
X=0:3.75:360; X=X(1:96);
Y=[-87.159   -83.479   -79.777   -76.070   -72.362   -68.652 -64.942   -61.232  -57.521   -53.810   -50.099   -46.389   -42.678   -38.967 -35.256   -31.545   -27.833  -24.122   -20.411   -16.700   -12.989    -9.278 -5.567    -1.856     1.856     5.567    9.278    12.989    16.700    20.411  24.122    27.833    31.545    35.256    38.967  42.678    46.389    50.099 53.810    57.521    61.232    64.942    68.652    72.362 76.070    79.777 83.479    87.159];

%Hielo de hadISST
ICE=ncread([path1,'HadISST_ice.nc'],'sic');
lon=ncread([path1,'HadISST_ice.nc'],'longitude'); % -180:+180
lat=ncread([path1,'HadISST_ice.nc'],'latitude'); % +89:-89

%cambiamos dimension de longitud de -180:180 a 0:360
lon = lon([181:end 1:180]); lon(181:end) = lon(181:end)+360;
ICE = ICE([181:end 1:180],:,:);

ICE=permute(ICE,[3 2 1]); %---> (tiempo,lat,lon)

%---Se desea obtener climatologia 1961-2010 (50 años)
tiempo = datenum(1870,1:1768,1);
% 01/1961 = posicion 1093 ; 12/2010 = posicion 1692
ICEclim = ICE(1093:1692,:,:);
ICEclim = reshape(ICEclim,[12 50 180 360]);
clim1 = squeeze(nanmean(ICEclim,2));

%Fills land with zeros.
clim1(isnan(clim1)) = 0;
% for j=1:12 %--> largo vector fecha
%     clim1(j,:,:)=inpaint_nans(squeeze(clim1(j,:,:)));
%     j
% end

% Rellenamos la longitud 0 con valores para evitar NaN al interpolar
clim1(:,:,2:361) = clim1(:,:,1:360);
clim1(:,:,1) = (clim1(:,:,2) + clim1(:,:,end))./2;

%---Interpolation
Lon = single([0 0.5:1:359.5])';
[clim]=interpol2(Lon,flipdim(lat,1),flipdim(clim1,2),X,Y); %--> hacer flipdim SST y lat

break

%------ Save New Clim 

clim2=flipdim(clim,2);   %flips Y dim
clim2=permute(clim2,[3 2 1]);
clim2=reshape(clim2,[4608,12]);
dummy=ones(4610,12)*variable(1);   %creates a matrix of 2.5829e-41
dummy(2:48*96+1,:)=clim2;          %puts good values within matrix
dummy=dummy(:);

fidw = fopen([path1,'Cond_borde/sic_clim6110Hadisst.t30.sea.grd'],'w')
count = fwrite(fidw,dummy,'float','b');

clear dummy

