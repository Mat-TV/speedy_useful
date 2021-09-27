%% Programa adaptado para calcular Cond. Borde para SPEEDY, de SST y SeaICE HadISST
% ESTA VEZ PARA RESOLUCION T63, por lo que cambia definiciones lon (X) y lat (Y)
%2018-04-24

%ruta a los archivos
path1 = '/home/danveloso/Documentos/Fondecyt_Garreaud/SPEEDY_SST-SIC/';

%----Climatology that comes with Speedy package

fid = fopen([path1,'t63/clim/sst_8190clim.t63.sea.grd'], 'r')
variable = fread(fid, 'float','b'); % archivo de dim 55320x1

% sstc=reshape(variable,4610,12); %4610*12=55320
% sstc=sstc(2:48*96+1,:); %48*96=4608 --> valor 1 y end es 2.5829e-41 (distinto a -9.999e+19 = NaN)
% sstc=reshape(sstc,96,48,12);
sstc=reshape(variable,18434,12); %18434*12=221208
sstc=sstc(2:96*192+1,:); %96*192=18432 --> valor 1 y end es 1.0331e-40 (distinto a -9.999e+19 = NaN)
sstc=reshape(sstc,192,96,12);

sstc=flipdim(sstc,2); % --> averiguar si parte desde +90° o desde -90°--> Parece que desde -90°
sstc=permute(sstc,[3 2 1 ]); %%%--------> IMPORTANTE! Dimensiones = (tiempo,latitud,longitud)

mm=min(min(min(sstc))); %-9.999e+19 = NaN (el valor + pequeño)
ii=find(sstc==mm); % ii contiene posiciones de los valores NAN,
SSTC = sstc; % para climatologia sin NaN
sstc(ii)=NaN;


%%% IMPORTANTE: Resolución y puntos de grilla usados en SPEEDY
X=0:1.875:360; X=X(1:192);
Y=[-88.572  -86.723  -84.862  -82.999  -81.135  -79.271  -77.406  -75.541  -73.676  -71.811  -69.946  -68.081  -66.216  -64.351  -62.486  -60.620  -58.755  -56.890  -55.025 ...
   -53.160  -51.294  -49.429  -47.564  -45.699  -43.833  -41.968  -40.103  -38.238  -36.372  -34.507  -32.642  -30.777  -28.911  -27.046  -25.181  -23.316  -21.450  -19.585 ...
   -17.720  -15.855  -13.989  -12.124  -10.259  -8.394   -6.528   -4.663   -2.798   -0.933    0.933    2.798    4.663    6.528    8.394    10.259   12.124   13.989   15.855 ...
    17.720   19.585   21.450   23.316   25.181   27.046   28.911   30.777   32.642   34.507   36.372   38.238   40.103   41.968   43.833   45.699   47.564   49.429   51.294  ...  
    53.160   55.025   56.890   58.755   60.620   62.486   64.351   66.216   68.081   69.946   71.811   73.676   75.541   77.406   79.271   81.135   82.999   84.862   86.723 ...   
    88.572];


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
% climaux = permute(clim1,[3 2 1]);
% figure
% contourf(Lon,lat,climaux(:,:,11)')
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
% climaux = permute(clim,[3 2 1]);
% figure
% contourf(X,Y,climaux(:,:,11)')
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
% contourf(X,Y,climaux(:,:,end)')
% colorbar

%------ Save New Clim 

clim2=flipdim(clim,2);   %flips Y dim
clim2=permute(clim2,[3 2 1]);
clim2=reshape(clim2,[18432,12]);
dummy=ones(18434,12)*variable(1);   %creates a matrix of 1.0331e-40
dummy(2:96*192+1,:)=clim2;          %puts good values within matrix
dummy=dummy(:);

fidw = fopen([path1,'Cond_borde/sst_clim6110Hadisst.t63.sea.grd'],'w')
count = fwrite(fidw,dummy,'float','b');

clear dummy

%------ Save New Anomalies

SSTa2=flipdim(SSTa,2);   %flips Y dim
% SSTa2=permute(SSTa2,[1 3 2]);
SSTa2=permute(SSTa2,[3 2 1]);
% SSTa2=SSTa(:,:)';
SSTa2=reshape(SSTa2,[18432,1768]);
dummy=ones(18434,1768)*variable(1);   %creates a matrix of 1.0331e-40
dummy(2:96*192+1,:)=SSTa2;            %puts good values within matrix
dummy=dummy(:);

fidw =fopen([path1,'Cond_borde/sst_anom6110Hadisst.t63.sea.grd'],'w')
count = fwrite(fidw,dummy,'float','b');

%% Climatologia para el hielo marino

%ruta a los archivos
path1 = '/home/danveloso/Escritorio/GEOFISICA/Fondecyt_Garreaud/SPEEDY_SST-SIC/';

fid = fopen([path1,'seaice_7908clim.t30.sea.grd'], 'r')
variable = fread(fid, 'float','b'); % archivo de dim 55320x1

icec=reshape(variable,4610,12); %4610*12=55320
icec=icec(2:48*96+1,:); %48*96=4608 --> valor 1 y end es 2.5829e-41 (di%%%% --------------- South Pacific Experiment ----------------
% limites: 30-90°S, 140°E-70°W
sstan1 = sstan;
%limites latitudinales y longitudinales
sstan1(:,32:38,:) = sstan1(:,32:38,:).*cosd(x);
sstan1(:,39:end,:) = 0;
%limites longidinales
sstan1(:,:,70:76) = sstan1(:,:,70:76).*sind(x);
sstan1(:,:,156:162) = ssta1(:,:,156:162).*cosd(x);
sstan1(:,:,1:69) = 0;
sstan1(:,:,163:end) = 0;

%X(76) = 140.625
%X(156) = 290.625 = 70.625W

% --- Guardamos nuevas anomalias ---
sstan1=flipdim(sstan1,2);   %flips Y dim
sstan1=permute(sstan1,[3 2 1]);
sstan1=reshape(sstan1,[18432,1768]);
dummy1=ones(18434,1768)*variable(1);   %creates a matrix of 1.0331e-40
dummy1(2:96*192+1,:) = sstan1;            %puts good values within matrix
dummy1=dummy1(:);

fidw1 =fopen([path1,'Cond_borde/sst_anom6110Hadisst_SPexp.t63.sea.grd'],'w')
count1 = fwrite(fidw1,dummy1,'float','b');stinto a -9.999e+19 = NaN)
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

%% -------------------------------------------------------------------------------------------------
%% Se editan las Condiciones de Borde, para realizar 3 nuevos experimentos en el Pacifico Sur
%% -------------------------------------------------------------------------------------------------
%2018-05-03

path1 = '/home/danveloso/Documentos/Fondecyt_Garreaud/SPEEDY_SST-SIC/';

fid = fopen([path1,'Cond_borde/sst_anom6110Hadisst.t63.sea.grd'],'r')

variable = fread(fid, 'float','b'); %archivo de dimension 32591312x1

sstan=reshape(variable,18434,1768); %18434*1768
sstan=sstan(2:96*192+1,:); %96*192=18432 --> valor 1 y end es  1.0331e-40 (distinto a -9.999e+19 = NaN)
sstan=reshape(sstan,192,96,1768);

sstan=flipdim(sstan,2);
sstan=permute(sstan,[3 2 1 ]); %%%--------> Dimensiones = (tiempo,latitud,longitud)

%%% IMPORTANTE: Resolución y puntos de grilla usados en SPEEDY
X=0:1.875:360; X=X(1:192);
Y=[-88.572  -86.723  -84.862  -82.999  -81.135  -79.271  -77.406  -75.541  -73.676  -71.811  -69.946  -68.081  -66.216  -64.351  -62.486  -60.620  -58.755  -56.890  -55.025 ...
   -53.160  -51.294  -49.429  -47.564  -45.699  -43.833  -41.968  -40.103  -38.238  -36.372  -34.507  -32.642  -30.777  -28.911  -27.046  -25.181  -23.316  -21.450  -19.585 ...
   -17.720  -15.855  -13.989  -12.124  -10.259  -8.394   -6.528   -4.663   -2.798   -0.933    0.933    2.798    4.663    6.528    8.394    10.259   12.124   13.989   15.855 ...
    17.720   19.585   21.450   23.316   25.181   27.046   28.911   30.777   32.642   34.507   36.372   38.238   40.103   41.968   43.833   45.699   47.564   49.429   51.294  ...  
    53.160   55.025   56.890   58.755   60.620   62.486   64.351   66.216   68.081   69.946   71.811   73.676   75.541   77.406   79.271   81.135   82.999   84.862   86.723 ...   
    88.572];

 %Angulos para aplicar pesos como funcion sinuosoidal entre 0 y 90°
x = linspace(0,90,7); 


%%%% --------------- South Pacific Experiment ----------------
% limites: 30-90°S, 140°E-70°W
sstan1 = sstan;
%limites latitudinales y longitudinales
sstan1(:,32:38,:) = sstan1(:,32:38,:).*cosd(x);
sstan1(:,39:end,:) = 0;
%limites longidinales
sstan1(:,:,70:76) = sstan1(:,:,70:76).*reshape(sind(x),[1 1 7]);
sstan1(:,:,156:162) = sstan1(:,:,156:162).*reshape(cosd(x),[1 1 7]);
sstan1(:,:,1:69) = 0;
sstan1(:,:,163:end) = 0;

%X(76) = 140.625
%X(156) = 290.625 = 70.625W
% Y(32) = -30.7770
% Y(38) = -19.5850

%%% Verificacion de que caja o region de aomalias este bien
% sstan0 = permute(sstan,[3 2 1]);
% figure, subplot(211),pcolor(X,Y(1:end/2),sstan0(:,1:end/2,end-1)'),colorbar, caxis([-2 2]), shading flat
% subplot(212),pcolor(X,Y(1:end/2),sstan3(:,1:end/2,end-1)'),colorbar, caxis([-2 2]), shading flat


% --- Guardamos nuevas anomalias ---
sstan1=flipdim(sstan1,2);   %flips Y dim
sstan1=permute(sstan1,[3 2 1]);
sstan1=reshape(sstan1,[18432,1768]);
dummy1=ones(18434,1768)*variable(1);   %creates a matrix of 1.0331e-40
dummy1(2:96*192+1,:) = sstan1;            %puts good values within matrix
dummy1=dummy1(:);

fidw1 =fopen([path1,'Cond_borde/sst_anom6110Hadisst_SPexp.t63.sea.grd'],'w')
count1 = fwrite(fidw1,dummy1,'float','b');



%%%% --------------- Southern Ocean Experiment ----------------

%limites: 50°S-90°S, 140°E-70°W
sstan2 = sstan;
%limites latitudinales y longitudinales
sstan2(:,21:27,:) = sstan2(:,21:27,:).*cosd(x);
sstan2(:,28:end,:) = 0;
%limites longidinales
sstan2(:,:,70:76) = sstan2(:,:,70:76).*reshape(sind(x),[1 1 7]);
sstan2(:,:,156:162) = sstan2(:,:,156:162).*reshape(cosd(x),[1 1 7]);
sstan2(:,:,1:69) = 0;
sstan2(:,:,163:end) = 0;

% Y(21) = -51.2940
% Y(27) = -40.1030

%%% Verificacion de que caja o region de aomalias este bien
% sstan0 = permute(sstan,[3 2 1]);
% sstan22 = permute(sstan2,[3 2 1]);
% figure, subplot(211),pcolor(X,Y(1:end/2),sstan22(:,1:end/2,end-1)'),colorbar, caxis([-2 2]), shading flat
% subplot(212),pcolor(X,Y(1:end/2),sstan0(:,1:end/2,end-1)'),colorbar, caxis([-2 2]), shading flat
% clear sstan22 sstan0


% --- Guardamos nuevas anomalias ---
sstan2=flipdim(sstan2,2);   %flips Y dim
sstan2=permute(sstan2,[3 2 1]);
sstan2=reshape(sstan2,[18432,1768]);
dummy2=ones(18434,1768)*variable(1);   %creates a matrix of 1.0331e-40
dummy2(2:96*192+1,:)=sstan2;            %puts good values within matrix
dummy2=dummy2(:);

fidw2 =fopen([path1,'Cond_borde/sst_anom6110Hadisst_SOexp.t63.sea.grd'],'w')
count2 = fwrite(fidw2,dummy2,'float','b');



%%%% --------------- Midlatitude South Pacific Experiment ----------------

%limites: 30°S-50°S, 140°E-70°W
sstan3 = sstan;
%limites latitudinales y longitudinales
sstan3(:,32:38,:) = sstan3(:,32:38,:).*cosd(x);
sstan3(:,16:22,:) = sstan3(:,16:22,:).*sind(x);
sstan3(:,1:15,:) = 0;
sstan3(:,39:end,:) = 0;
%limites longidinales
sstan3(:,:,70:76) = sstan3(:,:,70:76).*reshape(sind(x),[1 1 7]);
sstan3(:,:,156:162) = sstan3(:,:,156:162).*reshape(cosd(x),[1 1 7]);
sstan3(:,:,1:69) = 0;
sstan3(:,:,163:end) = 0;

% Y(16) = -60.6200
% Y(22) = -49.4290

%%% Verificacion de que caja o region de aomalias este bien
sstan0 = permute(sstan,[3 2 1]);
sstan33 = permute(sstan3,[3 2 1]);
for i=1:1768
subplot(211),pcolor(X,Y(1:end/2),sstan33(:,1:end/2,i)'),colorbar, caxis([-2 2]), shading flat
subplot(212),pcolor(X,Y(1:end/2),sstan0(:,1:end/2,i)'),colorbar, caxis([-2 2]), shading flat
title(num2str(i/12+1870))
pause(0.05)
end
clear sstan33 sstan0


% --- Guardamos nuevas anomalias ---
sstan3=flipdim(sstan3,2);   %flips Y dim
sstan3=permute(sstan3,[3 2 1]);
sstan3=reshape(sstan3,[18432,1768]);
dummy3=ones(18434,1768)*variable(1);   %creates a matrix of 1.0331e-40
dummy3(2:96*192+1,:)=sstan3;            %puts good values within matrix
dummy3=dummy3(:);

fidw3 =fopen([path1,'Cond_borde/sst_anom6110Hadisst_MLexp.t63.sea.grd'],'w')
count3 = fwrite(fidw3,dummy3,'float','b');
