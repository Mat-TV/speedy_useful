%% Programa para leer los archivos de salida de SPEEDY y calcular variables derivadas
aseo % https://github.com/Trufumut/MATLAB_useful/blob/main/aseo.m
%% Preliminares
%ruta a los archivos
cd '/home/matt/Documentos/SPEEDY/Salidas'
% nombre='attm500_1870.nc'; % modificado con anomalía = 1 
nombre='attm600_1870.nc'; % con anom del HadISST
% ncdisp(nombre)
fechas=double(ncread(nombre,'time'));
fecha0 = datenum('01-01-1870');
fechas = fechas./24 + fecha0; %
%
lat=ncread(nombre,'lat'); 
lon=ncread(nombre,'lon'); 
% Selección de SST en la zona NIÑO 3.4
SSTa=double(ncread(nombre,'ssta',[52 23 1],[14 4 Inf])); % (longitudes,latitudes,tiempos(meses))
%% Cálculo de índice El Niño
for i=1:length(SSTa(1,1,:))-1
    ENOS(i) = mean(mean(SSTa(:,:,i)));
end
%
figure
anomaly(fechas(1:end-1),ENOS(1:end))
datetick('x','yyyy','keepticks')
axis tight
grid minor
title('Anomalía de TSM en la zona 5°S<->5°N y 170°O<->120°O')
%
%% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 