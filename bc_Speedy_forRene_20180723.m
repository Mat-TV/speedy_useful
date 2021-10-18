%% Programa adaptado para calcular Cond. Borde para SPEEDY, de SST HadISST
% -------------------------------------------------------------------------------
%20180723 --> se obtienen condiciones de borde para experimentos con anomalias en 
% Subtropical SouthWestern Pacific (SSWP) de 0.5°C, 1.0°C, 1.5°C y 2.0°C (y 2.5°C ? )
% durante periodo de la Megasequia (2010-2014)
% -------------------------------------------------------------------------------
clear all

%% -------------------------------------------------------------------------------------------------
%% Se editan las Condiciones de Borde, para realizar 5 nuevos experimentos (Calenton en SSWP)
%% -------------------------------------------------------------------------------------------------

% path1 = '/home/danveloso/Documentos/Cond_borde_SPEEDY_revision/';
path1 = '/home/matt/speedy_ver41.5/data/bc/t30/';
path2 = '/home/matt/Documentos/SPEEDY/';
cd '/home/matt/Documentos/SPEEDY/' 


fid = fopen([path1,'clim/sst_clim6110Hadisst.t30.sea.grd'],'r')

variable = fread(fid, 'float','b'); % archivo de dim 55320x1
sstc=reshape(variable,4610,12); %4610*12=55320
sstc=sstc(2:48*96+1,:); %48*96=4608 --> valor 1 y end es 2.5829e-41 (distinto a -9.999e+19 = NaN)
sstc=reshape(sstc,96,48,12);

sstc=flipdim(sstc,2); % --> averiguar si parte desde +90° o desde -90°--> Parece que desde -90°
sstc=permute(sstc,[3 2 1 ]); %%%--------> IMPORTANTE! Dimensiones = (tiempo,latitud,longitud)

% mm=min(min(min(sstc))); %-9.999e+19 = NaN (el valor + pequeño)
% ii=find(sstc==mm); % ii contiene posiciones de los valores NAN,
% SSTC = sstc; % para climatologia sin NaN
% sstc(ii)=NaN;


%%% IMPORTANTE: Resolución y puntos de grilla usados en SPEEDY
X=0:3.75:360; X=X(1:96);
Y=[-87.159   -83.479   -79.777   -76.070   -72.362   -68.652 -64.942   -61.232  -57.521   -53.810   -50.099   -46.389   -42.678   -38.967 -35.256   -31.545   -27.833  -24.122   -20.411   -16.700   -12.989    -9.278 -5.567    -1.856     1.856     5.567    9.278    12.989    16.700    20.411  24.122    27.833    31.545    35.256    38.967  42.678    46.389    50.099 53.810    57.521    61.232    64.942    68.652    72.362 76.070    79.777 83.479    87.159];
% x = linspace(0,90,4); %Angulos para aplicar pesos como funcion sinuosoidal entre 0 y 90° 
x = linspace(0,1,4); % Funcion lineal simple


% ---------------------------------------------------------------
%%%% --------------- Experiment Clim + X °C SSWP   ----------------
% ---------------------------------------------------------------

% NOTA: Estos experimentos se han amortiguado en los bordes utilizando una
% funcion decreciente lineal, a diferencia de todos los anteriores donde se
% utilizo una funcion decreciente sinuosoidal. Notar que las esquinas de la
% caja, al igual que en los experimentos anteriores decrecen en magnitud a una tasa mayor

% limites: 25S-45S / 190E-210E

calenton = [0.5 1.0 1.5 2.0 2.5];

for k=1:length(calenton)
    %sumamos la anomalia al campo total
    sstc1{k} = zeros(size(sstc)) + calenton(k) ;

    %limites latitudinales y longitudinales
    % Y(12) = -46.3890; Y(18) = -24.1220
    % X(52) = 191.2500; X(57) = 210.00
    % sstc1(:,12:18,52:57) = sstc(:,12:18,52:57) + 0.5;

    %limites longidinales: suavizamos bordes
    sstc1{k}(:,:,49:52) = sstc1{k}(:,:,49:52).*reshape(x,[1 1 4]);
    sstc1{k}(:,:,57:60) = sstc1{k}(:,:,57:60).*reshape(flip(x),[1 1 4]);
    sstc1{k}(:,:,[1:48 61:end]) = 0;
    %limites latitudinales: suavizamos bordes
    sstc1{k}(:,9:12,:) = sstc1{k}(:,9:12,:).*reshape(x,[1 4 1]);
    sstc1{k}(:,18:21,:) = sstc1{k}(:,18:21,:).*reshape(flip(x),[1 4 1]);
    sstc1{k}(:,[1:8 22:end],:) = 0;

    %agregamos la climatologia al campo, obteniendo clim+0.5(SSWP)
    sstc1{k} = sstc1{k}+sstc;

    %%% Verificacion de que caja o region de anomalias este bien
    sstc0 = permute(sstc,[3 2 1]);
    sstc11{k} = permute(sstc1{k},[3 2 1]);
    figure(k)
    for i=1:12
        subplot(211),pcolor(X(end/3:5*end/6),Y(1:end/2+5),sstc11{k}(end/3:5*end/6,1:end/2+5,i)'-273.15),colorbar, caxis([0 30]), shading flat
        subplot(212),pcolor(X(end/3:5*end/6),Y(1:end/2+5),sstc0(end/3:5*end/6,1:end/2+5,i)'-273.15),colorbar, caxis([0 30]), shading flat
        title(num2str(i))
        pause(0.5)
    end

end
clear sstc11 sstc0
% close all
%%
% --- Guardamos nuevas anomalias ---
% path2 = '/home/danveloso/Documentos/Fondecyt_Garreaud/SPEEDY_SST-SIC/';
for k=1:length(calenton)
    sstc2{k}=flipdim(sstc1{k},2);   %flips Y dim
    sstc2{k}=permute(sstc2{k},[3 2 1]);
    sstc2{k}=reshape(sstc2{k},[4608,12]);
    dummy{k}=ones(4610,12)*variable(1);   %creates a matrix of 2.5829e-41
    dummy{k}(2:48*96+1,:) = sstc2{k};            %puts good values within matrix
    dummy{k}=dummy{k}(:);

    fidw =fopen([path1,'clim/sst_clim6110Hadisst_' num2str(calenton(k)) '_SSWPexp.t30.sea.grd'],'w')
    count{k} = fwrite(fidw,dummy{k},'float','b');
end


