#!/bin/bash

# ---------------------------------------------------------------------------------------
#Script para procesar varios archivos .grd, y convertirlos a .nc
# Daniel Veloso - 2018/01/29
# ---------------------------------------------------------------------------------------

#Primero definimos variables
DIRI='/home/danveloso/Documentos/Fondecyt_Garreaud/SSWP0.0' #Directorio con archivos .grd (Disco Duro 1)
DIRO='/media/danveloso/TOSHIBA EXT/NetCDF-SPEEDY/T30/Exp_SSWP_0.0' #Directorio de salida (Disco Duro 2)
exp=({721..750}) #Nombre de los directorios con los archivos a transformar

# --- Inicio del Programa ---

echo "Experimentos a procesar... ${exp[*]}"
echo "Total de experimentos a procesar... ${#exp[*]}"

for i in ${exp[*]}
do

   echo "-----------------------------------------------"
   echo "          Procesando exp $i"
   echo "-----------------------------------------------"

   # Concatenamos todos los archivos .grd de cada año en un solo .grd
   echo "Concatenando archivos del exp $i ..."
   cd "$DIRI"/exp_${i}/
   cat attm${i}_*.grd > "${DIRO}"/attm${i}_19602017.grd

   #Cambiamos el nombre del ^DSET del .ctl respectivo
   more attm${i}.ctl | sed 's/%y4/19602017/g' > "${DIRO}"/attm${i}_19602017.ctl

   # Hacemos conversion de .grd a.nc
   echo "Haciendo conversión a NetCDF del exp $i ...	"
   cd "$DIRO"
   cdo -f nc import_binary attm${i}_19602017.ctl attm${i}_19602017.nc
   rm attm${i}_19602017.ctl
   rm attm${i}_19602017.grd

done

echo " ---------------------------------------------------"
echo "                  FIN DEL PROGRAMA "
echo " ---------------------------------------------------"

