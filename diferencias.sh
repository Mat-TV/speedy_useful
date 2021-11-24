#!/bin/bash

# ---------------------------------------------------------------------------------------
# Script para procesar las diferencias entre simulaciones
# Mat Troncoso - 8/11/2021
# ---------------------------------------------------------------------------------------

# --- Preámbulo ---
#exp=106
read -p "Número de simulación a comparar: " exp
read -p "Número de simulación con la que comparar: " num

DIRI='/home/speedy/output/exp_'$exp #Directorio de origen
DIRO=${DIRI}'/../exp_'${num} #Directorio de la comparación


# --- Inicio del Programa ---
   echo "-----------------------------------------------"
   echo "        Comparando exp_$exp y exp_$num"
   echo "-----------------------------------------------"
echo "---------- Comparación del .ctl ----------"
diff $DIRI"/attm"$exp".ctl" $DIRO"/attm"$num".ctl"

echo "---------- Comparación de los años ----------"
echo "(en construcción)"
echo "poner aquí algo que cuente los elementos .grd y lea los años diferentes entre los 2 exp..."



echo " ---------------------------------------------------"
echo "               FIN DEL PROGRAMA "
echo " ---------------------------------------------------"

