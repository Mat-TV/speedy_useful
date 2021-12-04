#!/bin/bash

# ---------------------------------------------------------------------------------------
# Script para describir las simulaciones
# Se especifica la simulación a describir y luego simplemente se describe.
# Mat Troncoso - 4/12/2021
# ---------------------------------------------------------------------------------------

# --- Preámbulo ---
#exp=106
read -p "Número de simulación a describir: " exp

DIRI='/home/speedy/output/exp_'$exp #Directorio del experimento

# --- Inicio del Programa ---
echo "---------- Descripción ----------"
read -p "Descripción: " text

echo "---------- Describiendo ----------"
#sed -i '3s/^/DESCRIPTION    '"${text}"' \n/' $DIRI/description_$exp.txt
echo "DESCRIPTION    ${text}" >> $DIRI/description_$exp.txt

echo " ---------------------------------------------------"
echo "               FIN DEL PROGRAMA "
echo " ---------------------------------------------------"
