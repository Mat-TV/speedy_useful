#!/bin/bash

# ---------------------------------------------------------------------------------------
# Script para procesar las diferencias entre simulaciones
# Se ingresan las dos simulaciones a comparar y luego se imprimen las diferencias en orden
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
echo "---------- Descripciones ----------"
diff $DIRI"/description_"$exp".txt" $DIRO"/description_"$num".txt"

echo "---------- Comparación del .ctl ----------"
diff $DIRI"/attm"$exp".ctl" $DIRO"/attm"$num".ctl"

echo "---------- Comparación de los años ----------"
echo "Años simulados en exp_$exp:"
ls $DIRI/attm$exp*.grd | wc -l
echo "Años simulados en exp_$num:"
ls $DIRO/attm$num*.grd | wc -l

echo "(En construcción) poner aquí algo que lea los años diferentes entre los 2 exp..."


echo " ---------------------------------------------------"
echo "               FIN DEL PROGRAMA "
echo " ---------------------------------------------------"
