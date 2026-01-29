#!/bin/bash
# doALL.sh

# --- CONFIGURACION GENERAL ---
peptido=$1
psf="repart_${peptido}.psf"
dcd_eq="02_eq_${peptido}.0.dcd"
# -----------------------------

echo "Iniciando analisis automatico de RMSD para replicas 1 a 3..."

# Ciclo para procesar replicas 1 a 3
for replica in {1..3}
do
    echo "-----------------------------------"
    echo "Procesando REPLICA: $replica"
    
    # Definimos el DCD dinamicamente segun el numero de replica
    dcd_prod="03_prod_${peptido}_r${replica}.1.dcd"

    # Ejecutar VMD
    vmd -dispdev text -e doMergeDcd.tcl -args $peptido $replica $psf $dcd_eq $dcd_prod 
done

echo "-----------------------------------"
echo "Equilibrado y Produccion del sistema ${peptido}."