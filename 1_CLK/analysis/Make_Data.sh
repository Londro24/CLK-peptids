#!/bin/bash
#doALL.sh

# --- CONFIGURACION GENERAL ---
peptido=$1
psf="repart_${peptido}.psf"
# -----------------------------

echo "Iniciando analisis automatico de RMSD para replicas 1 a 3..."

# Ciclo para procesar replicas 1 a 3
for replica in {1..3}
do
    echo "-----------------------------------"
    echo "Procesando REPLICA: $replica"
    
    # Definimos el DCD dinamicamente segun el numero de replica
    dcd="00_merge_${peptido}_r${replica}.dcd"

    # Ejecutar VMD
    vmd -dispdev text -e doRMSD_merge.tcl -args $peptido $psf $dcd $replica
    vmd -dispdev text -e doRMSD_prod.tcl -args $peptido $psf $dcd $replica
    #vmd -dispdev text -e doRMSF.tcl -args $peptido $psf $dcd $replica
    vmd -dispdev text -e doContact.tcl -args $peptido $psf $dcd $replica
done

echo "-----------------------------------"
echo "Analisis completado para todas las replicas."