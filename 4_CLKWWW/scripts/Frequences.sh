#!/bin/bash

# ==========================================
# CONFIGURACIÓN DEL USUARIO
# ==========================================
# Nombre de tus archivos de entrada
# Se ejecuta desde /Peptido/scripts/
peptido="4_CLKWWW"
psf_file="../build/repart_${peptido}.psf"
dcd_file="../output/02_eq_${peptido}.0.dcd"

# Parámetros del análisis
gold_selection="type IAU"   ;# Como seleccionas tu oro en VMD
cutoff="3.5"                ;# Distancia en Angstroms para considerar contacto
output_name_1="../output/freq_contacto_${peptido}.dat"
output_name_2="../output/timeline_contacto_${peptido}.dat"

# ==========================================
# EJECUCIÓN
# ==========================================
# Llama a VMD en modo texto y pasa las 5 variables al script Tcl
vmd -dispdev text -e Frequency_contact.tcl -args "${psf_file}" "${dcd_file}" "${gold_selection}" "${cutoff}" "${output_name_1}"

echo "Analisis de frecuencia finalizado. Revisa el archivo: ${output_name_1}"

vmd -dispdev text -e Frequency_timeline.tcl -args "${psf_file}" "${dcd_file}" "${gold_selection}" "${output_name_2}"

echo "Analisis de timeline finalizado. Archivo generado: ${output_name_2}"