#!/bin/bash

# --- DEFINIR VARIABLE GLOBAL ---
PEPTIDO="4_CLKWWW"

echo "=== Iniciando proceso para: $PEPTIDO ==="

# 1. Ejecutar los scripts .sh (pasando la variable como $1)
echo "-> Ejecutando script bash 1..."
bash Merge.sh "$PEPTIDO"

echo "-> Ejecutando script bash 2..."
bash Make_Data.sh "$PEPTIDO"

# 2. Ejecutar los scripts .py (pasando la variable como argumento)
echo "-> Ejecutando script python 1..."
python3 doGraphics.py "$PEPTIDO"

echo "-> Ejecutando script python 2..."
python3 doMDAnalysis.py "$PEPTIDO"

echo "=== Todo finalizado ==="