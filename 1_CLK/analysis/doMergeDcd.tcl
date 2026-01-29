# doMergeDcd.tcl
# Script para juntar archivos .dcd

if { [llength $argv] < 5 } {
    puts "Error: Faltan archivos."
    puts "Uso: -args <peptido> <replica> <psf> <dcd_eq> <dcd_prod>"
    quit
}

set peptido         [lindex $argv 0]
set replica         [lindex $argv 1]
set psf_file        [lindex $argv 2]
set dcd_eq_file     [lindex $argv 3]
set dcd_prod_file   [lindex $argv 4]

# 1. Cargar archivos
set psf_file "../build/$psf_file"
set dcd_eq_file "../output/02_eq_${peptido}.0.dcd"
set dcd_prod_file "../output/03_prod_${peptido}_r${replica}.1.dcd"
set out_file "data/00_merge_${peptido}_r${replica}.dcd"

puts "\nProcesando Replica $replica..."

# Cargar molécula
if { [catch {mol new $psf_file type psf} err] } {
    puts "ERROR al cargar PSF: $err"
    continue
}
if { [catch {mol addfile $dcd_eq_file type dcd waitfor all} err] } {
    puts "ERROR al cargar DCD: $err (Quizas no existe?)"
    mol delete top
    continue
}
if { [catch {mol addfile $dcd_prod_file type dcd waitfor all} err] } {
    puts "ERROR al cargar DCD: $err (Quizas no existe?)"
    mol delete top
    continue
}

puts "--- Iniciando proceso VMD para $peptido Rep $replica ---"

# 3. Seleccion
set seleccion [atomselect top "all"]

# 4. Guardar archivo concatenado
puts "Guardando $out_file ..."

animate write dcd $out_file beg 0 end -1 waitfor all sel $seleccion

puts "Terminado."
exit