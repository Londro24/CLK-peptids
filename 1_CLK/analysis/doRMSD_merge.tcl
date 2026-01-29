# doRMSD_merge.tcl
# Script para calcular RMSD de todo el peptido

if { [llength $argv] < 4 } {
    puts "Error: Faltan archivos."
    puts "Uso: -args <peptido> <psf> <dcd> <replica>"
    quit
}

set peptido     [lindex $argv 0]
set psf_file    [lindex $argv 1]
set dcd_file    [lindex $argv 2]
set replica     [lindex $argv 3]

puts "------------------------------------------------"
puts "Iniciando analisis para: ${peptido}.r${replica} "
puts "------------------------------------------------"

set psf_file "../build/$psf_file"
set dcd_file "data/$dcd_file"
set out_file "data/01_RMSD_m_${peptido}_r${replica}.dat"

# Cargar molécula
if { [catch {mol new $psf_file type psf} err] } {
    puts "ERROR al cargar PSF: $err"
    quit
}
if { [catch {mol addfile $dcd_file type dcd waitfor all} err] } {
    puts "ERROR al cargar DCD: $err"
    mol delete top
    quit
}

# SELECCION
# "protein" selecciona todos los aminoacidos (backbone + sidechains)
set sel_definition "protein"

# Verificacion para el usuario
set check_sel [atomselect top $sel_definition]
if { [$check_sel num] == 0 } {
    puts "ERROR CRITICO: VMD no detecto ninguna proteina/peptido."
    puts "Asegurate de que tu PSF esta correcto."
    quit
}

puts "------------------------------------------------"
puts " Total de atomos seleccionados: [$check_sel num]"
puts "------------------------------------------------"

# Configuracion de Referencia (Frame 0)
set reference [atomselect top $sel_definition frame 0]
set compare   [atomselect top $sel_definition]
set num_steps [molinfo top get numframes]

# Guardar
set out_file [open $out_file w]

puts "Calculando RMSD frame a frame..."

# Bucle de calculo
for {set frame 0} {$frame < $num_steps} {incr frame} {
    $compare frame $frame
    
    # A. Alineamiento (Fitting)
    # Superpone el peptido actual con la referencia para eliminar rotacion/traslacion
    set trans_mat [measure fit $compare $reference]
    $compare move $trans_mat
    
    # B. Calculo RMSD
    set rmsd [measure rmsd $compare $reference]
    
    # Guardar: Frame RMSD
    puts $out_file "$frame \t $rmsd"
}

close $out_file
$check_sel delete
$reference delete
$compare delete

puts "Listo. Datos en '${out_file}'"
quit