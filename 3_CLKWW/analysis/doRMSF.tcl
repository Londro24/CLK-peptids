# doRMSF.tcl
# Script para calcular RMSF de los residuos

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
set out_file "data/02_RMSF_${peptido}_r${replica}.dat"

# Cargar molecula
if { [catch {mol new $psf_file type psf} err] } {
    puts "ERROR al cargar PSF: $err"
    quit
}
if { [catch {mol addfile $dcd_file type dcd waitfor all} err] } {
    puts "ERROR al cargar DCD: $err"
    mol delete top
    quit
}

set molid [molinfo top]
set num_frames [molinfo $molid get numframes]

# Alineamiento
puts "Alineando $num_frames frames..."
set ref [atomselect $molid "backbone" frame 0]
set sel [atomselect $molid "backbone"]
set all [atomselect $molid "all"]

for {set i 0} {$i < $num_frames} {incr i} {
    $sel frame $i
    $all frame $i
    set mat [measure fit $sel $ref]
    $all move $mat
}

# Calculo RMSF (Sobre Carbonos Alfa)
set ca_sel [atomselect $molid "protein and name CA"]
set rmsf_vals [measure rmsf $ca_sel]

set resids   [$ca_sel get resid]
set resnames [$ca_sel get resname]

# Guardar todos los residuos
set outfile [open $out_file "w"]

# Cabecera
puts $outfile "# Reporte RMSF - Replica $replica"
puts $outfile "# Resname Resid RMSF"

# Iterar sobre las listas. Como quitamos el 'if', guarda todo.
foreach name $resnames id $resids val $rmsf_vals {
    # Formato alineado:
    # %-8s : Nombre (ej. TRP, ALA)
    # %-6d : ID Residuo
    # %.4f : Valor RMSF
    puts $outfile [format "%-8s %-6d %.4f" $name $id $val]
}

close $outfile
puts "Listo. Datos en '${out_file}'"

# Limpieza
mol delete $molid
quit