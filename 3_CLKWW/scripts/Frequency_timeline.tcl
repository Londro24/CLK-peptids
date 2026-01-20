# ==================================================================
# SCRIPT TCL: MATRIZ DE DISTANCIAS (COM a COM)
# ==================================================================

if { [llength $argv] < 4 } {
    puts "ERROR: Faltan argumentos (PSF DCD ORO OUTPUT)."
    quit
}

set psf_file    [lindex $argv 0]
set dcd_file    [lindex $argv 1]
set gold_sel    [lindex $argv 2]
set out_file    [lindex $argv 3]

# Cargar archivos
mol new $psf_file
mol addfile $dcd_file waitfor all

set outfile [open $out_file w]
set num_frames [molinfo top get numframes]
set protein [atomselect top "protein"]
set resid_list [lsort -unique -integer [$protein get resid]]
$protein delete

# ------------------------------------------------------------------
# 1. ENCABEZADO
# ------------------------------------------------------------------
puts -nonewline $outfile "Frame\tTime(ns)"

foreach r $resid_list {
    set temp [atomselect top "protein and resid $r"]
    set rname [lindex [$temp get resname] 0]
    $temp delete
    # Encabezado tipo: MET_1  GLY_2 ...
    puts -nonewline $outfile "\t${rname}_${r}"
}
puts $outfile ""

# ------------------------------------------------------------------
# 2. BUCLE DE CÁLCULO
# ------------------------------------------------------------------
puts "Calculando distancias..."

for {set f 0} {$f < $num_frames} {incr f} {
    animate goto $f
    
    # Tiempo aproximado (asumiendo 20ns totales, ajusta si es necesario)
    set time_ns [expr ($f * 20.0) / $num_frames]
    puts -nonewline $outfile "$f\t[format "%.2f" $time_ns]"

    # Definir la selección de oro y calcular su CENTRO DE MASA (COM)
    set surface [atomselect top "$gold_sel"]
    set com_gold [measure center $surface weight mass]
    
    # Bucle por residuo
    foreach r $resid_list {
        set res_sel [atomselect top "protein and resid $r"]
        
        # Calcular el CENTRO DE MASA del residuo
        set com_res [measure center $res_sel weight mass]
        
        # Calcula solo la diferencia en el eje Z (altura)
        set dist [expr abs([lindex $com_res 2] - [lindex $com_gold 2])]
        
        # Guardar distancia con 2 decimales
        puts -nonewline $outfile [format "\t%.2f" $dist]
        
        $res_sel delete
    }
    $surface delete
    
    puts $outfile "" ;# Fin de la fila
    
    if {$f % 100 == 0} { puts " -> Frame $f procesado" }
}

close $outfile
puts "Terminado."
quit