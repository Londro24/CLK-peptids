# ==================================================================
# SCRIPT TCL: CÁLCULO DE FRECUENCIA DE CONTACTO (Recibe argumentos)
# ==================================================================

# 1. LEER ARGUMENTOS DESDE EL SCRIPT BASH
# VMD guarda los argumentos de -args en la lista $argv
if { [llength $argv] < 5 } {
    puts "ERROR: Faltan argumentos. Se requieren: PSF DCD SEL_ORO CUTOFF OUTPUT"
    quit
}

set psf_file    [lindex $argv 0]
set dcd_file    [lindex $argv 1]
set gold_sel    [lindex $argv 2]
set cutoff      [lindex $argv 3]
set out_file    [lindex $argv 4]

puts "------------------------------------------------"
puts "Iniciando analisis con:"
puts "  PSF: $psf_file"
puts "  DCD: $dcd_file"
puts "  Oro: $gold_sel"
puts "  Cutoff: $cutoff A"
puts "------------------------------------------------"

# 2. CARGAR TRAYECTORIA
mol new $psf_file
mol addfile $dcd_file waitfor all

# 3. PREPARAR VARIABLES
set outfile [open $out_file w]
puts $outfile "# ResidID\tResName\tFrecuencia_Contacto"

set num_frames [molinfo top get numframes]
set protein [atomselect top "protein"]
# Lista unica de IDs de residuos del peptido
set resid_list [lsort -unique -integer [$protein get resid]]
$protein delete

# Inicializar contadores a 0
foreach r $resid_list {
    set contact_counts($r) 0
}

# 4. BUCLE PRINCIPAL (Frames)
puts "Procesando $num_frames frames..."

for {set f 0} {$f < $num_frames} {incr f} {
    animate goto $f
    
    # Progreso en pantalla cada 100 frames
    if {$f % 100 == 0} { puts " -> Frame $f..." }

    # Seleccionamos la superficie en el frame actual
    set surface [atomselect top "$gold_sel"]

    # Bucle por cada residuo
    foreach r $resid_list {
        set res_sel [atomselect top "protein and resid $r"]
        
        # 'measure contacts' devuelve dos listas. Si la primera no esta vacia, hay contacto.
        set contact_check [measure contacts $cutoff $res_sel $surface]
        
        if {[llength [lindex $contact_check 0]] > 0} {
            incr contact_counts($r)
        }
        
        $res_sel delete
    }
    $surface delete
}

# 5. ESCRIBIR RESULTADOS
puts "Escribiendo resultados en $out_file..."

foreach r $resid_list {
    # Sacamos el nombre del residuo (ej. ALA, CYS) para que quede bonito
    set temp [atomselect top "protein and resid $r"]
    set rname [lindex [$temp get resname] 0]
    $temp delete

    # Calculo: (Veces que tocó) / (Total de frames)
    set freq [expr double($contact_counts($r)) / double($num_frames)]
    
    # Formato:  1   CYS   0.95
    puts $outfile [format "%d\t%s\t%.4f" $r $rname $freq]
}

close $outfile
puts "Terminado exitosamente."
quit