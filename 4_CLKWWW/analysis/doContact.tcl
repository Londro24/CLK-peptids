# doContact.tcl

if { [llength $argv] < 4 } {
    puts "Error: Faltan archivos."
    puts "Uso: -args <peptido> <psf> <dcd> <replica>"
    quit
}

set peptido     [lindex $argv 0]
set psf_file    [lindex $argv 1]
set dcd_file    [lindex $argv 2]
set replica     [lindex $argv 3]

puts "----------------------------------------------------"
puts "Analisis BINARIO (Contacto Indol TRP) para: ${peptido}.r${replica}"
puts "----------------------------------------------------"

set psf_file "../build/$psf_file"
set dcd_file "../output/03_prod_${peptido}_r${replica}.1.dcd"
set out_file "data/03_contact_${peptido}_r${replica}.dat" 

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

# =================================================================
# PARAMETROS
# =================================================================
# CUTOFFS (Distancias de corte)
set cut_C 3.5;   # Distancia estandar para Azufre-Oro (puedes ajustar a 3.0 o 3.5)
set cut_W 3.5;   # Distancia para el anillo Indol

# UMBRAL PARA TRIPTOFANO (IMPORTANTE)
# 1 = Basta con que 1 atomo del anillo toque para contar como contacto (Sensible).
# 7 = Requiere que casi todo el anillo este plano sobre la superficie (Estricto).
set min_atoms_W 5; 
# =================================================================

# =============================================================
# DETECCION DE TRIPTOFANOS
# =============================================================
set sel_all_trp [atomselect top "resname TRP"]
set trp_resids [lsort -integer -unique [$sel_all_trp get resid]]
$sel_all_trp delete
set num_trps [llength $trp_resids]

puts " -> TRP IDs encontrados: $trp_resids"

# =============================================================
# CABECERA
# =============================================================
set header "Frame\tCys_S"

for {set k 1} {$k <= $num_trps} {incr k} {
    set real_resid [lindex $trp_resids [expr $k-1]]
    append header "\tTRP${real_resid}"
}

set out_file [open $out_file w]
puts $out_file $header

# =============================================================
# SELECCIONES ESTATICAS
# =============================================================
set sel_oro "resname IAUM"
set sel_text_C "resname CYS and name SG"

set num_frames [molinfo top get numframes]

# 3. BUCLE DE FRAMES
for {set i 0} {$i < $num_frames} {incr i} {
    
    animate goto $i
    
    # --- A. CYS (Azufre) ---
    set sel [atomselect top "($sel_text_C) and within $cut_C of ($sel_oro)"]
    set n_C [$sel num]
    $sel delete
    set bin_C [expr {$n_C > 0 ? 1 : 0}]

    # --- B. TRP (Anillo Indol) ---
    set trp_values_str ""
    
    foreach resid $trp_resids {
        set sel_text_W_ind "resname TRP and resid $resid and name CG CD1 CD2 NE1 CE2 CE3 CZ2 CZ3 CH2"
        
        # Contamos cuantos de estos 9 atomos estan cerca del oro
        set sel [atomselect top "($sel_text_W_ind) and within $cut_W of ($sel_oro)"]
        set n_W [$sel num]
        $sel delete
        
        # LOGICA BINARIA
        if { $n_W >= $min_atoms_W } {
            set bin_W 1
        } else {
            set bin_W 0
        }
        
        append trp_values_str "\t$bin_W"
    }

    # Escribir linea
    puts $out_file "$i\t$bin_C${trp_values_str}"

    if { [expr $i % 100] == 0 } { puts "    ...frame $i" }
}

close $out_file
mol delete top
puts "Terminado. Archivo generado correctamente."
quit