#############################################################
## Pipeline de Preparación de Sistema para NAMD            ##
#############################################################

# 0. Crear Inicio
# Genera una carpeta "build" y copia el PDB original allí.
# Ademas de colocar la lamina generada por CHARMM-GUI (pdb y psf) en la carpeta "build"

# 1. Validación de Argumentos
# Uso esperado: vmd -dispdev text -e Build_System.tcl -args <nombre_input> <nombre_final>
if {[llength $argv] < 3} {
    puts "Error: Argumentos insuficientes."
    puts "Uso: vmd -dispdev text -e Build_System.tcl -args <name_pdb_original> <name_surface> <directory>"
    exit
}

# 2. Definir Variables
set name [lindex $argv 0]
set surface [lindex $argv 1]
set dir [lindex $argv 2]
set structures [list $name $surface]

# Cargar paquetes necesarios
package require psfgen
package require solvate
package require autoionize
package require pbctools


#############################################################
## PASO 1: Build PSF/PDB y agregar Hidrógenos              ##
#############################################################
puts "--- PASO 1: Generando PSF y agregando Hidrógenos ---"

topology ../parameters/top_all36_prot.rtf

# Aliases estándar para CHARMM
pdbalias residue HIS HSD
pdbalias atom SER HG HG1
pdbalias atom ILE CD1 CD
pdbalias atom CYS HG HG1
pdbalias atom ILE 1HD1 HD1
pdbalias atom ILE 2HD1 HD2
pdbalias atom ILE 3HD1 HD3

resetpsf

segment A { 
    first NTER
    last CTER
    pdb ${dir}/${name}.pdb 
}
coordpdb ${dir}/${name}.pdb A
regenerate angles dihedrals
guesscoord

writepdb ${dir}/1-${name}_h.pdb
writepsf ${dir}/1-${name}_h.psf

#############################################################
## PASO 2: Juntar estructuras (Peptido + Lamina)           ##
#############################################################
puts "--- PASO 2: Juntando estructuras ---"

resetpsf

foreach st $structures {
    if {$st eq $name} {
        mol load psf ${dir}/1-${st}_h.psf pdb ${dir}/1-${st}_h.pdb
    } 
    if {$st eq $surface} {
        mol load psf ${dir}/${st}.psf pdb ${dir}/${st}.pdb
    }

    set sel [atomselect top "all"]
    set center [measure center $sel]
    set m1 [trans origin $center]
    $sel move $m1
    
    if {$st eq $name} {
        $sel move [transaxis x 90]
        set distancia 30.0
        $sel moveby [list 0 0 $distancia]
    } 
    
    $sel writepdb ${dir}/2-${st}_c.pdb
    $sel writepsf ${dir}/2-${st}_c.psf

    mol delete all
}

foreach st $structures {
    readpsf ${dir}/2-${st}_c.psf
    coordpdb ${dir}/2-${st}_c.pdb
}

writepsf ${dir}/2-sist-${name}.psf
writepdb ${dir}/2-sist-${name}.pdb

#############################################################
## PASO 3: Solvatación y centrado (Caja de Agua)           ##
#############################################################
puts "--- PASO 3: Solvatando el sistema ---"

# centrado del sistema
mol load psf ../${dir}/2-sist-${name}.psf pdb ${dir}/2-sist-${name}.pdb
set sel [atomselect top "all"]
set center [measure center $sel]
set m1 [trans origin $center]
$sel move $m1

writepsf ${dir}/3-sist-${name}_c.psf
writepdb ${dir}/3-sist-${name}_c.pdb

mol delete all

# solvatación
solvate ${dir}/3-sist-${name}_c.psf ${dir}/3-sist-${name}_c.pdb -o ${dir}/3-sist-${name}_solv -s WT -b 2.4 -x 0 +x 0 -y 0 +y 0 -z 10 +z 10

mol delete all

#############################################################
## PASO 4: Ionización (Neutralizado)                       ##
#############################################################
puts "--- PASO 4: Ionizando el sistema ---"

# ionización
autoionize -psf ${dir}/3-sist-${name}_solv.psf -pdb ${dir}/3-sist-${name}_solv.pdb -o ${dir}/final-sist-${name}_ion -neutralize

#############################################################
## PASO 5: Generar archivo XSC para NAMD                   ##
#############################################################
puts "--- PASO 5: Generando archivo .xsc y archivos de minimización ---"

# Crear archivo XSC
pbc writexst ${dir}/init.xsc -now

puts "--- PROCESO COMPLETADO EXITOSAMENTE ---"
exit
