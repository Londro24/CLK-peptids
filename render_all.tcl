# ==========================================
# CONFIGURACIÓN GENERAL
# ==========================================
# Distancia de separación en el eje X (Angstroms)
set spacing 20 
set current_x 0

# Lista de sistemas (carpetas)
set systems {
    "1_CLK"
    "2_CLKW"
    "3_CLKWW"
    "4_CLKWWW"
}

# ==========================================
# BUCLE PRINCIPAL
# ==========================================
foreach sys $systems {
    
    # Rutas a los archivos (según tu estructura de carpetas)
    set psf_file "$sys/build/final-sist-${sys}_ion.psf"
    set pdb_file "$sys/build/final-sist-${sys}_ion.pdb"

    if { [file exists $psf_file] && [file exists $pdb_file] } {
        
        # 1. CARGAR SISTEMA
        mol new $psf_file type psf
        mol addfile $pdb_file type pdb
        mol rename top "$sys"

        # ---------------------------------------------------------
        # APLICAR VISUALIZACIÓN (Tu configuración)
        # ---------------------------------------------------------
        # Borramos la rep 0 (Lines/Name) que VMD crea por defecto
        mol delrep 0 top

        # Vista 1: AIUM | VDW | Opaque | ColorID 32
        mol selection {resname IAUM}
        mol representation VDW 1.0 12.0
        mol material Opaque
        mol color ColorID 32
        mol addrep top

        # Vista 2: Water | QuickSurf | Glass1 | ColorID 10
        mol selection {water}
        mol representation QuickSurf 1.0 0.5 1.0 1.0
        mol material Glass1
        mol color ColorID 10
        mol addrep top

        # Vista 3: Protein | Licorice | Opaque | Name
        mol selection {protein}
        mol representation Licorice 0.3 12.0 12.0
        mol material Opaque
        mol color Name
        mol addrep top
        # ---------------------------------------------------------

        # 2. MOVER SISTEMA (Separación en X)
        if {$current_x != 0} {
            set sel [atomselect top all]
            $sel moveby [list $current_x 0 0]
            $sel delete
        }
        
        puts ">> Cargado: $sys | Desplazado X: $current_x A"

        # Actualizar posición para el siguiente
        set current_x [expr {$current_x + $spacing}]

    } else {
        puts "ERROR CRÍTICO: No se encontró $psf_file"
    }
}

# ==========================================
# AJUSTES FINALES DE ESCENA
# ==========================================
display resetview
axes location off
color Display Background white