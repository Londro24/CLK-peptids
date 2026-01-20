# Seleccionamos la molécula activa
set molid top

# --- PASO CLAVE: BORRAR LA VISUALIZACIÓN POR DEFECTO ---
# Borramos la rep 0 (Lines/Name) que VMD crea al abrir
mol delrep 0 $molid

# ---------------------------------------------------------
# AHORA AÑADIMOS TUS 3 VISTAS
# ---------------------------------------------------------

# 1. AIUM | VDW | Opaque | ColorID 32
mol selection {resname IAUM}
mol representation VDW 1.0 12.0
mol material Opaque
mol color ColorID 32
mol addrep $molid

# 2. Water | QuickSurf | Glass1 | ColorID 10
mol selection {water}
mol representation QuickSurf 1.0 0.5 1.0 1.0
mol material Glass1
mol color ColorID 10
mol addrep $molid

# 3. Protein | Licorice | Opaque | Name
mol selection {protein}
mol representation Licorice 0.3 12.0 12.0
mol material Opaque
mol color Name
mol addrep $molid
