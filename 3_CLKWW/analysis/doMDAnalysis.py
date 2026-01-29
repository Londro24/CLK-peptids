import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import warnings
from MDAnalysis.analysis import distances
from pathlib import Path
import sys

# ==========================================
# CONFIGURACIÓN DE USUARIO
# ==========================================
if len(sys.argv) > 1:
    peptido = sys.argv[1] 
else:
    print("Error: Debes indicar el péptido como argumento.")
    sys.exit(1)

replicas = [1, 2, 3]
psf_file = f'../build/repart_{peptido}.psf'
ventana_suavizado = 20 
output_dir_data = "data/"
output_dir_graphics = "graphics/"
dt = 0.04  # ns por frame

TIEMPO_EQUILIBRIO = 20.0 

# Asegurar carpetas
Path(output_dir_data).mkdir(parents=True, exist_ok=True)
Path(output_dir_graphics).mkdir(parents=True, exist_ok=True)

# ==========================================
# ESTILO VISUAL
# ==========================================
def configurar_estilo():
    warnings.filterwarnings("ignore", message=".*DCDReader currently makes independent timesteps.*")
    try:
        sns.set_theme(style="whitegrid", context="paper", font_scale=1.2)
        sns.set_palette("deep")
    except Exception:
        plt.style.use("seaborn-whitegrid")

    plt.rcParams.update({
        "figure.dpi": 300,
        "lines.linewidth": 1.5,
        "axes.grid": True,
        "grid.alpha": 0.3,
        "legend.frameon": True,
        "legend.fancybox": True,
    })

configurar_estilo()

dfs_list = []
labels_detectados = [] # Para guardar nombres de residuos (CYS2, TRP4...)

print(f"Iniciando análisis para {peptido}...")

# =============================================================================
# 1. CÁLCULO DE DATOS
# =============================================================================
for r in replicas:
    print(f"\n--- Procesando Réplica {r} ---")
    dcd_file = f'data/00_merge_{peptido}_r{r}.dcd'
    
    if not Path(psf_file).exists() or not Path(dcd_file).exists():
        print(f"  [AVISO] Archivos no encontrados para R{r}. Saltando...")
        continue
    
    try:
        u = mda.Universe(psf_file, dcd_file)
    except Exception as e:
        print(f"  Error loading universe: {e}")
        continue

    oro = u.select_atoms("resname IAUM")
    peptido_sel = u.select_atoms("protein")
    
    # Identificar residuos de interés (TRP y CYS) - Solo cadenas laterales
    residuos_interes = u.select_atoms("resname TRP CYS").residues
    sidechain_groups = {}
    
    for res in residuos_interes:
        label = f"{res.resname}{res.resid}"
        if label not in labels_detectados:
            labels_detectados.append(label)
        # Excluir backbone
        sidechain_groups[label] = res.atoms.select_atoms("not name N CA C O")

    if len(oro) == 0:
        print("  [ERROR] No se encontró 'resname IAUM'.")
        continue

    # Bucle de Trayectoria
    data_rows = []
    for ts in u.trajectory:
        t_ns = ts.frame * dt
        
        # A. Distancia Global
        d_global = np.min(distances.distance_array(peptido_sel.positions, oro.positions))
        
        row = {"Tiempo": t_ns, "Replica": f"R{r}", "Distancia_Global": d_global}
        
        # B. Distancia por Residuo
        for label, ag in sidechain_groups.items():
            if len(ag) > 0:
                row[label] = np.min(distances.distance_array(ag.positions, oro.positions))
            else:
                row[label] = np.nan

        data_rows.append(row)

    # Crear DataFrame y Suavizar
    df_rep = pd.DataFrame(data_rows)
    
    # Suavizar todas las columnas numéricas
    cols_num = [c for c in df_rep.columns if c not in ["Tiempo", "Replica"]]
    for col in cols_num:
        df_rep[f"{col}_Smooth"] = df_rep[col].rolling(window=ventana_suavizado, center=True, min_periods=1).mean()

    # Guardar .DAT
    file_dat = f"{output_dir_data}05_distancias_{peptido}_r{r}.dat"
    df_rep.to_csv(file_dat, sep="\t", index=False, float_format="%.3f")
    dfs_list.append(df_rep)

# =============================================================================
# 2. GENERACIÓN DE GRÁFICOS
# =============================================================================
if dfs_list:
    df_total = pd.concat(dfs_list, ignore_index=True)
    labels_detectados = sorted(list(set(labels_detectados)))
    
    # -------------------------------------------------------------------------
    # GRÁFICO 1: 05_m_ (Completo 0-70ns, Absoluto)
    # -------------------------------------------------------------------------
    plt.figure(figsize=(10, 6))
    
    # Sombreado de Zonas
    plt.axvspan(0, TIEMPO_EQUILIBRIO, color='tab:red', alpha=0.1, zorder=0)
    plt.axvspan(TIEMPO_EQUILIBRIO, 70, color='tab:green', alpha=0.1, zorder=0)

    # Datos Promedio Global
    sns.lineplot(
        data=df_total, 
        x="Tiempo", 
        y="Distancia_Global_Smooth", 
        errorbar="sd", 
        color="darkblue",
        label="Promedio Global"
    )

    # Línea Vertical
    plt.axvline(x=TIEMPO_EQUILIBRIO, color='black', linestyle='--', alpha=0.5)

    # Textos Alineados a la línea de 20ns
    plt.text(TIEMPO_EQUILIBRIO - 0.5, 24, "Equilibrado", ha='right', va='top', fontsize=11, fontweight='bold', color='darkred')
    plt.text(TIEMPO_EQUILIBRIO + 0.5, 24, "Producción", ha='left', va='top', fontsize=11, fontweight='bold', color='darkgreen')

    plt.xlim(0, 70)
    plt.ylim(0, 25) 
    plt.xlabel("Tiempo de Simulación (ns)")
    plt.ylabel("Distancia Mínima (Å)")
    plt.title(f"Distancia Global Promedio (Completa) - {peptido}")
    plt.legend(loc='upper right')
    plt.tight_layout()
    plt.savefig(f"{output_dir_graphics}05_Distancia_Full_{peptido}.png")
    plt.close()
    print("Generado: 05_Distancia_Full (0-70ns)")

    # -------------------------------------------------------------------------
    # FILTRAR DATOS PARA PRODUCCIÓN Y AJUSTAR TIEMPO (0-50ns)
    # -------------------------------------------------------------------------
    df_prod = df_total[df_total["Tiempo"] >= TIEMPO_EQUILIBRIO].copy()
    
    # === CAMBIO CLAVE: RESTAR EL TIEMPO DE EQUILIBRIO ===
    # Esto hace que 20ns se convierta en 0ns, y 70ns en 50ns
    df_prod["Tiempo"] = df_prod["Tiempo"] - TIEMPO_EQUILIBRIO

    if not df_prod.empty:
        # ---------------------------------------------------------------------
        # GRÁFICO 2: 05_p_ (Solo Producción, Replicas Individuales)
        # ---------------------------------------------------------------------
        plt.figure(figsize=(10, 6))
        sns.lineplot(
            data=df_prod, 
            x="Tiempo", 
            y="Distancia_Global_Smooth", 
            hue="Replica", 
            palette="deep",
            linewidth=1.5
        )
        plt.xlim(0, 50) # Aseguramos eje exacto
        plt.ylim(0, 5)
        plt.xlabel("Tiempo de Producción (ns)")
        plt.ylabel("Distancia Mínima (Å)")
        plt.title(f"Distancia por Réplica (Producción) - {peptido}")
        plt.legend(loc='upper right', title="Réplica")
        plt.tight_layout()
        plt.savefig(f"{output_dir_graphics}05_Distancia_Replicas_{peptido}.png")
        plt.close()
        print("Generado: 05_Distancia_Replicas (0-50ns)")

        # ---------------------------------------------------------------------
        # GRÁFICO 3: 05_p_ (Solo Producción, Promedio + StdDev)
        # ---------------------------------------------------------------------
        plt.figure(figsize=(10, 6))
        sns.lineplot(
            data=df_prod, 
            x="Tiempo", 
            y="Distancia_Global_Smooth", 
            errorbar="sd", 
            color="teal", 
            label="Promedio ± DE"
        )
        plt.xlim(0, 50)
        plt.ylim(0, 5)
        plt.xlabel("Tiempo de Producción (ns)")
        plt.ylabel("Distancia Mínima (Å)")
        plt.title(f"Distancia Promedio (Producción) - {peptido}")
        plt.legend(loc='upper right')
        plt.tight_layout()
        plt.savefig(f"{output_dir_graphics}05_Distancia_Promedio_{peptido}.png")
        plt.close()
        print("Generado: 05_Distancia_Promedio (0-50ns)")

        # ---------------------------------------------------------------------
        # GRÁFICO 4: COMPARACIÓN RESIDUOS (Producción, Promedios + StdDev)
        # ---------------------------------------------------------------------
        if labels_detectados:
            plt.figure(figsize=(10, 6))
            
            for label in labels_detectados:
                col_smooth = f"{label}_Smooth"
                if col_smooth in df_prod.columns:
                    sns.lineplot(
                        data=df_prod, 
                        x="Tiempo", 
                        y=col_smooth, 
                        errorbar="sd", 
                        linewidth=2,
                        label=label,
                        alpha=0.8
                    )
            
            plt.xlim(0, 50)
            plt.ylim(0, 10)
            plt.xlabel("Tiempo de Producción (ns)")
            plt.ylabel("Distancia Mínima (Å)")
            plt.title(f"Comparación Residuos con DE (Producción) - {peptido}")
            plt.legend(loc='upper right', title="Residuo")
            plt.tight_layout()
            plt.savefig(f"{output_dir_graphics}05_Comparacion_Residuos_{peptido}.png")
            plt.close()
            print("Generado: 05_Comparacion_Residuos (0-50ns)")
        else:
            print("AVISO: No se encontraron residuos (TRP/CYS) para el gráfico 4.")

else:
    print("No se generaron datos.")

print("\nProceso finalizado.")