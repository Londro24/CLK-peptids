from pathlib import Path
from typing import List
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import re

#=======================
#CONFIGURACIÓN DE USUARIO
#=======================
if len(sys.argv) > 1:
    PEPTIDO = sys.argv[1]
else:
    print("Error: Debes indicar el péptido como argumento.")
    sys.exit(1)

RUTA_DATOS = Path("data")
CARPETA_SALIDA = Path("graphics")
DT = 0.04  # ns por frame
REPLICAS = 3
TIEMPO_INICIO_CONTACTOS = 5.0 # ns para filtrar contactos

# =======================
# ESTILO VISUAL
# =======================
def configurar_estilo():
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
        "axes.titlesize": 14,
        "axes.labelsize": 12,
        "legend.fontsize": 10,
    })

# =======================
# CARGA DE DATOS
# =======================
def cargar_rmsd(ruta: Path, peptido: str, n_replicas: int, dt: float) -> List[pd.DataFrame]:
    """Carga los archivos RMSD (usando el formato 01_RMSD_m_...)."""
    print(f"\n--- Cargando RMSD para {peptido} ---")
    dfs = []
    for i in range(1, n_replicas + 1):
        # Ajusta aquí si tu archivo tiene otro nombre (ej: _p_ en vez de _m_)
        archivo = ruta / f"01_RMSD_p_{peptido}_r{i}.dat"
        if not archivo.exists():
            print(f"  [AVISO] R{i}: No encontrado {archivo}")
            continue
        try:
            df = pd.read_csv(archivo, sep=r"\s+", header=None, names=["Frame", "RMSD"])
            df["Time_ns"] = df["Frame"] * dt
            df["Replica"] = f"R{i}"
            dfs.append(df)
            print(f"  R{i}: Cargado ({len(df)} frames).")
        except Exception as e:
            print(f"  [ERROR] R{i}: {e}")
    return dfs

def cargar_contactos(ruta: Path, peptido: str, n_replicas: int, dt: float) -> List[pd.DataFrame]:
    """Carga los archivos de contacto (.dat)."""
    print(f"\n--- Cargando Contactos para {peptido} ---")
    dfs = []
    for i in range(1, n_replicas + 1):
        archivo = ruta / f"03_contact_{peptido}_r{i}.dat"
        if not archivo.exists():
            print(f"  [AVISO] R{i}: No encontrado {archivo}")
            continue
        try:
            df = pd.read_csv(archivo, sep=r"\s+")
            if "Time_ns" not in df.columns:
                df["Time_ns"] = df["Frame"] * dt
            df["Replica"] = f"R{i}"
            dfs.append(df)
            print(f"  R{i}: Cargado.")
        except Exception as e:
            print(f"  [ERROR] R{i}: {e}")
    return dfs

# =======================
# GRÁFICOS RMSD
# =======================
def graficar_rmsd_comparativo(dfs: List[pd.DataFrame], carpeta: Path, peptido: str):
    """Gráfico 1: Todas las réplicas en una sola imagen."""
    plt.figure(figsize=(10, 6))
    
    colores = sns.color_palette("deep", len(dfs))
    
    for i, df in enumerate(dfs):
        rep = df["Replica"].iloc[0]
        plt.plot(df["Time_ns"], df["RMSD"], label=rep, color=colores[i], alpha=0.8)

    plt.ylim(0, 5) # Escala solicitada
    plt.xlabel("Tiempo (ns)")
    plt.ylabel(r"RMSD ($\AA$)")
    plt.title(f"Comparación RMSD entre Réplicas - {peptido}")
    plt.legend(loc="upper right", frameon=True, title="Réplica")
    plt.tight_layout()
    
    salida = carpeta / f"01_RMSD_Comparativo_{peptido}.png"
    plt.savefig(salida)
    plt.close()
    print(f"Generado: {salida.name}")

def graficar_rmsd_promedio_std(dfs: List[pd.DataFrame], carpeta: Path, peptido: str):
    """Gráfico 2: Promedio y Desviación Estándar."""
    df_concat = pd.concat(dfs)
    # Agrupar por tiempo para sacar media y std
    stats = df_concat.groupby("Time_ns")["RMSD"].agg(["mean", "std"]).reset_index()

    plt.figure(figsize=(10, 6))
    
    # Sombreado STD
    plt.fill_between(
        stats["Time_ns"],
        stats["mean"] - stats["std"],
        stats["mean"] + stats["std"],
        color="navy", alpha=0.2,
        label=r"Desviación Estándar ($\pm 1\sigma$)"
    )
    
    # Línea Promedio
    plt.plot(stats["Time_ns"], stats["mean"], color="navy", linewidth=2, label="Promedio Global")

    plt.ylim(0, 5) # Escala solicitada
    plt.xlabel("Tiempo (ns)")
    plt.ylabel(r"RMSD Promedio ($\AA$)")
    plt.title(f"Estabilidad Estructural Promedio - {peptido}")
    plt.legend(loc="upper right", frameon=True)
    plt.tight_layout()
    
    salida = carpeta / f"01_RMSD_Promedio_{peptido}.png"
    plt.savefig(salida)
    plt.close()
    print(f"Generado: {salida.name}")

# =======================
# GRÁFICO CONTACTOS
# =======================
def graficar_contactos_stacked(dfs: List[pd.DataFrame], carpeta: Path, peptido: str, tiempo_min: float):
    """
    Gráfico de barras apiladas corregido.
    - Ordena residuos naturalmente (CYS2 antes de CYS10).
    - Mantiene la diferenciación visual por réplica.
    - La altura total es la suma de frecuencias (0 a 3).
    """
    # 1. Filtrar datos (Asegúrate que tiempo_min sea consistente con el otro script, ej: 20.0)
    dfs_filtrados = [df[df["Time_ns"] >= tiempo_min] for df in dfs if not df.empty]
    if not dfs_filtrados:
        print("No hay datos de contacto suficientes.")
        return

    # 2. Recolectar frecuencias
    data_list = []
    
    # Identificar columnas de residuos dinámicamente
    # Tomamos el primer DF para ver las columnas, excluyendo metadatos
    cols_meta = ["Frame", "Time_ns", "Replica", "Time"]
    cols_residuos = [c for c in dfs_filtrados[0].columns if c not in cols_meta]

    # Calcular frecuencia promedio por réplica
    for df in dfs_filtrados:
        rep = df["Replica"].iloc[0] # Ej: "R1"
        medias = df[cols_residuos].mean() # Promedio (Frecuencia 0-1)
        for res, frec in medias.items():
            data_list.append({"Residuo": res, "Replica": rep, "Frecuencia": frec})

    df_plot = pd.DataFrame(data_list)

    # 3. Pivotar tabla: Filas=Residuos, Col=Réplicas
    df_pivot = df_plot.pivot(index="Residuo", columns="Replica", values="Frecuencia").fillna(0)

    # --- CORRECCIÓN DE ORDENAMIENTO (NATURAL SORT) ---
    def natural_key(text):
        return [int(c) if c.isdigit() else c.lower() for c in re.split(r'(\d+)', text)]
    
    # Reordenar el índice del pivot para que el eje X salga bien (CYS2, CYS3... CYS10)
    residuos_ordenados = sorted(df_pivot.index, key=natural_key)
    df_pivot = df_pivot.reindex(residuos_ordenados)

    # 4. Graficar
    plt.figure(figsize=(12, 6))
    ax = plt.gca()

    residuos = df_pivot.index
    replicas = df_pivot.columns # ["R1", "R2", "R3"]
    
    # Inicializar el "piso" de las barras en 0
    bottom_vals = pd.Series([0.0] * len(residuos), index=residuos)
    
    # Colores y patrones para diferenciar réplicas
    colores = ["#4c72b0", "#55a868", "#c44e52"] 
    patrones = ["", "///", ".."] 
    
    for i, rep in enumerate(replicas):
        valores = df_pivot[rep]
        
        ax.bar(
            residuos, 
            valores, 
            bottom=bottom_vals,  # Apilar sobre el valor anterior
            label=f"{rep}",
            color=colores[i % len(colores)],
            edgecolor="black",
            linewidth=0.8,
            hatch=patrones[i % len(patrones)],
            alpha=0.85
        )
        
        # Subir el piso para la siguiente réplica
        bottom_vals += valores

    # 5. Ajustes finales
    plt.ylim(0, 3.1) # Escala hasta 3 (suma de 3 réplicas)
    plt.yticks([0, 1, 2, 3])
    plt.ylabel("Frecuencia Acumulada (Suma de Réplicas)", fontweight='bold')
    plt.xlabel("Residuo", fontweight='bold')
    plt.title(f"Interacción Detallada por Réplica - {peptido}")
    
    plt.legend(loc="upper right", title="Réplica", frameon=True)
    plt.xticks(rotation=45)
    plt.grid(axis='y', linestyle='--', alpha=0.5)
    plt.tight_layout()

    salida = carpeta / f"03_Contacto_Acumulado_{peptido}.png"
    plt.savefig(salida)
    plt.close()
    print(f"Generado: {salida.name}")

# =======================
# MAIN
# =======================
def main():
    configurar_estilo()
    CARPETA_SALIDA.mkdir(exist_ok=True, parents=True)

    # Cargar
    rmsd_data = cargar_rmsd(RUTA_DATOS, PEPTIDO, REPLICAS, DT)
    contact_data = cargar_contactos(RUTA_DATOS, PEPTIDO, REPLICAS, DT)

    # Graficar RMSD
    if rmsd_data:
        graficar_rmsd_comparativo(rmsd_data, CARPETA_SALIDA, PEPTIDO)
        graficar_rmsd_promedio_std(rmsd_data, CARPETA_SALIDA, PEPTIDO)
    else:
        print("[AVISO] No hay datos RMSD para graficar.")

    # Graficar Contactos
    if contact_data:
        graficar_contactos_stacked(contact_data, CARPETA_SALIDA, PEPTIDO, TIEMPO_INICIO_CONTACTOS)
    else:
        print("[AVISO] No hay datos de Contactos para graficar.")

    print("\nProceso finalizado.")

if __name__ == "__main__":
    main()