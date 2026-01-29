import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import warnings
import numpy as np
import re
import sys
from matplotlib.lines import Line2D

# =======================
# CONFIGURACIÓN GENERAL
# =======================
REPLICAS = 3
DT = 0.04                   # ns por frame
TIEMPO_TOTAL = 70.0         # ns total simulación
TIEMPO_EQUILIBRIO = 20.0    # ns a descartar como equilibrado
VENTANA_SUAVIZADO = 20
CARPETA_SALIDA = Path("COMPARISON_GRAPHICS")

# =======================
# UTILIDADES
# =======================
def configurar_estilo():
    warnings.filterwarnings("ignore")
    try:
        sns.set_theme(style="whitegrid", context="paper", font_scale=1.2)
        sns.set_palette("bright")
    except:
        plt.style.use("seaborn-whitegrid")
    
    plt.rcParams.update({
        "figure.dpi": 300,
        "lines.linewidth": 2,
        "axes.grid": True,
        "grid.alpha": 0.3,
        "legend.frameon": True,
        "legend.fancybox": True,
        "axes.titlesize": 14,
        "axes.labelsize": 12,
    })

def natural_key(text):
    """Para ordenar residuos correctamente."""
    return [int(c) if c.isdigit() else c.lower() for c in re.split(r'(\d+)', text)]

def detectar_sistemas() -> list:
    """Busca carpetas que empiecen con un número."""
    root = Path(".")
    sistemas = [p.name for p in root.iterdir() if p.is_dir() and p.name[0].isdigit() and "_" in p.name]
    sistemas.sort(key=lambda x: int(x.split('_')[0]))
    return sistemas

# =======================
# CARGA DE DATOS
# =======================

def cargar_rmsd_p_directo(sistemas: list) -> pd.DataFrame:
    print("--- Cargando RMSD (Archivos _p_) ---")
    dfs = []
    for sistema in sistemas:
        for r in range(1, REPLICAS + 1):
            archivo = Path(f"{sistema}/analysis/data/01_RMSD_p_{sistema}_r{r}.dat")
            if not archivo.exists(): continue
            try:
                df = pd.read_csv(archivo, sep=r"\s+", header=None, names=["Frame", "RMSD"])
                df["Time_ns"] = df["Frame"] * DT
                min_time = df["Time_ns"].min()
                df["Time_Prod"] = df["Time_ns"] - min_time
                df["Sistema"] = sistema  # Cambiado a "Sistema"
                df["Replica"] = f"R{r}"
                dfs.append(df)
            except Exception as e:
                print(f"  [Error] {archivo}: {e}")
    return pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame()

def cargar_contactos_acumulados(sistemas: list) -> pd.DataFrame:
    print("--- Cargando Contactos (Acumulado) ---")
    data_final = []

    for sistema in sistemas:
        ruta_datos = Path(f"{sistema}/analysis/data")
        dfs_replicas = []
        
        for r in range(1, REPLICAS + 1):
            archivo = ruta_datos / f"03_contact_{sistema}_r{r}.dat"
            if not archivo.exists(): continue
            try:
                df = pd.read_csv(archivo, sep=r"\s+")
                if "Time_ns" not in df.columns:
                    df["Time_ns"] = df["Frame"] * DT
                
                # Filtrar Equilibrio (Usamos TIEMPO_EQUILIBRIO - 20 si es necesario ajustar, 
                # o TIEMPO_EQUILIBRIO directo según tu lógica previa)
                df = df[df["Time_ns"] >= TIEMPO_EQUILIBRIO - 15] 
                if not df.empty:
                    dfs_replicas.append(df)
            except: pass
        
        if not dfs_replicas: continue

        frecuencias_sistema = {}
        cols_meta = ["Frame", "Time_ns", "Replica", "Time"]
        cols_residuos = [c for c in dfs_replicas[0].columns if c not in cols_meta]

        for df in dfs_replicas:
            promedios = df[cols_residuos].mean()
            for res, valor in promedios.items():
                frecuencias_sistema[res] = frecuencias_sistema.get(res, 0.0) + valor

        for res, val in frecuencias_sistema.items():
            data_final.append({
                "Residuo": res,
                "Frecuencia_Acumulada": val,
                "Sistema": sistema
            })

    return pd.DataFrame(data_final)

def cargar_distancias(sistemas: list) -> pd.DataFrame:
    print("--- Cargando Distancias ---")
    dfs = []
    for sistema in sistemas:
        for r in range(1, REPLICAS + 1):
            archivo = Path(f"{sistema}/analysis/data/05_distancias_{sistema}_r{r}.dat")
            if not archivo.exists(): continue
            try:
                df = pd.read_csv(archivo, sep=r"\s+")
                if "Tiempo" in df.columns: df.rename(columns={"Tiempo": "Time_ns"}, inplace=True)
                
                cols_residuos = [c for c in df.columns if c not in ["Time_ns", "Replica", "Distancia_Global", "Distancia_Global_Smooth"]]
                
                if "Distancia_Global" in df.columns:
                    # Suavizado Global
                    df["Dist_Global_Smooth"] = df["Distancia_Global"].rolling(VENTANA_SUAVIZADO, center=True, min_periods=1).mean()
                    
                    # Suavizado Residuos
                    for col in cols_residuos:
                        if not col.endswith("_Smooth"):
                            df[f"{col}_Smooth"] = df[col].rolling(VENTANA_SUAVIZADO, center=True, min_periods=1).mean()
                    
                    df["Sistema"] = sistema # Cambiado a "Sistema"
                    df["Replica"] = f"R{r}"
                    dfs.append(df)
            except Exception as e:
                print(f"  [Error] {archivo}: {e}")
    return pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame()

# =======================
# FUNCIONES DE GRAFICADO
# =======================

def plot_rmsd_produccion(df_prod):
    if df_prod.empty: return
    plt.figure(figsize=(10, 6))
    sns.lineplot(
        data=df_prod, x="Time_Prod", y="RMSD", 
        hue="Sistema", linewidth=2, errorbar="sd"
    )
    plt.xlim(0, 50)
    plt.ylim(0, 5)
    plt.xlabel("Tiempo de Producción (ns)")
    plt.ylabel(r"RMSD ($\AA$)")
    plt.title("Comparación Estabilidad Estructural (RMSD) - Producción")
    plt.legend(loc="upper right", title="Sistema")
    plt.tight_layout()
    salida = CARPETA_SALIDA / "01_RMSD_Comparativo_Produccion.png"
    plt.savefig(salida)
    plt.close()
    print(f"  -> Generado: {salida.name}")

def plot_contactos_barras(df_contactos):
    if df_contactos.empty: return
    
    residuos_unicos = sorted(df_contactos["Residuo"].unique(), key=natural_key)
    
    plt.figure(figsize=(14, 7))
    ax = sns.barplot(
        data=df_contactos,
        x="Residuo",
        y="Frecuencia_Acumulada",
        hue="Sistema",
        order=residuos_unicos,
        edgecolor="black",
        linewidth=0.8,
        alpha=0.9
    )
    
    plt.ylim(0, 3.1)
    plt.ylabel("Frecuencia Acumulada (Suma 3 Réplicas)", fontweight='bold')
    plt.xlabel("Residuo", fontweight='bold')
    plt.title("Comparación de Interacciones por Residuo entre Sistemas", fontweight='bold', pad=40)
    
    # --- TRUCO: "Sistemas" en línea ---
    handles, labels = ax.get_legend_handles_labels()
    dummy_handle = Line2D([], [], color='none', label='Sistemas')
    handles = [dummy_handle] + handles
    labels = ["Sistemas"] + labels
    
    plt.legend(
        handles=handles,
        labels=labels,
        loc='lower center',
        bbox_to_anchor=(0.5, 1.01),
        ncol=len(labels),
        frameon=True,
        shadow=True,
        fancybox=True,
        handletextpad=0.5
    )

    plt.xticks(ha='right')
    plt.grid(axis='y', linestyle='--', alpha=0.5)
    plt.tight_layout()
    
    salida = CARPETA_SALIDA / "03_Contactos_Barras_Suma_Comparativo.png"
    plt.savefig(salida)
    plt.close()
    print(f"  -> Generado: {salida.name}")

def plot_distancia_1_m_custom(df):
    """
    Gráfico Full: -15 a 25.
    Sombreado Dorado entre -10 y 0 con etiqueta de texto.
    Sombreado Azul en el resto (-15 a -10 y 0 a 25).
    """
    if df.empty: return
    plt.figure(figsize=(10, 6))
    
    # --- SOMBREADOS ---
    # Zona Azul Inferior
    plt.axhspan(-15, -10, color='lightsteelblue', alpha=0.2, zorder=0, linewidth=0)
    
    # Zona Dorada Central
    plt.axhspan(-10, 0, color='gold', alpha=0.25, zorder=0, linewidth=0)
    
    # Zona Azul Superior
    plt.axhspan(0, 25, color='lightsteelblue', alpha=0.2, zorder=0, linewidth=0)
    
    # --- ETIQUETA SUPERFICIE (NUEVO) ---
    # Ubicada al centro del tiempo (TIEMPO_TOTAL/2) y al centro de la banda (-5)
    plt.text(TIEMPO_TOTAL / 2, -5, "Superficie de Oro (Au111)", 
             ha='center', va='center', 
             color='#665500', # Un tono dorado oscuro/marrón para contraste y estética
             fontweight='bold', 
             fontsize=12,
             alpha=0.9)

    # --- LÍNEAS DE DATOS ---
    sns.lineplot(data=df, x="Time_ns", y="Dist_Global_Smooth", hue="Sistema", linewidth=2, errorbar=None)
    
    # --- REFERENCIA EQUILIBRIO ---
    plt.axvline(x=TIEMPO_EQUILIBRIO, color='black', linestyle='--', alpha=0.5)
    
    y_text = 22
    plt.text(TIEMPO_EQUILIBRIO - 0.5, y_text, "Equilibrado", ha='right', fontweight='bold', color='darkred')
    plt.text(TIEMPO_EQUILIBRIO + 0.5, y_text, "Producción", ha='left', fontweight='bold', color='darkgreen')
    
    # --- AJUSTES FINALES ---
    plt.ylim(-15, 25)
    plt.xlim(0, TIEMPO_TOTAL)
    plt.xlabel("Tiempo (ns)")
    plt.ylabel(r"Distancia Mínima Global ($\AA$)")
    plt.title("Comparación Distancia Global (Completa)")
    plt.legend(loc="upper right", title="Sistema")
    plt.tight_layout()
    
    salida = CARPETA_SALIDA / "05_Distancia_Full_Comparativo.png"
    plt.savefig(salida)
    plt.close()
    print(f"  -> Generado: {salida.name}")

def plot_distancia_promedio_std_final(df_prod):
    """
    Gráfico con sombra de desviación estándar, pero titulado Promedio.
    Sustituye al anterior 'Promedio' plano y al 'Std'.
    """
    if df_prod.empty: return
    plt.figure(figsize=(10, 6))
    
    # Usamos errorbar='sd' para mostrar la desviación
    sns.lineplot(
        data=df_prod, x="Time_Prod", y="Dist_Global_Smooth", 
        hue="Sistema", linewidth=2, errorbar="sd", alpha=0.9
    )
    
    plt.xlim(0, 50)
    plt.ylim(0, 5) # Ajustar si es necesario
    plt.xlabel("Tiempo de Producción (ns)")
    plt.ylabel(r"Distancia Mínima Global ($\AA$)")
    # Título sin "+- DE"
    plt.title("Comparación Distancia Promedio (Producción)")
    plt.legend(loc="upper right", title="Sistema")
    plt.tight_layout()
    
    # Guardamos con el nombre del Promedio para que sea el principal
    salida = CARPETA_SALIDA / "05_Distancia_Promedio_Sistemas.png"
    plt.savefig(salida)
    plt.close()
    print(f"  -> Generado: {salida.name}")

def plot_distancia_por_residuo_individual(df_prod):
    """
    Genera UN gráfico POR CADA RESIDUO comparando los sistemas.
    """
    # Detectar columnas de residuos suavizados (ej: CYS2_Smooth)
    cols_residuos = [c for c in df_prod.columns if ("CYS" in c or "TRP" in c) and "_Smooth" in c]
    if not cols_residuos: return
    
    print("  -> Generando gráficos individuales por residuo...")
    
    for col in cols_residuos:
        nombre_residuo = col.replace("_Smooth", "")
        
        plt.figure(figsize=(8, 5))
        sns.lineplot(
            data=df_prod, 
            x="Time_Prod", 
            y=col, 
            hue="Sistema", 
            linewidth=2, 
            errorbar=None
        )
        
        plt.xlim(0, 50)
        # Ajustamos eje Y dinámicamente o fijo según preferencia. Fijo suele ser mejor para comparar.
        plt.ylim(0, 15) 
        
        plt.xlabel("Tiempo de Producción (ns)")
        plt.ylabel(r"Distancia Mínima ($\AA$)")
        plt.title(f"Distancia Residuo: {nombre_residuo}")
        plt.legend(loc="upper right", title="Sistema")
        plt.tight_layout()
        
        salida = CARPETA_SALIDA / f"05_Distancia_Residuo_{nombre_residuo}.png"
        plt.savefig(salida)
        plt.close()

def plot_distancia_por_residuo_general(df_prod):
    # Detectar columnas de residuos suavizados (ej: CYS2_Smooth)
    cols_residuos = [c for c in df_prod.columns if ("CYS" in c or "TRP" in c) and "_Smooth" in c]
    if not cols_residuos: return
    
    df_long = df_prod.melt(id_vars=["Time_Prod", "Sistema", "Replica"], value_vars=cols_residuos, var_name="Residuo", value_name="Distancia")
    df_long["Residuo"] = df_long["Residuo"].str.replace("_Smooth", "")
    
    plt.figure(figsize=(12, 6))
    sns.lineplot(data=df_long, x="Time_Prod", y="Distancia", hue="Residuo", style="Sistema", linewidth=2, errorbar=None)
    plt.xlim(0, 50)
    plt.ylim(0, 15)
    plt.xlabel("Tiempo de Producción (ns)")
    plt.ylabel(r"Distancia Mínima Residuo ($\AA$)")
    plt.title("Comparación Detallada de Residuos entre Sistemas")
    plt.legend(loc="upper right", bbox_to_anchor=(1.15, 1))
    plt.tight_layout()
    salida = CARPETA_SALIDA / "05_Distancia_Residuo_Comparativo.png"
    plt.savefig(salida)
    plt.close()
    print(f"  -> Generado: {salida.name}")

# =======================
# MAIN
# =======================
def main():
    configurar_estilo()
    CARPETA_SALIDA.mkdir(exist_ok=True, parents=True)
    
    # 1. Detectar sistemas
    sistemas = detectar_sistemas()
    if not sistemas:
        print("Error: No se detectaron carpetas de sistemas (ej: 1_CLK).")
        return
    print(f"Sistemas detectados: {sistemas}\n")

    # 2. Procesar RMSD
    df_rmsd = cargar_rmsd_p_directo(sistemas)
    if not df_rmsd.empty:
        plot_rmsd_produccion(df_rmsd)
    else:
        print("[AVISO] No se encontraron datos de RMSD.")

    # 3. Procesar Contactos
    df_contactos = cargar_contactos_acumulados(sistemas)
    if not df_contactos.empty:
        plot_contactos_barras(df_contactos)
    else:
        print("[AVISO] No se encontraron datos de Contactos.")

    # 4. Procesar Distancias
    df_dist = cargar_distancias(sistemas)
    if not df_dist.empty:
        # Gráfico Full personalizado
        plot_distancia_1_m_custom(df_dist)
        
        # Filtrar solo producción
        df_prod_dist = df_dist[df_dist["Time_ns"] >= TIEMPO_EQUILIBRIO].copy()
        df_prod_dist["Time_Prod"] = df_prod_dist["Time_ns"] - TIEMPO_EQUILIBRIO
        
        # Gráfico Promedio (con STD visual pero nombre limpio)
        plot_distancia_promedio_std_final(df_prod_dist)
        
        # Gráficos individuales por residuo
        plot_distancia_por_residuo_individual(df_prod_dist)

        # Gráficos Generales por residuo
        plot_distancia_por_residuo_general(df_prod_dist)
    else:
        print("[AVISO] No se encontraron datos de Distancias.")

    print("\n¡Análisis comparativo completado!")

if __name__ == "__main__":
    main()