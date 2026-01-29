from pathlib import Path
from typing import List, Dict, Optional
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys

# =======================
# CONFIGURACIÓN DE USUARIO
# =======================
# Obtener péptido desde argumento de línea de comandos
if len(sys.argv) > 1:
    PEPTIDO = sys.argv[1] # Toma el valor enviado por bash
else:
    print("Error: Debes indicar el péptido como argumento.")
    sys.exit(1)

RUTA_DATOS = Path("data")
CARPETA_SALIDA = Path("graphics")
DT = 0.04  # ns por frame
REPLICAS = 3
VENTANA_SUAVIZADO = 20
TIEMPO = 5.0

# =======================
# ESTILO VISUAL
# =======================
def configurar_estilo() -> None:
    """Configura el estilo global de las figuras."""
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
    })

# =======================
# CARGA RMSD
# =======================
def cargar_datos_rmsd_1(ruta: Path, peptido: str, n_replicas: int, dt: float) -> List[pd.DataFrame]:
    """Carga archivos RMSD de todas las réplicas."""
    print(f"\n--- Cargando RMSD para {peptido} ---")
    dfs = []

    for i in range(1, n_replicas + 1):
        archivo = ruta / f"01_RMSD_m_{peptido}_r{i}.dat"
        if archivo is None:
            print(f"  [AVISO] R{i}: RMSD no encontrado ({archivo}).")
            continue

        try:
            df = pd.read_csv(archivo, sep=r"\s+", header=None, names=["Frame", "RMSD"])
            df["Time_ns"] = df["Frame"] * dt
            df["Replica"] = f"R{i}"
            dfs.append(df)
            print(f"  R{i}: OK ({len(df)} frames).")
        except Exception as e:
            print(f"  [ERROR] R{i}: Fallo al leer {archivo.name}: {e}")

    return dfs

def cargar_datos_rmsd_2(ruta: Path, peptido: str, n_replicas: int, dt: float) -> List[pd.DataFrame]:
    """Carga archivos RMSD de todas las réplicas."""
    print(f"\n--- Cargando RMSD para {peptido} ---")
    dfs = []

    for i in range(1, n_replicas + 1):
        archivo = ruta / f"01_RMSD_p_{peptido}_r{i}.dat"
        if archivo is None:
            print(f"  [AVISO] R{i}: RMSD no encontrado ({archivo}).")
            continue

        try:
            df = pd.read_csv(archivo, sep=r"\s+", header=None, names=["Frame", "RMSD"])
            df["Time_ns"] = df["Frame"] * dt
            df["Replica"] = f"R{i}"
            dfs.append(df)
            print(f"  R{i}: OK ({len(df)} frames).")
        except Exception as e:
            print(f"  [ERROR] R{i}: Fallo al leer {archivo.name}: {e}")

    return dfs

# =======================
# CARGA RMSF
# =======================
def cargar_datos_rmsf(ruta: Path, peptido: str, n_replicas: int) -> List[Dict[str, Optional[pd.DataFrame]]]:
    """
    Carga RMSF desde el nuevo formato único (Resname, Resid, RMSF).
    Filtra automáticamente CYS y TRP del dataframe principal.
    """
    print(f"\n--- Cargando RMSF para {peptido} ---")
    datos = []

    for i in range(1, n_replicas + 1):
        # Nombre del archivo generado por el nuevo Tcl
        archivo = ruta / f"02_RMSF_{peptido}_r{i}.dat"
        
        # Estructura para guardar los datos
        info = {"Replica": f"R{i}", "ALL": None, "CYS": None, "TRP": None}

        if not archivo.exists():
            print(f"  [AVISO] R{i}: Archivo RMSF no encontrado ({archivo.name}).")
            continue

        try:
            # Leer archivo ignorando comentarios (#)
            # Formato esperado: Resname Resid RMSF
            df = pd.read_csv(archivo, sep=r"\s+", comment="#", header=None, names=["Resname", "Resid", "RMSF"])
            
            info["ALL"] = df

            # Filtrar subgrupos de interés directamente del DF principal
            df_cys = df[df["Resname"] == "CYS"]
            if not df_cys.empty:
                info["CYS"] = df_cys
            
            df_trp = df[df["Resname"] == "TRP"]
            if not df_trp.empty:
                info["TRP"] = df_trp
            
            print(f"  R{i}: OK ({len(df)} residuos).")
            datos.append(info)

        except Exception as e:
            print(f"  [ERROR] R{i}: Error leyendo RMSF: {e}")

    return datos

# =======================
# CARGA CONTACTOS
# =======================
def cargar_datos_contacto(ruta: Path, peptido: str, n_replicas: int, dt: float) -> List[pd.DataFrame]:
    """Carga los datos de contacto reales (.dat) de todas las réplicas."""
    print(f"\n--- Cargando Contactos para {peptido} ---")
    dfs = []

    for i in range(1, n_replicas + 1):
        # Busca archivos con nombres comunes generados por los scripts previos
        base = ruta / f"03_contact_{peptido}_r{i}.dat"
        
        if base.exists():
            try:
                # Lee el base detectando espacios o tabs
                df = pd.read_csv(base, sep=r"\s+")
                
                # Agrega columna de tiempo si no existe
                if "Time_ns" not in df.columns:
                    df["Time_ns"] = df["Frame"] * dt
                
                df["Replica"] = f"R{i}"
                dfs.append(df)
                print(f"  R{i}: OK ({base.name}).")
                
            except Exception as e:
                print(f"  [ERROR] R{i}: No se pudo leer {base.name}. Error: {e}")
        else:
            print(f"  [AVISO] R{i}: Archivo de contactos no encontrado en {ruta}")

    return dfs

# =======================
# CARGA ENERGÍA DESDE LOG
# =======================
def cargar_energia_desde_log(ruta_dummy: Path, peptido: str, n_replicas: int, dt: float) -> List[pd.DataFrame]:
    """
    Lee los archivos LOG de NAMD desde la carpeta relative '../namd/'.
    """
    print(f"\n--- Cargando Energía desde LOGs (../namd/) para {peptido} ---")
    dfs = []
    
    # Definimos la ruta específica solicitada
    ruta_namd = Path("../namd")

    if not ruta_namd.exists():
        print(f"  [ERROR] La carpeta '{ruta_namd.resolve()}' no existe.")
        return []

    for i in range(1, n_replicas + 1):
        # Construcción del nombre exacto solicitado:
        nombre_archivo = f"03_prod_pamam_{peptido}_r{i}.1.log"
        archivo = ruta_namd / nombre_archivo

        if not archivo.exists():
            print(f"  [AVISO] R{i}: Archivo no encontrado: {archivo}")
            continue

        try:
            data_rows = []
            cols = []
            
            with open(archivo, 'r') as f:
                for line in f:
                    # Detectar cabecera
                    if line.startswith("ETITLE:"):
                        cols = line.split()[1:] 
                    
                    # Detectar datos de energía
                    elif line.startswith("ENERGY:"):
                        parts = line.split()
                        # Validar que la línea tenga datos y coincida con las columnas
                        if len(parts) > 1 and len(cols) > 0:
                            # parts[0] es la etiqueta "ENERGY:", los datos empiezan en [1:]
                            vals = parts[1:]
                            try:
                                # Convertir strings a float
                                vals_float = [float(v) for v in vals]
                                # Solo agregar si coincide la longitud con las columnas
                                if len(vals_float) == len(cols):
                                    data_rows.append(vals_float)
                            except ValueError:
                                continue

            if not data_rows:
                print(f"  [ERROR] R{i}: El log existe pero no tiene líneas 'ENERGY:' válidas.")
                continue

            # Crear DataFrame
            df = pd.DataFrame(data_rows, columns=cols)
            
            # --- AJUSTE DE TIEMPO ---
            # Si el log tiene pasos (TS), usamos eso, o generamos basado en frames.
            # Aquí asumimos que cada línea de energía corresponde a un frame guardado.
            df["Frame"] = range(len(df))
            df["Time_ns"] = df["Frame"] * dt
            df["Replica"] = f"R{i}"
            
            # Extraer solo Energía Potencial
            if "POTENTIAL" in df.columns:
                df_clean = df[["Frame", "Time_ns", "POTENTIAL", "Replica"]].copy()
                df_clean.rename(columns={"POTENTIAL": "Potential_Energy"}, inplace=True)
                dfs.append(df_clean)
                print(f"  R{i}: OK ({len(df)} registros). E. Potencial Media: {df_clean['Potential_Energy'].mean():.2f}")
            else:
                print(f"  [ERROR] R{i}: Columna POTENTIAL no encontrada en el log.")

        except Exception as e:
            print(f"  [ERROR] R{i}: Fallo leyendo {archivo.name}: {e}")

    return dfs

# =======================
# GRÁFICOS RMSD
# =======================
def graficar_rmsd_individuales(dfs, carpeta, peptido):
    for df in dfs:
        rep = df["Replica"].iloc[0]
        max_time = df["Time_ns"].max()
        
        plt.figure(figsize=(12, 6))
        
        # --- ZONAS DE FONDO ---
        # Equilibrado (0 a 20 ns): Rojo tenue
        plt.axvspan(0, 20, color='tab:red', alpha=0.1, label='Equilibrado (0-20ns)', zorder=0)
        
        # Producción (20 ns en adelante): Verde tenue
        if max_time > 20:
            plt.axvspan(20, max_time, color='tab:green', alpha=0.1, label='Producción (>20ns)', zorder=0)

        # Gráfico de línea
        plt.plot(df["Time_ns"], df["RMSD"], label=rep, color="teal", zorder=2)
        
        plt.xlabel("Tiempo (ns)")
        plt.ylabel(r"RMSD ($\AA$)")
        plt.title(f"Desviación estructural – {peptido} {rep}")        
        plt.legend(loc='lower right')
        plt.savefig(carpeta / f"01_RMSD_{peptido}_{rep}.png")
        plt.close()

def graficar_rmsd_promedio(dfs, carpeta, peptido, bool):
    df = pd.concat(dfs)
    stats = df.groupby("Time_ns")["RMSD"].agg(["mean", "std"]).reset_index()
    max_time = stats["Time_ns"].max()

    plt.figure(figsize=(12, 6))
    if bool:
        # --- ZONAS DE FONDO ---
        plt.axvspan(0, 20, color='tab:red', alpha=0.1, label='Equilibrado', zorder=0)
        if max_time > 20:
            plt.axvspan(20, max_time, color='tab:green', alpha=0.1, label='Producción', zorder=0)

    # Línea promedio
    plt.plot(stats["Time_ns"], stats["mean"], label="Promedio", color="darkblue", linewidth=2, zorder=3)
    
    # Sombreado de Desviación Estándar
    plt.fill_between(
        stats["Time_ns"],
        stats["mean"] - stats["std"],
        stats["mean"] + stats["std"],
        color="blue", alpha=0.15,
        label=r"Desviación Estándar ($\pm 1\sigma$)",
        zorder=2
    )

    plt.xlabel("Tiempo (ns)")
    plt.ylabel(r"RMSD Promedio ($\AA$)")

    if bool:
        plt.title(f"Desviación estructural Global – {peptido}")
        plt.legend(loc='lower right')
        plt.savefig(carpeta / f"01_RMSD_m_Promedio_{peptido}.png")
        plt.close()
    else:
        plt.title(f"Desviación estructural Global en Produccion– {peptido}")
        plt.legend(loc='best')
        plt.savefig(carpeta / f"01_RMSD_p_Promedio_{peptido}.png")
    plt.close()

# =======================
# GRÁFICOS RMSF
# =======================
def graficar_rmsf_individuales(datos, carpeta, peptido):
    for d in datos:
        if d["ALL"] is None: continue
        
        df = d["ALL"].copy()
        df["Label"] = df["Resid"].astype(str) + " - " + df["Resname"]
        
        plt.figure(figsize=(12, 6))
        
        sns.scatterplot(
            data=df, 
            x="Resid", 
            y="RMSF", 
            hue="Resname", 
            palette="tab10", 
            s=100, 
            edgecolor="black",
            zorder=3,
            legend=False  # <--- SIN LEYENDA
        )
        
        plt.plot(df["Resid"], df["RMSF"], color="gray", alpha=1, linewidth=1, zorder=1)

        plt.xlabel("Residuo")
        plt.ylabel(r"RMSF ($\AA$)")
        plt.title(f"Fluctuación estructural - {peptido} {d['Replica']}")
        
        plt.xticks(
            ticks=df["Resid"], 
            labels=df["Label"], 
            rotation=45,
            fontsize=10
        )
        
        plt.grid(True, linestyle='--', alpha=0.5)
        plt.tight_layout()
        plt.savefig(carpeta / f"02_RMSF_{peptido}_{d['Replica']}.png")
        plt.close()

def graficar_rmsf_promedio(datos, carpeta, peptido):
    dfs = [d["ALL"] for d in datos if d["ALL"] is not None]
    if not dfs: return
    df_concat = pd.concat(dfs)
    stats = df_concat.groupby(["Resid", "Resname"])["RMSF"].agg(["mean", "std"]).reset_index()
    stats["Label"] = stats["Resid"].astype(str) + " - " + stats["Resname"]
    plt.figure(figsize=(12, 6))
    # Área de sombra
    plt.fill_between(
        stats["Resid"],
        stats["mean"] - stats["std"],
        stats["mean"] + stats["std"],
        color="gray", alpha=0.15
        # Sin label para que no genere necesidad de leyenda
    )
    plt.plot(stats["Resid"], stats["mean"], color="gray", alpha=0.4, linewidth=1.5, zorder=1)
    sns.scatterplot(
        data=stats,
        x="Resid",
        y="mean",
        hue="Resname",
        palette="tab10",
        s=120,
        edgecolor="black",
        zorder=5,
        legend=False # <--- SIN LEYENDA
    )

    plt.xlabel("Residuo")
    plt.ylabel(r"RMSF Promedio ($\AA$)")
    plt.title(f"Fluctuación estructural Global – {peptido}")
    
    plt.xticks(
        ticks=stats["Resid"], 
        labels=stats["Label"], 
        rotation=45,
        fontsize=10
    )
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.savefig(carpeta / f"02_RMSF_Promedio_{peptido}.png")
    plt.close()

# =======================
# GRÁFICOS CONTACTOS
# =======================
def graficar_contactos_individuales(dfs: List[pd.DataFrame], carpeta: Path, peptido: str, tiempo_inicio: float):
    """
    Grafica la frecuencia de contacto para cada réplica individualmente.
    Maneja DataFrames vacíos y filtra por tiempo.
    """
    for df in dfs:
        # 1. SEGURIDAD: Verificar si el DF original está vacío antes de intentar leer nada
        if df.empty:
            print("  [AVISO] Se saltó una réplica porque el DataFrame estaba vacío.")
            continue

        rep = df["Replica"].iloc[0]
        
        # 2. FILTRADO POR TIEMPO
        df_filtrado = df[df["Time_ns"] >= tiempo_inicio].copy()
        
        # 3. SEGURIDAD: Verificar si quedó vacío DESPUÉS de filtrar el tiempo
        if df_filtrado.empty:
            print(f"  [AVISO] {rep}: No hay datos después de {tiempo_inicio} ns.")
            continue

        cols_datos = [c for c in df_filtrado.columns if c not in ["Frame", "Time_ns", "Replica"]]
        cols_datos.sort()

        if not cols_datos:
            print(f"  [AVISO] No se encontraron datos de interacción en {rep}")
            continue

        # Transformar a formato largo
        df_long = df_filtrado.melt(
            id_vars=["Time_ns", "Replica"], 
            value_vars=cols_datos,
            var_name="Interacción", 
            value_name="Frecuencia/Contacto"
        )

        plt.figure(figsize=(10, 6)) 
        
        sns.barplot(
            data=df_long, 
            x="Interacción", 
            y="Frecuencia/Contacto", 
            hue="Interacción", 
            palette="tab10", 
            legend=False, 
            errorbar=None
        )

        # Título con corrección de LaTeX (\\geq)
        plt.title(f"Perfil de Interacción - {peptido} {rep} (t $\\geq$ {tiempo_inicio} ns)")
        plt.ylabel("Frecuencia de Contacto (0-1)")
        plt.xlabel("Residuo")
        plt.ylim(0, 1.1)
        plt.xticks(rotation=45) 
        plt.grid(axis='y', linestyle='--', alpha=0.5)
        
        plt.tight_layout()
        plt.savefig(carpeta / f"03_Contact_{peptido}_{rep}.png")
        plt.close()

def graficar_contactos_promedio(dfs: List[pd.DataFrame], carpeta: Path, peptido: str, tiempo_inicio: float):
    """
    Grafica la frecuencia ACUMULADA (suma) de contactos considerando todas las réplicas.
    """
    if not dfs: return

    # Filtramos por tiempo, ignorando DataFrames vacíos
    dfs_filtrados = [df[df["Time_ns"] >= tiempo_inicio] for df in dfs if not df.empty]
    
    # Verificar si quedaron datos útiles
    if not any(not d.empty for d in dfs_filtrados):
        print(f"  [AVISO] No hay datos suficientes después de {tiempo_inicio} ns para el gráfico acumulado.")
        return

    df_concat = pd.concat(dfs_filtrados)

    if df_concat.empty:
        return

    cols_datos = [c for c in df_concat.columns if c not in ["Frame", "Time_ns", "Replica"]]
    cols_datos.sort()

    if not cols_datos: return

    # 1. Transformar a formato largo
    df_long = df_concat.melt(
        id_vars=["Time_ns", "Replica"], 
        value_vars=cols_datos,
        var_name="Interacción", 
        value_name="Contacto" # 1 si hay contacto, 0 si no
    )

    # 2. Calcular la frecuencia promedio POR RÉPLICA (temporal)
    # Esto nos da un valor entre 0 y 1 para cada réplica (ej: R1=0.5, R2=0.8)
    df_frec_replica = df_long.groupby(["Replica", "Interacción"])["Contacto"].mean().reset_index()

    # 3. Calcular la SUMA (Acumulada) de las frecuencias de todas las réplicas
    # Si tenemos 3 réplicas, el valor máximo posible será 3.0
    df_acumulado = df_frec_replica.groupby("Interacción")["Contacto"].sum().reset_index()

    plt.figure(figsize=(12, 6))

    sns.barplot(
        data=df_acumulado, 
        x="Interacción", 
        y="Contacto", 
        hue="Interacción", 
        palette="tab10",
        legend=False,
        edgecolor="black"
    )

    plt.title(f"Frecuencia Acumulada de Interacción - {peptido} (t $\\geq$ {tiempo_inicio} ns)")
    plt.ylabel("Frecuencia Acumulada (Suma de Réplicas)")
    plt.xlabel("Residuos Interactuantes")
    
    # Importante: Quitamos el límite fijo (0, 1.1) porque la suma puede ser mayor a 1
    # Si quieres ver el techo teórico, puedes usar el número de réplicas:
    # plt.ylim(0, len(dfs_filtrados) + 0.1) 
    
    plt.xticks(rotation=45)
    plt.grid(axis='y', linestyle='--', alpha=0.5)
    plt.tight_layout()
    
    # Guardamos con un nombre distinto o el mismo según prefieras
    plt.savefig(carpeta / f"03_Contact_Acumulado_{peptido}.png")
    plt.close()

# =======================
# GRÁFICOS ENERGÍA
# =======================
def graficar_energia_individuales(dfs, carpeta, peptido):
    """
    Grafica la energía de cada réplica.
    Muestra los datos crudos (transparentes) de fondo y una línea suavizada encima.
    """
    for df in dfs:
        rep = df["Replica"].iloc[0]
        
        # Calcular media móvil (suavizado)
        df["Energy_Smooth"] = df["Potential_Energy"].rolling(window=VENTANA_SUAVIZADO, center=True).mean()

        plt.figure(figsize=(10, 6))
        
        # 1. Gráfica de fondo: Datos crudos (ruidosos) muy transparentes
        sns.lineplot(
            data=df, 
            x="Time_ns", 
            y="Potential_Energy", 
            color="blue", 
            alpha=0.2,
            linewidth=0.5,
            legend=False,
            zorder=1
        )

        # 2. Gráfica principal: Datos suavizados (Tendencia)
        sns.lineplot(
            data=df, 
            x="Time_ns", 
            y="Energy_Smooth", 
            color="blue", 
            linewidth=1.5,
            label=f"{rep} (Suavizado)"
        )

        plt.xlabel("Tiempo (ns)")
        plt.ylabel("Energía Potencial (kcal/mol)")
        plt.title(f"Energía Potencial - {peptido} {rep}")
        plt.grid(True, linestyle='--', alpha=0.5)
        plt.legend()
        plt.tight_layout()
        plt.savefig(carpeta / f"04_Energy_{peptido}_{rep}.png")
        plt.close()

def graficar_energia_promedio(dfs, carpeta, peptido):
    """
    Grafica el promedio de la energía suavizada de todas las réplicas.
    Esto genera bandas de desviación estándar más limpias y útiles.
    """
    if not dfs: return

    # Aplicamos el suavizado a cada réplica INDIVIDUALMENTE antes de promediar
    dfs_smooth = []
    for df in dfs:
        d = df.copy()
        d["Energy_Smooth"] = d["Potential_Energy"].rolling(window=(VENTANA_SUAVIZADO + 30), center=True).mean()
        dfs_smooth.append(d)

    df_concat = pd.concat(dfs_smooth, ignore_index=True)

    plt.figure(figsize=(12, 6))

    # Graficamos usando la columna suavizada
    sns.lineplot(
        data=df_concat,
        x="Time_ns",
        y="Energy_Smooth",
        color="blue",
        errorbar="sd", 
        label="Promedio Global (Suavizado) ± DE"
    )

    plt.title(f"Estabilidad Energética – {peptido} (n={len(dfs)})")
    plt.xlabel("Tiempo (ns)")
    plt.ylabel("Energía Potencial (kcal/mol)")
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.legend()
    plt.tight_layout()
    plt.savefig(carpeta / f"04_Energy_Promedio_{peptido}.png")
    plt.close()

# =======================
# MAIN
# =======================
def main():
    configurar_estilo()
    CARPETA_SALIDA.mkdir(exist_ok=True, parents=True)

    # 1. Cargar datos
    rmsd_1 = cargar_datos_rmsd_1(RUTA_DATOS, PEPTIDO, REPLICAS, DT)
    rmsd_2 = cargar_datos_rmsd_2(RUTA_DATOS, PEPTIDO, REPLICAS, DT)
    #rmsf = cargar_datos_rmsf(RUTA_DATOS, PEPTIDO, REPLICAS)
    contactos = cargar_datos_contacto(RUTA_DATOS, PEPTIDO, REPLICAS, DT)
    #energia = cargar_energia_desde_log(RUTA_DATOS, PEPTIDO, REPLICAS, DT)

    # 2. Generar gráficos si hay datos
    if rmsd_1:
        graficar_rmsd_individuales(rmsd_1, CARPETA_SALIDA, PEPTIDO)
        graficar_rmsd_promedio(rmsd_1, CARPETA_SALIDA, PEPTIDO, True)
    else:
        print("[AVISO] No se generaron gráficos RMSD (faltan datos).")

    if rmsd_2:
        graficar_rmsd_promedio(rmsd_2, CARPETA_SALIDA, PEPTIDO, False)
    else:
        print("[AVISO] No se generaron gráficos RMSD (faltan datos).")

    #if rmsf:
    #    graficar_rmsf_individuales(rmsf, CARPETA_SALIDA, PEPTIDO)
    #    graficar_rmsf_promedio(rmsf, CARPETA_SALIDA, PEPTIDO)
    #else:
    #    print("[AVISO] No se generaron gráficos RMSF (faltan datos).")

    if contactos:
        graficar_contactos_individuales(contactos, CARPETA_SALIDA, PEPTIDO, TIEMPO)
        graficar_contactos_promedio(contactos, CARPETA_SALIDA, PEPTIDO, TIEMPO)
    else:
        print("[AVISO] No se generaron gráficos de Contactos (faltan datos).")

    #if energia:
    #    graficar_energia_individuales(energia, CARPETA_SALIDA, PEPTIDO)
    #    graficar_energia_promedio(energia, CARPETA_SALIDA, PEPTIDO)
    #else:
    #    print("[AVISO] No se generaron gráficos de Energía (no se encontraron logs).")

    print(f"\nProceso finalizado. Gráficos guardados en: {CARPETA_SALIDA.resolve()}")

if __name__ == "__main__":
    main()
