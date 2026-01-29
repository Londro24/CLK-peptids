# Proyecto Nano AuNP // CLP
_Alejandro A. Ide Figueroa; Alumno de 5to año de Ingenieria Civil en Bioinformática_

## Introducción
Tras el trabajo realizado por _Galaz-Araya C._ en _Molecular Dynamics Study on the Influence of the CLK-Motif on the Structural Stability of Collagen-Like Peptides Adsorbed on Gold Nanosurfaces_, aparece la posibilidad de realizar un trabajo realicionado al mismo.
Con este conococimiento previo, en conjunto a la posibilidad de generar una mejor estabilidad del peptido por medio de la modificación de la estructura de los Collagen-Like Peptides (CLP) con aminoacidos como el triptofano (W), se realiza este proyecto.

# Descripción de tareas realizadas

## Desarrollo del peptido
Para poder empezar, se decidio generar diferentes modelos de peptidos para poder visualizar la interaccion del **W** desde la base de la secuencia aminoacidica Cisteína-Leucina-Lisina (CLK) generando los siguientes modelos:
1) CLK
2) CLKW
3) CLKWW
4) CLKWWW

El primer modelo reflejaria la interacción base que deberia realizar, es decir, funcionaria de control. El resto de modelos reflejaria a su vez los interacción entre si mismo y el oro.

### Seleccion del modelo del peptido
En primera instancia se se modifico la secuencia aminoacidica, de los archivos _.pdb_ del trabajo anteriormente mencionado, para obtener los 4 modelos en ChimeraX. Luego se construyo desde 0 en PyMol-openssl los 4 modelos mencionados. Posteriormente se uso el servicio web de CABSflex para obtener la estructura del modelo 2 al 4. Y por ultimo, se predijo la estructura de los modelos 3 y 4 por el servicio web de PEP-FOLD.
Se revisaron los 13 modelos hechos, a lo cual se decidio por literatura y eliminar errores provocados, por la construccion y modifición de una cadena pre-hecha, usar el modelo 4 creado por PEP-FOLD y eliminar los **W** hasta obtener los 4 modelos.

## Construcción de los sistemas

### Superficie de oro
La superficie de oro es generada por medio de CHARM-GUI, con MM 111, tamaño de 50x50x10 Å³ y condiciones periodicas de borde en los ejes X e Y.

### Sistema
Con los modelos de los peptidos hechos y la superficie de oro, se realizo la construcción de los sistemas para cada uno de los peptidos con las caracteristicas similares a los descrito por el paper mencionado al inicio. Las caracteristicas que contemplan el sistema son: El peptido a 15 Å de distancia de la superficie de la lamina de oro, en una caja de agua neutralizada a pH 7.0 de tamaño generado por ambas estructuras, agregando 15 Å más de tamaño en el eje Z, centrado en el origen.

## Minimizacion y equilibrado

### Minimizado
El minimizado realizado para estos sistemas se llevo acabo con NAMD3_Multicore a 2 hilos en mi computador personal. Este proceso por cada sistema contemplo 20,000 pasos a _Time Step_ (ts) 4, con temperatura de 300K y 1 atmósfera de presión.

### Equilibrado
El equilibradro realizado para estos sistemas se llevo acabo con NAMD3_CUDAMulticore a 8 hilos en un computador del equipo de trabajo. Este proceso por cada sistema contemplo 20 ns (5,000,000 pasos) a 4 ts, con temperatura de 300K y 1 atmósfera de presión.

### Estado de ejecucion (Equilibrado)
| Sistema  | Pasos de ejecución | Estado |
|----------|--------------------|--------|
| 1_CLK    | 5,000,000 (20 ns)  | Listo  |
| 2_CLKW   | 5,000,000 (20 ns)  | Listo  |
| 3_CLKWW  | 5,000,000 (20 ns)  | Listo  |
| 4_CLKWWW | 5,000,000 (20 ns)  | Listo  |

## Producción
* Para este proceso se llevaron acabo 150 ns de producción, el cual fue distribuido equitativamente en 3 replicas de 50 ns respectivamnte

### Estado de ejecucion (Producción)
| Sistema  | Pasos de ejecución | Estado de replica        |
|----------|--------------------|--------------------------|
| 1_CLK    | 12,500,000 (50 ns) | Listos (Replica 1,2 y 3) |
| 2_CLKW   | 12,500,000 (50 ns) | Listos (Replica 1,2 y 3) |
| 3_CLKWW  | 12,500,000 (50 ns) | Listos (Replica 1,2 y 3) |
| 4_CLKWWW | 12,500,000 (50 ns) | Listos (Replica 1,2 y 3) |

## Análisis
Para este punto de desarrollaron scripts para llevarlo acabo entorno a el comportamiento y estabilidad de los peptidos en contacto con la lamina, especificamente buscando el RMSD, la distancia minima que tiene entre el peptido y la superficie, y el contacto que se puede llevar a cabo entre la C y el/los W

---
