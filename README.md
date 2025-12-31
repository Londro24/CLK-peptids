# Proyecto Nano AuNP // CLP
_Alejandro Ide Figueroa; Alumno de 5to año de Ingenieria Civil en Bioinformática_
## Introducción
Tras el trabajo realizado por _Galaz-Araya C._ en _Molecular Dynamics Study on the Influence of the CLK-Motif on the Structural Stability of Collagen-Like Peptides Adsorbed on Gold Nanosurfaces_, aparece la posibilidad de realizar un trabajo realicionado al mismo.
Con este conococimiento previo, en conjunto a la posibilidad de generar una mejor estabilidad del peptido por medio de la modificación de la estructura de los Collagen-Like Peptides (CLP) con aminoacidos como el triptofano (W), se realiza este proyecto.
## Desarrollo del peptido
Para poder empezar, se decidio generar diferentes modelos de peptidos para poder visualizar la interaccion del **W** desde la base de la secuencia aminoacidica Cisteína-Leucina-Lisina (CLK) generando los siguientes modelos
1) CLK
2) CLKW
3) CLKWW
4) CLKWWW
El primer modelo reflejaria la interacción base que deberia realizar, es decir, funcionaria de control. El resto de modelos reflejaria a su vez los interacción entre si mismo y el oro.
### Seleccion del modelo del peptido
En primera instancia se tomo desde el trabajo realizado anteriormente el archivo _pdb_ y se modifico la secuencia aminoacidica para obtener los 4 modelos en ChimeraX. Luego se construyo desde 0 en PyMol-openssl los 4 modelos mencionados. Posteriormente se uso el modelo de CABSflex para obtener la estructura del modelo 2 al 4. Y por ultimo se predijo la estructura de los modelos 3 y 4 por PEP-FOLD.
Se revisaron los 13 modelos hechos, a lo cual se decidio por literatura y eliminar errores provocados, por la construccion y modifición de una cadena pre-hecha, usar el modelo 4 creado por PEP-FOLD y eliminar los **W** hasta obtener los 4 modelos.
## Minimizacion y equilibrado
* Integrar datos y justifiacion de los mismos *
