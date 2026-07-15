# TFM: Análisis Espaciotemporal Bayesiano de la Mortalidad por Ictus en EEUU

Repositorio creado con el fin de que mi Trabajo de Fin de Máster sea reproducible.

## Resumen
Este repositorio contiene el código y los datos necesarios para reproducir los análisis del Trabajo de Fin de Máster titulado **"Análisis espaciotemporal bayesiano de la mortalidad por ictus en Estados Unidos (2001-2020)"**. 

El trabajo aborda la distribución geográfica y temporal de la mortalidad por enfermedad cerebrovascular en los condados continentales del país a lo largo de dos décadas. Se adopta un marco de modelización jerárquica bayesiana basado en la verosimilitud de Poisson, utilizando la reparametrización BYM2 para la componente espacial y una interacción espaciotemporal de tipo no estructurado (Tipo I de Knorr-Held). La inferencia se realiza mediante la aproximación integrada de Laplace anidada (INLA) combinada con previas de complejidad penalizada.

## Estructura del Repositorio
* **`datos/`**: Directorio destinado a los archivos de datos espaciales, datos de mortalidad extraídos de CDC WONDER.
* **`figuras/`**: Carpeta donde se exportan los mapas de riesgo relativo y gráficos generados y el código para obtenerlos.
* **`Codigo.R`**: Script principal que contiene la carga de librerías, procesamiento espacial, modelización jerárquica bayesiana y extracción de resultados.
* **`INLA info.R`**: Script complementario con configuraciones o detalles técnicos específicos sobre el ajuste de los modelos gaussianos latentes.

## Requisitos y Entorno de Cálculo
El análisis ha sido desarrollado bajo el siguiente entorno:
* **R**: Versión 4.6.1.
* **Paquete principal**: `INLA` (versión 25.10.19).
* **Paquetes de manipulación y análisis espacial**: `tidyverse`, `sf`, `tigris`, `spdep`, `tidycensus`, `patchwork`, `Matrix`.

*Nota sobre el rendimiento:* La ejecución íntegra del análisis lleva aproximadamente 30 minutos en un equipo con 16 GB de RAM.

## Reproducibilidad
Dado que la metodología INLA aproxima las distribuciones analíticamente de forma determinista, la ejecución del script generará los mismos resultados sin necesidad de configurar una semilla aleatoria. Asegúrate de tener los paquetes en las versiones mencionadas para evitar problemas de compatibilidad.

## Autora
**María Pallares Diez**
Máster en Bioestadística - Universitat de València
