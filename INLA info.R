# =============================================================================
# PURPOSE:      Verificar la versión de R-INLA y la configuración realmente
#               utilizada en los modelos del TFM, para poder afirmarlo en el
#               texto con seguridad ante el tribunal.
# INPUT:        todos_los_modelos.RData  (objetos INLA ajustados)
# OUTPUT:       Consola: versión, defaults y estrategia por modelo.
#               Fichero verificacion_inla.txt con la salida completa.
# USAGE NOTES:  Ejecutar en el mismo equipo y con la misma versión de INLA
#               con la que se ajustaron los modelos.
#               Copiar la salida de los apartados 1 y 4 al texto del TFM.
# 1st AUTHOR    DATE: Maria Pallares, 2026
# =============================================================================

library(INLA)

sink("verificacion_inla.txt", split = TRUE)

cat("=====================================================================\n")
cat(" 1. VERSION EXACTA DE R-INLA\n")
cat("=====================================================================\n")
# Esto es lo que hay que poner en el texto donde ahora pone VERSION_INLA
inla.version()
cat("\nCadena corta para el texto: ", as.character(packageVersion("INLA")), "\n\n")


cat("=====================================================================\n")
cat(" 2. DEFAULTS DE control.inla EN ESTA VERSION\n")
cat("=====================================================================\n")
defaults <- inla.set.control.inla.default()

cat("strategy      =", defaults$strategy, "\n")
cat("int.strategy  =", defaults$int.strategy, "\n\n")

# Comprobacion explicita de lo que afirma el TFM
ok_strategy <- identical(defaults$strategy, "simplified.laplace")
ok_int      <- identical(defaults$int.strategy, "auto")

cat("El TFM afirma: strategy por defecto = 'simplified.laplace'  -> ",
    ifelse(ok_strategy, "CORRECTO", "REVISAR: es distinto"), "\n")
cat("El TFM afirma: int.strategy por defecto = 'auto'            -> ",
    ifelse(ok_int, "CORRECTO", "REVISAR: es distinto"), "\n\n")


cat("=====================================================================\n")
cat(" 3. TODOS LOS DEFAULTS (por si el tribunal pregunta por alguno)\n")
cat("=====================================================================\n")
str(defaults, max.level = 1)
cat("\n")


cat("=====================================================================\n")
cat(" 4. QUE USO REALMENTE CADA MODELO AJUSTADO\n")
cat("=====================================================================\n")
cat("La estrategia de integracion depende del numero de hiperparametros:\n")
cat("  int.strategy='auto'  ->  'grid' si n_hiper <= 2,  'ccd' si n_hiper > 2\n\n")

load("todos_los_modelos.RData")

# Ajustar estos nombres a los objetos reales del RData
modelos <- Filter(function(x) inherits(get(x), "inla"), ls())

if (length(modelos) == 0) {
  cat("No se han encontrado objetos de clase 'inla'. Revisar nombres con ls().\n")
} else {
  for (m in modelos) {
    fit    <- get(m)
    nhyper <- nrow(fit$summary.hyperpar)
    integr <- if (nhyper <= 2) "grid" else "ccd"
    
    cat("---------------------------------------------------------------\n")
    cat("Modelo: ", m, "\n")
    cat("  hiperparametros (", nhyper, "):\n", sep = "")
    for (h in rownames(fit$summary.hyperpar)) cat("     -", h, "\n")
    cat("  -> integracion efectiva: ", integr, "\n")
    
    # strategy que quedo registrada en el objeto ajustado
    st <- fit$.args$control.inla$strategy
    if (!is.null(st)) cat("  -> strategy registrada:  ", st, "\n")
    
    # diagnostico kld: discrepancia entre aproximaciones
    if (!is.null(fit$summary.random)) {
      klds <- unlist(lapply(fit$summary.random, function(d)
        if ("kld" %in% names(d)) max(d$kld, na.rm = TRUE) else NA))
      klds <- klds[!is.na(klds)]
      if (length(klds)) {
        cat("  -> kld maximo en efectos aleatorios: ",
            format(max(klds), scientific = TRUE, digits = 3), "\n")
        cat("     (valores muy pequenos indican que la aproximacion es fiable)\n")
      }
    }
    cat("\n")
  }
}


cat("=====================================================================\n")
cat(" 5. sessionInfo() COMPLETO (para el repositorio)\n")
cat("=====================================================================\n")
print(sessionInfo())

sink()

cat("\n\nListo. Revisar 'verificacion_inla.txt'.\n")
cat("Pasos siguientes:\n")
cat("  1. Copiar la version del apartado 1 y sustituir VERSION_INLA en el Word.\n")
cat("  2. Comprobar que el apartado 2 dice CORRECTO en las dos lineas.\n")
cat("  3. Comprobar que el apartado 4 da 2 hiperparametros en el modelo\n")
cat("     espacial y 3 en el espaciotemporal y el conjunto.\n")