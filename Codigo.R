# =============================================================================
# TFM Análisis espacio-temporal de mortalidad por ictus en EE. UU.
# Máster en Bioestadística · Universitat de València
# María Pallarés Díez · 2026
#
# Modelos ajustados:
#   1. BYM2 espacial puro (mujeres y hombres)
#   2. BYM2 + interacción Tipo I (mujeres y hombres)  -- MODELO FINAL
#   3. Modelo conjunto con covariables socioeconómicas e interacciones sexo×cov
#
# Decisión metodológica: no se incluyen los modelos con interacción Tipo II,
# III y IV de Knorr-Held por problemas de identificabilidad numérica y por
# la limitada información para estructuras RW1 con T=4 períodos (Goicoa et
# al., 2018). El modelo Tipo I se justifica como modelo final por su mejora
# clara frente al BYM2 espacial puro (DIC y WAIC) sin necesidad de
# estructuras temporales más complejas.
# =============================================================================


# ── 0. LIBRERÍAS Y AJUSTES GLOBALES ───────────────────────────────────────────

library(tidyverse)
library(sf)
library(tigris)
library(spdep)
library(INLA)
library(Matrix)
library(patchwork)
library(tidycensus)

options(tigris_use_cache = TRUE)
inla.setOption(num.threads = "4:1")


# ── 1. CARGA DE DATOS CDC WONDER ──────────────────────────────────────────────

path_datos <- file.path(dirname(rstudioapi::getSourceEditorContext()$path), "datos")

archivos_q <- list.files(path = path_datos, pattern = "stroke_.*\\.tsv",
                         full.names = TRUE)
periodos <- c("2001-2005", "2006-2010", "2011-2015", "2016-2020")

datos_q <- map2_dfr(archivos_q, periodos, function(archivo, periodo) {
  read.delim(archivo, sep = "\t", header = TRUE, stringsAsFactors = FALSE) %>%
    filter(!is.na(County) & County != "") %>%
    mutate(periodo = periodo)
}) %>%
  mutate(
    Deaths     = as.numeric(Deaths),
    Population = as.numeric(Population),
    GEOID      = str_pad(as.character(County.Code), 5, pad = "0")
  ) %>%
  filter(
    !str_sub(GEOID, 1, 2) %in% c("02", "15"),
    !is.na(Deaths), !is.na(Population)
  )

cat("Filas cargadas:", nrow(datos_q), "\n")


# ── 2. TASAS NACIONALES DE REFERENCIA Y CASOS ESPERADOS ──────────────────────

ref_q_sex <- datos_q %>%
  group_by(Ten.Year.Age.Groups.Code, Sex.Code, periodo) %>%
  summarise(D_ref = sum(Deaths), P_ref = sum(Population), .groups = "drop") %>%
  mutate(tasa_ref = D_ref / P_ref)

smr_q_sex <- datos_q %>%
  filter(!GEOID %in% c("46113", "51515")) %>%
  left_join(ref_q_sex,
            by = c("Ten.Year.Age.Groups.Code", "Sex.Code", "periodo")) %>%
  mutate(esperados = Population * tasa_ref) %>%
  group_by(GEOID, County, periodo, Sex) %>%
  summarise(O = sum(Deaths), E = sum(esperados), .groups = "drop") %>%
  mutate(SMR = O / E,
         ID_tiempo = as.integer(factor(periodo, levels = periodos)))

cat("Condados únicos con datos:", n_distinct(smr_q_sex$GEOID), "\n")


# ── 3. CARTOGRAFÍA Y LISTA MAESTRA DE CONDADOS ───────────────────────────────

counties_sf <- counties(cb = TRUE, year = 2019, resolution = "20m") %>%
  filter(!STATEFP %in% c("02", "15", "72")) %>%
  mutate(GEOID = paste0(STATEFP, COUNTYFP))

condados_todos <- counties_sf %>%
  st_drop_geometry() %>%
  select(GEOID) %>%
  filter(!GEOID %in% c("46113", "51515")) %>%
  arrange(GEOID) %>%
  mutate(ID_area = row_number())

counties_modelo <- counties_sf %>%
  inner_join(condados_todos, by = "GEOID") %>%
  arrange(ID_area)

cat("Condados en lista maestra:", nrow(condados_todos), "\n")


# ── 4. GRAFO DE VECINDAD ─────────────────────────────────────────────────────

vecindad <- poly2nb(counties_modelo, queen = TRUE)
vecindad <- make.sym.nb(
  union.nb(vecindad,
           knn2nb(knearneigh(st_centroid(counties_modelo), k = 4)))
)
nb2INLA("mapa_todos.graph", vecindad)

grafo <- inla.read.graph("mapa_todos.graph")
stopifnot(grafo$n == nrow(condados_todos))


# ── 5. DATOS DEL MODELO (CUADRÍCULA COMPLETA) ────────────────────────────────

datos_modelo <- crossing(
  condados_todos,
  periodo = periodos,
  Sex     = c("Female", "Male")
) %>%
  mutate(ID_tiempo = as.integer(factor(periodo, levels = periodos))) %>%
  left_join(smr_q_sex %>% select(GEOID, periodo, Sex, O, E, SMR),
            by = c("GEOID", "periodo", "Sex")) %>%
  mutate(ID_st = (ID_tiempo - 1) * max(ID_area) + ID_area)

datos_mujeres <- datos_modelo %>% filter(Sex == "Female")
datos_hombres <- datos_modelo %>% filter(Sex == "Male")


# ── 6. COVARIABLES SOCIOECONÓMICAS (ACS) ─────────────────────────────────────

# census_api_key("TU_KEY", install = TRUE)  

covariables <- get_acs(
  geography = "county", year = 2019, survey = "acs5",
  variables = c(
    pct_pobreza    = "B17001_002",  total_pobreza  = "B17001_001",
    pct_sin_seguro = "B27010_017",  total_seguro   = "B27010_001",
    pct_negro      = "B02001_003",  total_raza     = "B02001_001",
    mediana_renta  = "B19013_001"
  ),
  output = "wide"
) %>%
  select(GEOID, ends_with("E")) %>%
  rename_with(~str_remove(., "E$")) %>%
  mutate(
    prop_pobreza    = pct_pobreza    / total_pobreza,
    prop_sin_seguro = pct_sin_seguro / total_seguro,
    prop_negro      = pct_negro      / total_raza
  ) %>%
  select(GEOID, mediana_renta, prop_pobreza, prop_sin_seguro, prop_negro) %>%
  mutate(across(c(mediana_renta, prop_pobreza, prop_sin_seguro, prop_negro),
                ~scale(.)[, 1], .names = "{.col}_z"))

datos_mujeres_cov <- datos_mujeres %>% left_join(covariables, by = "GEOID")
datos_hombres_cov <- datos_hombres %>% left_join(covariables, by = "GEOID")


# =============================================================================
# MODELOS
# =============================================================================

# ── 7. HIPER-PRIORS Y OPCIONES COMUNES ───────────────────────────────────────

hyper_bym2 <- list(
  phi  = list(prior = "pc",      param = c(0.5, 0.5)),
  prec = list(prior = "pc.prec", param = c(0.3, 0.01))
)
hyper_prec <- list(prec = list(prior = "pc.prec", param = c(0.3, 0.01)))

ctrl_compute   <- list(dic = TRUE, waic = TRUE, cpo = TRUE,
                       return.marginals.predictor = TRUE)
ctrl_predictor <- list(compute = TRUE)


# ── 8. MODELO BYM2 ESPACIAL PURO ─────────────────────────────────────────────

formula_bym2 <- O ~ factor(ID_tiempo) +
  f(ID_area, model = "bym2", graph = "mapa_todos.graph",
    scale.model = TRUE, constr = TRUE, hyper = hyper_bym2)

cat("\nAjustando BYM2 espacial - mujeres...\n")
t0 <- Sys.time()
modelo_mujeres <- inla(formula_bym2, family = "poisson",
                       data = datos_mujeres, E = E,
                       control.compute = ctrl_compute, control.predictor = ctrl_predictor)
cat("  Tiempo:", round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n")

cat("Ajustando BYM2 espacial - hombres...\n")
t0 <- Sys.time()
modelo_hombres <- inla(formula_bym2, family = "poisson",
                       data = datos_hombres, E = E,
                       control.compute = ctrl_compute, control.predictor = ctrl_predictor)
cat("  Tiempo:", round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n")


# ── 9. MODELO BYM2 + INTERACCIÓN TIPO I (MODELO FINAL) ───────────────────────

formula_t1 <- O ~ factor(ID_tiempo) +
  f(ID_area, model = "bym2", graph = "mapa_todos.graph",
    scale.model = TRUE, constr = TRUE, hyper = hyper_bym2) +
  f(ID_st, model = "iid", constr = TRUE, hyper = hyper_prec)

cat("\nAjustando Tipo I - mujeres...\n")
t0 <- Sys.time()
modelo_st_mujeres <- inla(formula_t1, family = "poisson",
                          data = datos_mujeres, E = E,
                          control.compute = ctrl_compute, control.predictor = ctrl_predictor)
cat("  Tiempo:", round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n")

cat("Ajustando Tipo I - hombres...\n")
t0 <- Sys.time()
modelo_st_hombres <- inla(formula_t1, family = "poisson",
                          data = datos_hombres, E = E,
                          control.compute = ctrl_compute, control.predictor = ctrl_predictor)
cat("  Tiempo:", round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n")


# ── 10. MODELO CONJUNTO CON COVARIABLES E INTERACCIONES POR SEXO ─────────────

datos_ambos <- bind_rows(datos_mujeres_cov, datos_hombres_cov) %>%
  mutate(sexo_hombre = as.integer(Sex == "Male"))

formula_conjunta <- O ~
  factor(ID_tiempo) + sexo_hombre +
  mediana_renta_z + prop_pobreza_z + prop_sin_seguro_z + prop_negro_z +
  sexo_hombre:prop_pobreza_z + sexo_hombre:prop_negro_z +
  sexo_hombre:prop_sin_seguro_z +
  f(ID_area, model = "bym2", graph = "mapa_todos.graph",
    scale.model = TRUE, constr = TRUE, hyper = hyper_bym2) +
  f(ID_st, model = "iid", constr = TRUE, hyper = hyper_prec)

cat("\nAjustando modelo conjunto con covariables...\n")
t0 <- Sys.time()
modelo_conjunto <- inla(formula_conjunta, family = "poisson",
                        data = datos_ambos, E = E,
                        control.compute = ctrl_compute, control.predictor = ctrl_predictor)
cat("  Tiempo:", round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n")


# =============================================================================
# EXTRACCIÓN DE RESULTADOS PARA LAS TABLAS DEL TFM
# =============================================================================

# ── 11. TABLA COMPARATIVA DE MODELOS (Tablas 4.1, 4.2 y 4.3 del TFM) ─────────

resumen_modelo <- function(m, etiqueta, sexo) {
  phi <- tryCatch(m$summary.hyperpar["Phi for ID_area", "mean"],
                  error = function(e) NA_real_)
  sigma_bym2 <- tryCatch({
    prec_mean <- m$summary.hyperpar["Precision for ID_area", "mean"]
    1 / sqrt(prec_mean)
  }, error = function(e) NA_real_)
  tibble(
    Modelo = etiqueta,
    Sexo   = sexo,
    DIC    = round(m$dic$dic, 1),
    p_DIC  = round(m$dic$p.eff, 1),
    WAIC   = round(m$waic$waic, 1),
    p_WAIC = round(m$waic$p.eff, 1),
    LCPO   = round(-mean(log(m$cpo$cpo), na.rm = TRUE), 4),
    sigma  = round(sigma_bym2, 3),
    Phi    = round(phi, 3)
  )
}

tabla_modelos <- bind_rows(
  resumen_modelo(modelo_mujeres,    "BYM2 espacial",  "Mujeres"),
  resumen_modelo(modelo_st_mujeres, "BYM2 + Tipo I",  "Mujeres"),
  resumen_modelo(modelo_hombres,    "BYM2 espacial",  "Hombres"),
  resumen_modelo(modelo_st_hombres, "BYM2 + Tipo I",  "Hombres")
)

print(tabla_modelos, n = 10)
write.csv(tabla_modelos, "tabla_comparativa_modelos.csv", row.names = FALSE)


# ── 12. COEFICIENTES DEL MODELO CONJUNTO (Tabla 4.4 del TFM) ─────────────────

tabla_covariables <- as_tibble(modelo_conjunto$summary.fixed, rownames = "term") %>%
  filter(!str_starts(term, "factor"), term != "(Intercept)") %>%
  mutate(
    significativo = ifelse(sign(`0.025quant`) == sign(`0.975quant`), "Si", "No"),
    mean = round(mean, 4),
    q025 = round(`0.025quant`, 4),
    q975 = round(`0.975quant`, 4)
  ) %>%
  select(term, mean, q025, q975, significativo)

print(tabla_covariables)
write.csv(tabla_covariables, "tabla_covariables.csv", row.names = FALSE)


# ── 13. TASA NACIONAL POR QUINQUENIO Y SEXO (Figura 4.1) ─────────────────────

tasa_nacional <- datos_q %>%
  group_by(periodo, Sex) %>%
  summarise(
    muertes   = sum(Deaths, na.rm = TRUE),
    personas  = sum(Population, na.rm = TRUE),
    tasa_100k = muertes / personas * 100000,
    .groups   = "drop"
  )

print(tasa_nacional)
write.csv(tasa_nacional, "tasa_nacional.csv", row.names = FALSE)


# ── 14. RR Y P(RR > 1) DEL MODELO FINAL (BYM2 + Tipo I) Y DEL BYM2 ───────────

calcular_prob_exceso <- function(modelo, datos, prefijo = "") {
  datos[[paste0("RR", prefijo)]] <- modelo$summary.fitted.values$mean
  datos[[paste0("prob", prefijo)]] <- sapply(
    modelo$marginals.fitted.values,
    function(x) 1 - inla.pmarginal(1, x)
  )
  datos
}

datos_mujeres <- calcular_prob_exceso(modelo_st_mujeres, datos_mujeres, "_st")
datos_hombres <- calcular_prob_exceso(modelo_st_hombres, datos_hombres, "_st")
datos_mujeres <- calcular_prob_exceso(modelo_mujeres,    datos_mujeres, "_bym2")
datos_hombres <- calcular_prob_exceso(modelo_hombres,    datos_hombres, "_bym2")


# =============================================================================
# GENERACIÓN DE FIGURAS PARA EL TFM (A4 con márgenes 3 cm)
# =============================================================================

dir.create("figuras", showWarnings = FALSE)
ANCHO    <- 6.3
DPI      <- 300
PROJ_USA <- 5070

carto <- counties_modelo %>% st_transform(PROJ_USA)
carto_outline <- carto %>% summarise(geometry = st_union(geometry))

theme_tfm <- function(base_size = 9) {
  theme_void(base_size = base_size, base_family = "serif") +
    theme(
      plot.title    = element_text(face = "bold", hjust = 0.5,
                                   margin = margin(b = 4)),
      plot.subtitle = element_text(size = base_size - 1, hjust = 0.5,
                                   color = "grey30", margin = margin(b = 6)),
      strip.text    = element_text(face = "bold", size = base_size),
      legend.position   = "bottom",
      legend.key.width  = unit(1.2, "cm"),
      legend.key.height = unit(0.3, "cm")
    )
}

escala_rr <- scale_fill_gradientn(
  colours = c("#2166AC", "#67A9CF", "#D1E5F0", "#F7F7F7",
              "#FDDBC7", "#EF8A62", "#B2182B"),
  values = scales::rescale(log(c(0.5, 0.7, 0.85, 1, 1.15, 1.4, 2))),
  limits = c(0.5, 2), oob = scales::squish, trans = "log",
  breaks = c(0.5, 0.75, 1, 1.5, 2),
  labels = c("0.5", "0.75", "1", "1.5", ">=2"),
  name = "RR", na.value = "grey85"
)

escala_prob <- scale_fill_stepsn(
  colours = c("#2166AC", "#67A9CF", "#F7F7F7", "#EF8A62", "#B2182B"),
  breaks = c(0, 0.2, 0.8, 0.95, 1), limits = c(0, 1),
  name = "P(RR > 1)", na.value = "grey85"
)


# --- Figura 2.1: grafo de vecindad ------------------------------------------

centroides <- carto %>% st_centroid()
edges <- vector("list", 0); k <- 1
for (i in seq_along(vecindad)) {
  v <- vecindad[[i]]
  if (length(v) > 0 && v[1] != 0) {
    for (j in v[v > i]) {
      edges[[k]] <- st_linestring(rbind(
        st_coordinates(centroides$geometry[i]),
        st_coordinates(centroides$geometry[j])))
      k <- k + 1
    }
  }
}
aristas <- st_as_sf(data.frame(id = seq_along(edges)),
                    geometry = st_sfc(edges, crs = st_crs(carto)))

p_grafo <- ggplot() +
  geom_sf(data = carto_outline, fill = "grey97", color = NA) +
  geom_sf(data = aristas, color = "#5B8FA8", linewidth = 0.08, alpha = 0.4) +
  geom_sf(data = centroides, color = "#1F3864", size = 0.15, alpha = 0.6) +
  geom_sf(data = carto_outline, fill = NA, color = "grey20", linewidth = 0.3) +
  theme_tfm() +
  labs(title    = "Grafo de vecindad de los 3 108 condados continentales",
       subtitle = "Vecindad queen + k=4 vecinos mas proximos")
ggsave("figuras/fig_2_1_grafo.png", p_grafo,
       width = ANCHO, height = 4.0, dpi = DPI)


# --- Figura 4.1: evolución temporal -----------------------------------------

tasa_nacional_plot <- tasa_nacional %>%
  mutate(periodo_mid = case_when(
    periodo == "2001-2005" ~ 2003, periodo == "2006-2010" ~ 2008,
    periodo == "2011-2015" ~ 2013, periodo == "2016-2020" ~ 2018
  ))

p_tasa <- ggplot(tasa_nacional_plot,
                 aes(periodo_mid, tasa_100k, color = Sex, group = Sex)) +
  geom_line(linewidth = 0.6) + geom_point(size = 2.2) +
  scale_color_manual(values = c("Female" = "#C04050", "Male" = "#1F3864"),
                     labels = c("Mujeres", "Hombres"), name = NULL) +
  scale_x_continuous(breaks = c(2003, 2008, 2013, 2018),
                     labels = c("2001-2005","2006-2010","2011-2015","2016-2020")) +
  scale_y_continuous(limits = c(0, NA)) +
  labs(title    = "Tasa nacional de mortalidad por ictus (>= 35 anos)",
       subtitle = "Por 100 000 personas-ano, CDC WONDER, I60-I69",
       x = NULL, y = "Tasa por 100 000") +
  theme_minimal(base_size = 10, base_family = "serif") +
  theme(plot.title    = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5, color = "grey30"),
        legend.position = "top", panel.grid.minor = element_blank())
ggsave("figuras/fig_4_1_tasa.png", p_tasa,
       width = ANCHO, height = 3.6, dpi = DPI)


# --- Figuras 4.2 a 4.6: mapas facetados -------------------------------------

mapa_facetado <- function(datos_unidos, var, escala, titulo, subtitulo) {
  ggplot(datos_unidos %>% filter(!is.na(periodo))) +
    geom_sf(aes(fill = .data[[var]]), color = NA) +
    geom_sf(data = carto_outline, fill = NA, color = "grey20", linewidth = 0.2) +
    facet_grid(Sexo ~ periodo, switch = "y") +
    escala + theme_tfm() + theme(strip.placement = "outside") +
    labs(title = titulo, subtitle = subtitulo)
}

construir_panel <- function(d_f, d_h, vars) {
  bind_rows(
    d_f %>% mutate(Sexo = "Mujeres"),
    d_h %>% mutate(Sexo = "Hombres")
  ) %>% select(GEOID, periodo, Sexo, all_of(vars)) %>%
    mutate(Sexo = factor(Sexo, levels = c("Mujeres", "Hombres")))
}

# 4.2 SMR crudas
smr_panel <- bind_rows(
  smr_q_sex %>% filter(Sex == "Female") %>% mutate(Sexo = "Mujeres"),
  smr_q_sex %>% filter(Sex == "Male")   %>% mutate(Sexo = "Hombres")
) %>% mutate(Sexo = factor(Sexo, levels = c("Mujeres", "Hombres")))
map_smr <- carto %>%
  left_join(smr_panel %>% select(GEOID, periodo, Sexo, SMR), by = "GEOID")
p_smr <- mapa_facetado(map_smr %>% rename(RR = SMR), "RR", escala_rr,
                       "SMR crudas por condado, sexo y quinquenio",
                       "Mortalidad por ictus, estandarizacion indirecta")
ggsave("figuras/fig_4_2_smr.png", p_smr,
       width = ANCHO, height = 4.2, dpi = DPI)

# 4.3 RR del BYM2 espacial
panel_bym2 <- construir_panel(datos_mujeres, datos_hombres,
                              c("RR_bym2", "prob_bym2"))
map_bym2 <- carto %>% left_join(panel_bym2, by = "GEOID")
p_rr_bym2 <- mapa_facetado(map_bym2 %>% rename(RR = RR_bym2), "RR", escala_rr,
                           "Riesgo relativo suavizado, modelo BYM2 espacial",
                           "Mortalidad por ictus en EE. UU. continental")
ggsave("figuras/fig_4_3_rr_bym2.png", p_rr_bym2,
       width = ANCHO, height = 4.2, dpi = DPI)

# 4.4 P(RR>1) del BYM2 espacial
p_prob_bym2 <- mapa_facetado(map_bym2 %>% rename(prob = prob_bym2),
                             "prob", escala_prob,
                             "Probabilidad de exceso de riesgo, BYM2 espacial",
                             "Umbrales de Richardson et al. (2004), 0,2 y 0,8")
ggsave("figuras/fig_4_4_prob_bym2.png", p_prob_bym2,
       width = ANCHO, height = 4.2, dpi = DPI)

# 4.5 RR del modelo final (BYM2 + Tipo I)
panel_st <- construir_panel(datos_mujeres, datos_hombres,
                            c("RR_st", "prob_st"))
map_st <- carto %>% left_join(panel_st, by = "GEOID")
p_rr_st <- mapa_facetado(map_st %>% rename(RR = RR_st), "RR", escala_rr,
                         "Riesgo relativo suavizado, modelo final BYM2 + Tipo I",
                         "Mortalidad por ictus en EE. UU. continental")
ggsave("figuras/fig_4_5_rr_final.png", p_rr_st,
       width = ANCHO, height = 4.2, dpi = DPI)

# 4.6 P(RR>1) del modelo final
p_prob_st <- mapa_facetado(map_st %>% rename(prob = prob_st),
                           "prob", escala_prob,
                           "Probabilidad de exceso de riesgo, modelo final",
                           "BYM2 + Tipo I, umbrales 0,2 y 0,8")
ggsave("figuras/fig_4_6_prob_final.png", p_prob_st,
       width = ANCHO, height = 4.2, dpi = DPI)


# ── 15. GUARDAR TODOS LOS MODELOS ────────────────────────────────────────────

save(modelo_mujeres, modelo_hombres,
     modelo_st_mujeres, modelo_st_hombres,
     modelo_conjunto,
     tabla_modelos, tabla_covariables, tasa_nacional,
     file = "todos_los_modelos.RData")


# ── 16. RESUMEN FINAL PARA RELLENAR LA MEMORIA ───────────────────────────────

cat("\n========================================================================\n")
cat("VALORES PARA RELLENAR LAS TABLAS DEL TFM\n")
cat("========================================================================\n\n")

cat("--- Tabla 4.1 (BYM2 espacial puro) ---\n")
cat("  Mujeres:  Phi =", round(modelo_mujeres$summary.hyperpar["Phi for ID_area", "mean"], 3),
    "  sigma =", round(1/sqrt(modelo_mujeres$summary.hyperpar["Precision for ID_area", "mean"]), 3),
    "  p_eff =", round(modelo_mujeres$dic$p.eff, 1), "\n")
cat("  Hombres:  Phi =", round(modelo_hombres$summary.hyperpar["Phi for ID_area", "mean"], 3),
    "  sigma =", round(1/sqrt(modelo_hombres$summary.hyperpar["Precision for ID_area", "mean"]), 3),
    "  p_eff =", round(modelo_hombres$dic$p.eff, 1), "\n\n")

cat("--- Tabla 4.2 y 4.3 (BYM2 espacial vs BYM2 + Tipo I) ---\n")
print(tabla_modelos, n = 10)
cat("\nNota: Tipos II, III y IV se descartan por problemas de identificabilidad\n")
cat("con T=4 periodos. Justificacion en seccion 4.3.2 y 5.5 del TFM.\n\n")

cat("--- Tabla 4.4 (coeficientes del modelo conjunto) ---\n")
print(tabla_covariables)
cat("\n")

cat("--- Figura 4.1 (tasa nacional) ---\n")
print(tasa_nacional)
cat("\n")

cat("Archivos generados:\n")
cat("  tabla_comparativa_modelos.csv\n")
cat("  tabla_covariables.csv\n")
cat("  tasa_nacional.csv\n")
cat("  todos_los_modelos.RData\n")
cat("  figuras/fig_2_1_grafo.png\n")
cat("  figuras/fig_4_1_tasa.png\n")
cat("  figuras/fig_4_2_smr.png\n")
cat("  figuras/fig_4_3_rr_bym2.png\n")
cat("  figuras/fig_4_4_prob_bym2.png\n")
cat("  figuras/fig_4_5_rr_final.png\n")
cat("  figuras/fig_4_6_prob_final.png\n")
