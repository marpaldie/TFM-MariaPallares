# =============================================================================
# GENERACIÓN DE FIGURAS DE MAPAS CON DISPOSICIÓN VERTICAL 4×2
# (4 quinquenios en filas × 2 sexos en columnas)
# =============================================================================


# ── 0. LIBRERÍAS ─────────────────────────────────────────────────────────────

library(tidyverse)
library(sf)
library(tigris)
library(INLA)           

options(tigris_use_cache = TRUE)


# ── 1. RUTAS ─────────────────────────────────────────────────────────────────

if (requireNamespace("rstudioapi", quietly = TRUE) &&
    rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
}

path_datos   <- "datos"
path_figuras <- "figuras"
path_rdata   <- "todos_los_modelos.RData"

dir.create(path_figuras, showWarnings = FALSE)

if (!file.exists(path_rdata)) {
  stop(
    "No encuentro '", path_rdata, "' en el directorio de trabajo.\n",
    "Ejecuta primero el script original completo (Codigo.R) para generarlo,\n",
    "o copia el .RData a esta carpeta."
  )
}


# ── 2. PARÁMETROS GLOBALES ───────────────────────────────────────────────────

periodos <- c("2001-2005", "2006-2010", "2011-2015", "2016-2020")
PROJ_USA <- 5070    # Albers Equal Area (CONUS)
DPI      <- 300

ANCHO_VERT <- 5.9
ALTO_VERT  <- 8.7


# ── 3. DATOS CRUDOS (necesarios para la Fig 4.2 de SMR) ─────────────────────

cat("Leyendo datos CDC WONDER...\n")
archivos_q <- list.files(path = path_datos, pattern = "stroke_.*\\.tsv",
                         full.names = TRUE)
stopifnot(length(archivos_q) == 4)

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


# ── 4. CARTOGRAFÍA ───────────────────────────────────────────────────────────

cat("Cargando cartografía...\n")
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

carto         <- counties_modelo %>% st_transform(PROJ_USA)
carto_outline <- carto %>% summarise(geometry = st_union(geometry))


# ── 5. CARGA DE MODELOS Y CÁLCULO DE RR / P(RR>1) ───────────────────────────

cat("Cargando '", path_rdata, "'...\n", sep = "")
load(path_rdata)

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

calcular_prob_exceso <- function(modelo, datos, prefijo = "") {
  datos[[paste0("RR", prefijo)]] <- modelo$summary.fitted.values$mean
  datos[[paste0("prob", prefijo)]] <- sapply(
    modelo$marginals.fitted.values,
    function(x) 1 - inla.pmarginal(1, x)
  )
  datos
}

cat("Calculando RR y P(RR>1) en cada celda condado-periodo...\n")
datos_mujeres <- calcular_prob_exceso(modelo_st_mujeres, datos_mujeres, "_st")
datos_hombres <- calcular_prob_exceso(modelo_st_hombres, datos_hombres, "_st")
datos_mujeres <- calcular_prob_exceso(modelo_mujeres,    datos_mujeres, "_bym2")
datos_hombres <- calcular_prob_exceso(modelo_hombres,    datos_hombres, "_bym2")


# ── 6. TEMA Y ESCALAS GRÁFICAS PARA VERSIÓN VERTICAL ────────────────────────

theme_tfm_vert <- function(base_size = 11) {
  theme_void(base_size = base_size, base_family = "serif") +
    theme(
      plot.background  = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.title    = element_text(face = "bold", hjust = 0.5,
                                   margin = margin(b = 4)),
      plot.subtitle = element_text(size = base_size - 1, hjust = 0.5,
                                   color = "grey30", margin = margin(b = 8)),
      strip.text         = element_text(face = "bold", size = base_size),
      strip.text.y.left  = element_text(face = "bold", size = base_size,
                                        angle = 90),
      strip.placement    = "outside",
      legend.position   = "bottom",
      legend.key.width  = unit(1.4, "cm"),
      legend.key.height = unit(0.35, "cm"),
      legend.text  = element_text(size = base_size - 1),
      legend.title = element_text(size = base_size - 1)
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


# ── 7. FUNCIONES DE COMPOSICIÓN VERTICAL ────────────────────────────────────

mapa_facetado_vert <- function(datos_unidos, var, escala, titulo, subtitulo) {
  datos_unidos <- datos_unidos %>%
    filter(!is.na(periodo)) %>%
    mutate(periodo = factor(periodo, levels = periodos))
  ggplot(datos_unidos) +
    geom_sf(aes(fill = .data[[var]]), color = NA) +
    geom_sf(data = carto_outline, fill = NA, color = "grey20",
            linewidth = 0.2) +
    facet_grid(periodo ~ Sexo, switch = "y") +
    escala + theme_tfm_vert() +
    labs(title = titulo, subtitle = subtitulo)
}

construir_panel <- function(d_f, d_h, vars) {
  bind_rows(
    d_f %>% mutate(Sexo = "Mujeres"),
    d_h %>% mutate(Sexo = "Hombres")
  ) %>% select(GEOID, periodo, Sexo, all_of(vars)) %>%
    mutate(Sexo = factor(Sexo, levels = c("Mujeres", "Hombres")))
}


# ── 8. GENERAR LAS 5 FIGURAS ────────────────────────────────────────────────

# 4.2 SMR crudas -------------------------------------------------------------
cat("\nGenerando fig_4_2_smr.png ...\n")
smr_panel <- bind_rows(
  smr_q_sex %>% filter(Sex == "Female") %>% mutate(Sexo = "Mujeres"),
  smr_q_sex %>% filter(Sex == "Male")   %>% mutate(Sexo = "Hombres")
) %>% mutate(Sexo = factor(Sexo, levels = c("Mujeres", "Hombres")))
map_smr <- carto %>%
  left_join(smr_panel %>% select(GEOID, periodo, Sexo, SMR), by = "GEOID")
p_smr <- mapa_facetado_vert(map_smr %>% rename(RR = SMR), "RR", escala_rr,
                            "SMR crudas por condado, sexo y quinquenio",
                            "Mortalidad por ictus, estandarizacion indirecta")
ggsave(file.path(path_figuras, "fig_4_2_smr.png"), p_smr,
       width = ANCHO_VERT, height = ALTO_VERT, dpi = DPI, bg = "white")

# 4.3 RR del BYM2 espacial ---------------------------------------------------
cat("Generando fig_4_3_rr_bym2.png ...\n")
panel_bym2 <- construir_panel(datos_mujeres, datos_hombres,
                              c("RR_bym2", "prob_bym2"))
map_bym2 <- carto %>% left_join(panel_bym2, by = "GEOID")
p_rr_bym2 <- mapa_facetado_vert(map_bym2 %>% rename(RR = RR_bym2), "RR",
                                escala_rr,
                                "Riesgo relativo suavizado, modelo BYM2 espacial",
                                "Mortalidad por ictus en EE. UU. continental")
ggsave(file.path(path_figuras, "fig_4_3_rr_bym2.png"), p_rr_bym2,
       width = ANCHO_VERT, height = ALTO_VERT, dpi = DPI, bg = "white")

# 4.4 P(RR>1) del BYM2 espacial ---------------------------------------------
cat("Generando fig_4_4_prob_bym2.png ...\n")
p_prob_bym2 <- mapa_facetado_vert(map_bym2 %>% rename(prob = prob_bym2),
                                  "prob", escala_prob,
                                  "Probabilidad de exceso de riesgo, BYM2 espacial",
                                  "Umbrales de Richardson et al. (2004), 0,2 y 0,8")
ggsave(file.path(path_figuras, "fig_4_4_prob_bym2.png"), p_prob_bym2,
       width = ANCHO_VERT, height = ALTO_VERT, dpi = DPI, bg = "white")

# 4.5 RR del modelo final ----------------------------------------------------
cat("Generando fig_4_5_rr_final.png ...\n")
panel_st <- construir_panel(datos_mujeres, datos_hombres,
                            c("RR_st", "prob_st"))
map_st <- carto %>% left_join(panel_st, by = "GEOID")
p_rr_st <- mapa_facetado_vert(map_st %>% rename(RR = RR_st), "RR", escala_rr,
                              "Riesgo relativo suavizado, modelo final BYM2 + Tipo I",
                              "Mortalidad por ictus en EE. UU. continental")
ggsave(file.path(path_figuras, "fig_4_5_rr_final.png"), p_rr_st,
       width = ANCHO_VERT, height = ALTO_VERT, dpi = DPI, bg = "white")

# 4.6 P(RR>1) del modelo final ----------------------------------------------
cat("Generando fig_4_6_prob_final.png ...\n")
p_prob_st <- mapa_facetado_vert(map_st %>% rename(prob = prob_st),
                                "prob", escala_prob,
                                "Probabilidad de exceso de riesgo, modelo final",
                                "BYM2 + Tipo I, umbrales 0,2 y 0,8")
ggsave(file.path(path_figuras, "fig_4_6_prob_final.png"), p_prob_st,
       width = ANCHO_VERT, height = ALTO_VERT, dpi = DPI, bg = "white")


cat("\n========================================================================\n")
cat("Figuras regeneradas en '", path_figuras, "/' (versión vertical 4x2)\n", sep = "")
cat("========================================================================\n")
cat("Archivos a reemplazar en el Word del TFM:\n")
cat("  fig_4_2_smr.png\n")
cat("  fig_4_3_rr_bym2.png\n")
cat("  fig_4_4_prob_bym2.png\n")
cat("  fig_4_5_rr_final.png\n")
cat("  fig_4_6_prob_final.png\n\n")
cat("Sugerencia para Word: clic derecho sobre la imagen anterior > 'Cambiar\n")
cat("imagen' > seleccionar el PNG nuevo. Mantiene el pie de figura y la\n")
cat("posicion en el documento.\n")