
# Introducción -------------------------------------------------------------


## El conjunto de datos ----------------------------------------------------


# Carga de paquetes y datos ------------------------------------------------

# [setup]

# Aseguramos codificación UTF-8, para evitar problemas con tildes y símbolos
# en los títulos de los gráficos generados por funciones de R base.
for (loc in c("es_AR.UTF-8", "es_ES.UTF-8", "en_US.UTF-8", "C.UTF-8", "C.utf8")) {
  if (Sys.setlocale("LC_CTYPE", loc) != "") break
}

library(factoextra)   # visualización y extracción de resultados de ACP
library(pcaPP)        # componentes principales robustas
library(MultiGroupO)  # componentes principales comunes
library(tidyverse)    # al final, para que dplyr::select() no quede enmascarado

theme_set(theme_minimal(base_size = 11))


# [carga-datos]

url_indice <- "https://raw.githubusercontent.com/pabort/MAE-Multivariado/main/indice.csv"
indice <- read_csv(url_indice, show_col_types = FALSE)

# Los nombres de empresa pasan a ser nombres de fila (para etiquetar los gráficos)
indice <- as.data.frame(indice)
rownames(indice) <- indice$EMPRESA
indice <- indice[, -1]

# Nos quedamos con los 8 primeros indicadores (LIQACID ... RENTECO)
indice8 <- indice[, 2:9]

dim(indice8)
head(indice8)


## Exploración inicial -----------------------------------------------------

# [descriptivas]

indice8 |>
  pivot_longer(everything(), names_to = "variable", values_to = "valor") |>
  group_by(variable) |>
  summarise(
    media    = mean(valor),
    de       = sd(valor),
    varianza = var(valor),
    minimo   = min(valor),
    maximo   = max(valor)
  )


# [matriz-correlacion]
#| fig-width: 7
#| fig-height: 6

R <- cor(indice8)
round(R, 3)


# [corrplot]
#| fig-width: 7
#| fig-height: 6

R |>
  as_tibble(rownames = "var1") |>
  pivot_longer(-var1, names_to = "var2", values_to = "correlacion") |>
  ggplot(aes(var1, var2, fill = correlacion)) +
  geom_tile() +
  geom_text(aes(label = round(correlacion, 2)), size = 3) +
  scale_fill_gradient2(low = "#B2182B", mid = "white", high = "#2166AC",
                       midpoint = 0, limits = c(-1, 1)) +
  labs(x = NULL, y = NULL, fill = "r") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Las componentes paso a paso ----------------------------------------------


## Matriz de covarianzas y descomposición espectral ------------------------

# [covarianzas-eigen]

S  <- cov(indice8)
eS <- eigen(S)

eS$values                    # autovalores
round(eS$vectors, 3)         # autovectores (por columna)


## Variabilidad explicada --------------------------------------------------

# [variabilidad-explicada]

p  <- ncol(indice8)
VT <- sum(eS$values)

tibble(
  componente = paste0("CP", 1:p),
  autovalor  = eS$values,
  prop       = 100 * autovalor / VT,
  acumulada  = cumsum(prop)
)


# [verificacion-traza]

c(traza_S = sum(diag(S)), suma_autovalores = VT)


## Cálculo de las componentes ----------------------------------------------

# [calculo-componentes]

Xdes <- scale(indice8, center = TRUE, scale = FALSE)  # matriz en desvíos
Y    <- Xdes %*% eS$vectors                            # componentes principales

round(Y[1:8, 1:3], 3)   # primeras 8 empresas, primeras 3 componentes


# [propiedades-componentes]

round(colMeans(Y), 10)   # medias ~ 0

c(varianza_total_X   = sum(diag(S)),
  varianza_total_Y   = sum(diag(var(Y))))

c(varianza_generalizada_X = det(S),
  varianza_generalizada_Y = det(var(Y)))


# [componentes-incorrelacionadas]

round(var(Y), 4)


# ACP sobre la matriz de covarianzas ---------------------------------------

# [prcomp-cov]

pca_cov <- prcomp(indice8)
pca_cov


# [resumen-prcomp-cov]

summary(pca_cov)


# [scree-cov]
#| fig-width: 7
#| fig-height: 4.5

fviz_eig(pca_cov, addlabels = TRUE,
         main = "Gráfico de sedimentación (matriz de covarianzas)")


# [aporte-varianzas]

tibble(
  variable = names(indice8),
  varianza = diag(S),
  aporte_pct = 100 * varianza / sum(varianza)
) |>
  arrange(desc(aporte_pct))


# ACP sobre la matriz de correlación ---------------------------------------

# [prcomp-cor]

pca_cor <- prcomp(indice8, scale. = TRUE)
summary(pca_cor)


# [comparacion-cov-cor]

tibble(
  componente     = paste0("CP", 1:p),
  sobre_S_pct    = 100 * eS$values / sum(eS$values),
  sobre_R_pct    = 100 * pca_cor$sdev^2 / sum(pca_cor$sdev^2)
) |>
  mutate(across(where(is.numeric), ~round(.x, 2)))


## Selección del número de componentes -------------------------------------

# [scree-cor]
#| fig-width: 7
#| fig-height: 4.5

fviz_eig(pca_cor, addlabels = TRUE, ncp = p,
         main = "Gráfico de sedimentación (matriz de correlación)") +
  geom_hline(yintercept = 100 / p, linetype = "dashed", colour = "firebrick")


# [criterios-seleccion]

tibble(
  componente = paste0("CP", 1:p),
  autovalor  = pca_cor$sdev^2,
  prop_pct   = 100 * autovalor / p,
  acumulada  = cumsum(prop_pct),
  kaiser     = ifelse(autovalor >= 1, "retener", "descartar")
)


## Interpretación de las cargas --------------------------------------------

# [cargas]

round(pca_cor$rotation[, 1:4], 3)


# [grafico-cargas]
#| fig-width: 8
#| fig-height: 6

pca_cor$rotation[, 1:4] |>
  as_tibble(rownames = "variable") |>
  pivot_longer(-variable, names_to = "componente", values_to = "carga") |>
  ggplot(aes(x = reorder(variable, carga), y = carga, fill = carga > 0)) +
  geom_col() +
  coord_flip() +
  facet_wrap(~componente) +
  scale_fill_manual(values = c("TRUE" = "#2166AC", "FALSE" = "#B2182B")) +
  labs(x = NULL, y = "Carga", fill = "Signo positivo") +
  theme(legend.position = "none")


## Correlaciones entre variables y componentes -----------------------------

# [correlaciones-var-comp]

U  <- pca_cor$rotation
Ds <- diag(pca_cor$sdev)      # raíces de los autovalores

CorZY <- U %*% Ds
dimnames(CorZY) <- list(names(indice8), paste0("CP", 1:p))

round(CorZY[, 1:4], 3)


# [cos2-contribuciones]

CorZY2 <- CorZY^2

# Calidad de representación acumulada en las primeras 3 componentes
tibble(
  variable = names(indice8),
  CP1 = CorZY2[, 1],
  CP2 = CorZY2[, 2],
  CP3 = CorZY2[, 3],
  acumulado_3CP = CP1 + CP2 + CP3,
  suma_fila_total = rowSums(CorZY2)
) |>
  mutate(across(where(is.numeric), ~round(.x, 3)))


# [suma-columnas]

# La suma por columnas devuelve los autovalores
tibble(
  componente = paste0("CP", 1:p),
  suma_columna = colSums(CorZY2),
  autovalor    = pca_cor$sdev^2
) |>
  mutate(across(where(is.numeric), ~round(.x, 4)))


# [factoextra-var]

var_info <- get_pca_var(pca_cor)

round(var_info$cos2[, 1:3], 3)      # calidad de representación
round(var_info$contrib[, 1:3], 2)   # contribuciones (%)


# [fviz-contrib]
#| fig-width: 8
#| fig-height: 4

fviz_contrib(pca_cor, choice = "var", axes = 1) +
  labs(title = "Contribución de las variables a la CP1")


# [circulo-correlaciones]
#| fig-width: 6.5
#| fig-height: 6

fviz_pca_var(pca_cor, col.var = "cos2", repel = TRUE,
             gradient.cols = c("#B2182B", "#E8A33D", "#2166AC")) +
  labs(title = "Círculo de correlaciones")


## Biplot ------------------------------------------------------------------

# [biplot]
#| fig-width: 8
#| fig-height: 7

fviz_pca_biplot(pca_cor, repel = TRUE,
                col.var = "#B2182B", col.ind = "grey40",
                labelsize = 3) +
  labs(title = "Biplot — CP1 y CP2")


# [biplot-grupos]
#| fig-width: 8
#| fig-height: 7

grupo <- factor(indice$CONDICIO)

fviz_pca_biplot(pca_cor, repel = TRUE,
                col.ind = grupo, palette = c("#2166AC", "#B2182B"),
                addEllipses = TRUE, ellipse.type = "confidence",
                legend.title = "Condición",
                label = "var", col.var = "black", labelsize = 3) +
  labs(title = "Biplot por condición")


# [biplot-otros-planos]
#| fig-width: 8
#| fig-height: 5

fviz_pca_biplot(pca_cor, axes = c(1, 3), repel = TRUE,
                col.var = "#B2182B", col.ind = "grey40", labelsize = 3) +
  labs(title = "Biplot — CP1 y CP3")


# Componentes principales robustas -----------------------------------------

# [pca-robusto]

pca_rob <- PCAproj(indice8, k = 6, method = "mad",
                   CalcMethod = "eachobs",
                   center = l1median_NLM, scale = "mad")

summary(pca_rob)


# [scree-robusto]
#| fig-width: 7
#| fig-height: 4.5

screeplot(pca_rob, main = "Sedimentacion - componentes robustas")


# [cargas-robustas]

round(unclass(pca_rob$loadings)[, 1:4], 3)


# [biplot-robusto]
#| fig-width: 8
#| fig-height: 7

biplot(pca_rob, main = "Biplot - componentes robustas")


# [comparacion-robusto-clasico]

tibble(
  variable = names(indice8),
  clasica  = pca_cor$rotation[, 1],
  robusta  = unclass(pca_rob$loadings)[, 1]
) |>
  mutate(across(where(is.numeric), ~round(.x, 3)))


# Componentes principales comunes ------------------------------------------

# [cp-comunes]

x     <- indice8
grupo <- factor(indice$CONDICIO)

table(grupo)


# [mgpca]

# Matriz a diagonalizar simultáneamente
mat_diag <- new.cov(x, cls = grupo, A = diag(ncol(x)))

res_mgpca <- mgpca(mat_diag, as.matrix(x), cls = grupo,
                   Plot = FALSE, ncomp = 2,
                   center = TRUE, scale = TRUE)

names(res_mgpca)


# [mgpca-varianza]

# Porcentaje de varianza recuperada por cada dimensión común
res_mgpca$prop_expl_var


# [mgpca-cargas]

# Cargas comunes a todos los grupos
cargas_comunes <- res_mgpca$loadings
dimnames(cargas_comunes) <- list(names(indice8), paste0("Dim", 1:ncol(cargas_comunes)))

round(cargas_comunes[, 1:2], 3)


# [mgpca-grafico]
#| fig-width: 7
#| fig-height: 6

# Coordenadas de las observaciones sobre las dos primeras dimensiones comunes
tibble(
  Dim1  = res_mgpca$variates[, 1],
  Dim2  = res_mgpca$variates[, 2],
  grupo = grupo
) |>
  ggplot(aes(Dim1, Dim2, colour = grupo)) +
  geom_point(size = 2, alpha = 0.8) +
  stat_ellipse(level = 0.95) +
  scale_colour_manual(values = c("#2166AC", "#B2182B")) +
  labs(title = "Componentes principales comunes",
       x = "Dimensión común 1", y = "Dimensión común 2",
       colour = "Condición")

