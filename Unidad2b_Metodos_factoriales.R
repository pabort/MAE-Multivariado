
# Introducción -------------------------------------------------------------

# [setup]

for (loc in c("es_AR.UTF-8", "es_ES.UTF-8", "en_US.UTF-8", "C.UTF-8", "C.utf8")) {
  if (Sys.setlocale("LC_CTYPE", loc) != "") break
}

library(factoextra)    # visualización de ACP y correspondencias
library(FactoMineR)    # correspondencias simples y múltiples
library(ca)             # correspondencias, como contraste de FactoMineR
library(MASS)           # isoMDS (escalado no métrico)
library(psych)          # métodos de estimación del análisis factorial
library(GPArotation)    # rotación quartimax para factanal()
library(tidyverse)      # al final, para que dplyr::select() no quede enmascarado

theme_set(theme_minimal(base_size = 11))


# Carga de paquetes y datos ------------------------------------------------


## Los conjuntos de datos --------------------------------------------------

# [carga-datos]

url_base <- "https://raw.githubusercontent.com/pabort/MAE-Multivariado/main/"

ciudades <- read_csv(paste0(url_base, "ciudades_espania.csv"), show_col_types = FALSE) |>
  as.data.frame()
rownames(ciudades) <- ciudades$ciudad
ciudades <- ciudades[, -1]

productos <- read_csv(paste0(url_base, "productos.csv"), show_col_types = FALSE) |>
  as.data.frame()
rownames(productos) <- productos$producto
productos <- productos[, -1]

eph <- read_csv(paste0(url_base, "eph.csv"), show_col_types = FALSE) |>
  mutate(across(everything(), factor))

epf <- read_csv(paste0(url_base, "epf.csv"), show_col_types = FALSE) |>
  select(-provincia)

ciudades
productos
head(eph)
epf


# Escalado multidimensional ------------------------------------------------


## Objetivo y ejemplo motivador --------------------------------------------

# [mds-ciudades-datos]

ciudades


# [mds-ciudades-cmdscale]
#| fig-width: 6
#| fig-height: 5

mds_ciudades <- cmdscale(as.dist(ciudades), k = 2)
colnames(mds_ciudades) <- c("dim1", "dim2")

as_tibble(mds_ciudades, rownames = "ciudad") |>
  ggplot(aes(dim1, dim2, label = ciudad)) +
  geom_point(size = 2, colour = "#2166AC") +
  ggrepel::geom_text_repel() +
  coord_equal() +
  labs(title = "Escalado multidimensional — 6 ciudades españolas",
       x = "Coordenada principal 1", y = "Coordenada principal 2")


# [mds-eurodist-orientacion]

mds_eurodist <- cmdscale(eurodist, k = 2)

# Referencia para orientar los ejes: Estocolmo (norte) vs. Atenas (sur).
# Antes de invertir, Estocolmo tiene la segunda coordenada más negativa que
# Atenas, es decir que "sur" queda arriba: hay que invertir el eje 2 para
# que el norte quede hacia arriba, como en un mapa convencional.
orientacion_ns <- mds_eurodist["Stockholm", 2] < mds_eurodist["Athens", 2]
orientacion_ns


# [mds-eurodist-plot]
#| fig-width: 8
#| fig-height: 7

mds_eurodist_orientado <- mds_eurodist
if (orientacion_ns) mds_eurodist_orientado[, 2] <- -mds_eurodist_orientado[, 2]
colnames(mds_eurodist_orientado) <- c("dim1", "dim2")

as_tibble(mds_eurodist_orientado, rownames = "ciudad") |>
  ggplot(aes(dim1, dim2, label = ciudad)) +
  geom_point(size = 2, colour = "#2166AC") +
  ggrepel::geom_text_repel(size = 3) +
  coord_equal() +
  labs(title = "Escalado multidimensional — 21 ciudades europeas (eje norte-sur invertido)",
       x = "Coordenada principal 1 (aprox. oeste-este)",
       y = "Coordenada principal 2 (aprox. sur-norte)")


## Construcción de las coordenadas principales, paso a paso {#sec-mds-paso-a-paso} ----

# [mds-Q-matricial]

D2 <- as.matrix(ciudades)^2
n_ciudades <- nrow(D2)

P <- diag(n_ciudades) - matrix(1 / n_ciudades, n_ciudades, n_ciudades)
Q <- -0.5 * P %*% D2 %*% P

round(Q, 1)


# [mds-Q-elemento-a-elemento]

d2_fila <- rowMeans(D2)
d2_col  <- colMeans(D2)
d2_gral <- mean(D2)

Q_elemento <- outer(1:n_ciudades, 1:n_ciudades,
                     Vectorize(\(i, j) -0.5 * (D2[i, j] - d2_fila[i] - d2_col[j] + d2_gral)))
dimnames(Q_elemento) <- dimnames(Q)

chequeo_Q <- isTRUE(all.equal(Q, Q_elemento))
stopifnot(chequeo_Q)
chequeo_Q


# [mds-Y-manual]

eQ <- eigen(Q)

Y_manual <- eQ$vectors[, 1:2] %*% diag(sqrt(eQ$values[1:2]))
dimnames(Y_manual) <- list(rownames(ciudades), c("Dim1", "Dim2"))

round(Y_manual, 1)


# [mds-verificacion-cmdscale]

chequeo_cmdscale <- isTRUE(all.equal(abs(as.numeric(Y_manual)),
                                      abs(as.numeric(mds_ciudades)),
                                      tolerance = 1e-6))
stopifnot(chequeo_cmdscale)
chequeo_cmdscale


# [mds-verificacion-acp-epf]

mds_epf <- cmdscale(dist(scale(epf)), k = 2)
pca_epf <- prcomp(scale(epf))

chequeo_mds_acp <- isTRUE(all.equal(abs(as.numeric(mds_epf)),
                                     abs(as.numeric(pca_epf$x[, 1:2])),
                                     tolerance = 1e-6))
stopifnot(chequeo_mds_acp)
chequeo_mds_acp


## Compatibilidad con la métrica euclídea ----------------------------------

# [mds-autovalores-ciudades]

round(eQ$values, 1)
c(negativos = sum(eQ$values < -1e-6))


# [mds-autovalores-eurodist]

eig_eurodist <- cmdscale(eurodist, k = 2, eig = TRUE)$eig

round(eig_eurodist, 0)
c(negativos = sum(eig_eurodist < -1e-6))


# [mds-gof]

gof_ciudades  <- cmdscale(as.dist(ciudades), k = 2, eig = TRUE)$GOF
gof_eurodist  <- cmdscale(eurodist, k = 2, eig = TRUE)$GOF

tibble(medida = c("GOF[1] (denom. abs)", "GOF[2] (denom. positivos)"),
       ciudades = gof_ciudades, eurodist = gof_eurodist)


# [mds-compatibilidad-conclusion]
#| output: asis

neg_ciudades <- sum(eQ$values < -1e-6)
neg_eurodist <- sum(eig_eurodist < -1e-6)

cat(glue::glue(
  "Ninguna de las dos matrices es exactamente compatible con una métrica euclídea, pero en distinta ",
  "medida: `ciudades_espania` tiene {neg_ciudades} autovalor negativo, contra {neg_eurodist} en `eurodist`. ",
  "Tiene sentido que eurodist se aparte más: son distancias por carretera entre ciudades de toda Europa, ",
  "separadas por mares y cordilleras que ninguna ruta cruza en línea recta, mientras que las seis ciudades ",
  "españolas están todas en la misma península. Aun con autovalores negativos, `cmdscale()` sigue dando ",
  "una solución útil: descarta esas dimensiones y se queda con las de mayor autovalor positivo. Ahí se ve ",
  "la diferencia entre los dos ejemplos: la bondad de ajuste con 2 dimensiones (GOF[2]) es del ",
  "{round(100 * gof_ciudades[2])}% para las ciudades españolas, contra un {round(100 * gof_eurodist[2])}% ",
  "más bajo para eurodist — coherente con que esta última matriz se aparta más de la geometría euclídea."
))


## Escalado no métrico -----------------------------------------------------

# [mds-no-metrico-datos]

productos


# [mds-isomds]

d_productos <- as.dist(as.matrix(productos))
mds_no_metrico <- isoMDS(d_productos, k = 2, trace = FALSE)

mds_no_metrico$stress
mds_no_metrico$points


# [mds-no-metrico-plot]
#| fig-width: 6
#| fig-height: 5

mds_no_metrico_puntos <- mds_no_metrico$points
colnames(mds_no_metrico_puntos) <- c("dim1", "dim2")

as_tibble(mds_no_metrico_puntos, rownames = "producto") |>
  ggplot(aes(dim1, dim2, label = producto)) +
  geom_point(size = 2, colour = "#2166AC") +
  ggrepel::geom_text_repel() +
  coord_equal() +
  labs(title = "Escalado no métrico — 7 productos",
       x = "Dimensión 1", y = "Dimensión 2")


# [mds-shepard]
#| fig-width: 6
#| fig-height: 5

d_ajustada <- dist(mds_no_metrico$points)

shepard <- tibble(disimilitud = as.numeric(d_productos),
                   distancia_ajustada = as.numeric(d_ajustada)) |>
  arrange(disimilitud)

ajuste_monotono <- isoreg(shepard$disimilitud, shepard$distancia_ajustada)

ggplot(shepard, aes(disimilitud, distancia_ajustada)) +
  geom_point(colour = "#2166AC") +
  geom_step(aes(y = ajuste_monotono$yf), colour = "firebrick") +
  labs(title = "Diagrama de Shepard", subtitle = paste0("STRESS = ", round(mds_no_metrico$stress, 1), "%"),
       x = "Disimilitud original", y = "Distancia ajustada (2 dimensiones)")


# Análisis de correspondencias ---------------------------------------------


## Correspondencias simples, paso a paso -----------------------------------

# [ca-tabla]

tabla_ca <- table(eph$ESTCIV, eph$NIVEL_ED)
tabla_ca


# [ca-frecuencias]

n_ca <- sum(tabla_ca)
F_mat <- tabla_ca / n_ca

f_fila <- rowSums(F_mat)
f_col  <- colSums(F_mat)

round(F_mat, 3)


# [ca-perfiles-fila]

Df_inv <- diag(1 / f_fila)
R_perfiles <- Df_inv %*% F_mat
dimnames(R_perfiles) <- dimnames(tabla_ca)

round(R_perfiles, 3)


# [ca-matrices-transformadas]

Dc_inv_medio <- diag(1 / sqrt(f_col))
Df_inv_medio <- diag(1 / sqrt(f_fila))

Zf <- Df_inv %*% F_mat %*% Dc_inv_medio
Zc <- diag(1 / f_col) %*% t(F_mat) %*% Df_inv_medio

round(Zf, 3)


# [ca-autodescomposicion]

Z <- Df_inv_medio %*% F_mat %*% Dc_inv_medio

eZ <- eigen(t(Z) %*% Z)
round(eZ$values, 4)


# [ca-coordenadas-filas]

A <- eZ$vectors[, 2:3]
lambda_ca <- eZ$values[2:3]

Cf_manual <- Df_inv %*% F_mat %*% Dc_inv_medio %*% A
dimnames(Cf_manual) <- list(rownames(tabla_ca), c("Dim1", "Dim2"))

round(Cf_manual, 3)


# [ca-verificacion-factominer]

res_ca <- CA(tabla_ca, graph = FALSE)

chequeo_ca <- isTRUE(all.equal(abs(as.numeric(Cf_manual)), abs(as.numeric(res_ca$row$coord)),
                                tolerance = 1e-6))
stopifnot(chequeo_ca)
chequeo_ca


# [ca-verificacion-paquete-ca]

res_ca_pkg <- ca(tabla_ca)

chequeo_ca_pkg <- isTRUE(all.equal(res_ca_pkg$sv, sqrt(lambda_ca), tolerance = 1e-6))
stopifnot(chequeo_ca_pkg)
chequeo_ca_pkg


## Relación entre chi-cuadrado e inercia total -----------------------------

# [ca-chi2-inercia]

chi2 <- chisq.test(tabla_ca)$statistic
inercia_total <- sum(eZ$values[-1])

chequeo_chi2 <- isTRUE(all.equal(as.numeric(chi2) / n_ca, inercia_total, tolerance = 1e-6))
stopifnot(chequeo_chi2)

c(chi2 = as.numeric(chi2), inercia_total_x_n = inercia_total * n_ca)


## Interpretación y gráfico simétrico --------------------------------------

# [ca-resumen]

res_ca$eig
res_ca$row$contrib
res_ca$row$cos2


# [ca-biplot]
#| fig-width: 7
#| fig-height: 6

fviz_ca_biplot(res_ca, repel = TRUE) +
  labs(title = "Análisis de correspondencias — Estado civil × Nivel educativo")


# [ca-haireye]

tabla_haireye <- apply(HairEyeColor, c(1, 2), sum)
tabla_haireye

res_haireye <- CA(tabla_haireye, graph = FALSE)
res_haireye$eig


# [ca-haireye-biplot]
#| fig-width: 7
#| fig-height: 6

fviz_ca_biplot(res_haireye, repel = TRUE) +
  labs(title = "Análisis de correspondencias — Color de pelo × Color de ojos")


## Correspondencias múltiples ----------------------------------------------

# [mca-disyuntiva]

Z_disyuntiva <- model.matrix(~ SEXO + ESTCIV + NIVEL_ED + ESTADO - 1, data = eph)

dim(Z_disyuntiva)
head(Z_disyuntiva)


# [mca-burt]

burt <- t(Z_disyuntiva) %*% Z_disyuntiva

dim(burt)
burt[1:6, 1:6]


# [mca-ajuste]

res_mca <- MCA(eph, graph = FALSE)

res_mca$eig[1:5, ]


# [mca-biplot]
#| fig-width: 7
#| fig-height: 6

fviz_mca_var(res_mca, repel = TRUE) +
  labs(title = "Análisis de correspondencias múltiples — EPH")


# [transicion-verificacion]

Cc_ca <- res_ca$col$coord
Cf_desde_columnas <- sweep(Df_inv %*% F_mat %*% Cc_ca, 2, sqrt(lambda_ca), "/")

chequeo_transicion <- isTRUE(all.equal(as.matrix(Cf_desde_columnas), as.matrix(res_ca$row$coord),
                                        check.attributes = FALSE, tolerance = 1e-6))
stopifnot(chequeo_transicion)
chequeo_transicion


# Análisis Factorial Exploratorio ------------------------------------------

# [fa-datos]

epf


## El modelo y en qué difiere del ACP --------------------------------------


## Estructura de la varianza -----------------------------------------------

# [fa-comunalidad]

fa_ml <- factanal(epf, factors = 2, scores = "regression")

comunalidad <- rowSums(fa_ml$loadings[, ]^2)
especifica  <- fa_ml$uniquenesses

tibble(variable = names(epf), comunalidad, especifica, suma = comunalidad + especifica)


# [fa-comunalidad-verificacion]

# factanal() estima Lambda y Psi por optimización numérica (Newton-Raphson),
# no por una fórmula cerrada, así que la igualdad se cumple solo hasta la
# tolerancia de convergencia del algoritmo, no de forma exacta.
chequeo_comunalidad <- isTRUE(all.equal(as.numeric(comunalidad + especifica),
                                         rep(1, ncol(epf)), tolerance = 1e-4))
stopifnot(chequeo_comunalidad)
chequeo_comunalidad


## Métodos de estimación ---------------------------------------------------

# [fa-metodos-estimacion]

fa_cp <- principal(epf, nfactors = 2, rotate = "none")
fa_fp <- fa(epf, nfactors = 2, fm = "pa", rotate = "none")
# fa_ml ya se calculó arriba, con factanal()

tibble(
  variable          = names(epf),
  comunalidad_CP     = fa_cp$communality,
  comunalidad_FP     = fa_fp$communality,
  comunalidad_ML     = comunalidad
) |>
  mutate(across(where(is.numeric), ~round(.x, 3)))


## Número de factores ------------------------------------------------------

# [fa-numero-factores]

numero_factores <- tibble(
  m = 1:3,
  ajuste = map(m, \(m) factanal(epf, factors = m))
) |>
  mutate(
    chi2   = map_dbl(ajuste, \(a) a$STATISTIC),
    gl     = map_dbl(ajuste, \(a) a$dof),
    p_valor = map_dbl(ajuste, \(a) a$PVAL)
  ) |>
  select(-ajuste)

numero_factores


## Rotación ----------------------------------------------------------------

# [fa-rotaciones]

fa_varimax   <- factanal(epf, factors = 2, rotation = "varimax")
fa_quartimax <- factanal(epf, factors = 2, rotation = "quartimax")

round(cbind(sin_rotar = fa_ml$loadings[, 1],
            varimax   = fa_varimax$loadings[, 1],
            quartimax = fa_quartimax$loadings[, 1]), 3)


# [fa-rotacion-interpretacion]

round(fa_varimax$loadings[, 1:2], 3)


# [fa-rotacion-verificacion]

com_varimax   <- rowSums(fa_varimax$loadings[, ]^2)
com_quartimax <- rowSums(fa_quartimax$loadings[, ]^2)

chequeo_rotacion <- isTRUE(all.equal(comunalidad, com_varimax, tolerance = 1e-6)) &&
  isTRUE(all.equal(comunalidad, com_quartimax, tolerance = 1e-6))
stopifnot(chequeo_rotacion)
chequeo_rotacion


## Puntuaciones factoriales ------------------------------------------------

# [fa-puntuaciones]

fa_bartlett <- factanal(epf, factors = 2, scores = "Bartlett")

puntuaciones <- tibble(
  provincia       = seq_len(nrow(epf)),
  bartlett_f1     = fa_bartlett$scores[, 1],
  regresion_f1    = fa_ml$scores[, 1]
)

puntuaciones


# [fa-puntuaciones-correlacion]

cor(puntuaciones$bartlett_f1, puntuaciones$regresion_f1)


# [fa-puntuaciones-plot]
#| fig-width: 6
#| fig-height: 5

ggplot(puntuaciones, aes(bartlett_f1, regresion_f1)) +
  geom_point(colour = "#2166AC") +
  geom_abline(linetype = "dashed", colour = "grey50") +
  labs(title = "Puntuaciones factoriales — Bartlett vs. regresión",
       x = "Factor 1 (Bartlett)", y = "Factor 1 (regresión)")

