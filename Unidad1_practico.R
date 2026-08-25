
# Introducción -------------------------------------------------------------


## El conjunto de datos: Palmer Penguins -----------------------------------


# Carga de paquetes y datos ------------------------------------------------

# [setup]

# Aseguramos codificación UTF-8, para evitar problemas con tildes y símbolos (², χ², etc.)
# en salidas generadas dinámicamente (cat(), paste()) en sistemas con locale distinto.
# Probamos varias variantes de locale UTF-8, sin frenar si ninguna está instalada.
for (loc in c("es_AR.UTF-8", "es_ES.UTF-8", "en_US.UTF-8", "C.UTF-8", "C.utf8")) {
  if (Sys.setlocale("LC_CTYPE", loc) != "") break
}

# Paquetes para tests multivariados específicos (se cargan primero)
library(psych)      # test de Mardia (mardia()), liviano
library(biotools)  # test de Box M (carga MASS internamente)
library(ICSNP)     # T2 de Hotelling

# Paquetes para gráficos y manejo de datos (al final, por el conflicto de select())
library(GGally)     # ggpairs()
library(tidyverse)  # dplyr, ggplot2, readr, etc.

theme_set(theme_minimal(base_size = 11))


# [carga-datos]

url_penguins <- "https://raw.githubusercontent.com/allisonhorst/palmerpenguins/main/inst/extdata/penguins.csv"
penguins <- read_csv(url_penguins, show_col_types = FALSE)

glimpse(penguins)


# [limpieza]

datos_completos <- penguins |>
  select(species, bill_length_mm, bill_depth_mm, flipper_length_mm, body_mass_g) |>
  drop_na()

datos_completos |> count(species)


# [subset-adelie]

Datos <- datos_completos |>
  filter(species == "Adelie") |>
  select(-species) |>
  as.data.frame()

n <- nrow(Datos)
p <- ncol(Datos)

cat("n =", n, " p =", p, "\n")
Datos


# Análisis descriptivo -----------------------------------------------------


## Medidas descriptivas univariadas ----------------------------------------

# [descriptivas-univariadas]

Datos |>
  summarise(across(
    everything(),
    list(
      media    = mean,
      de       = sd,
      cv       = ~ sd(.x) / mean(.x),
      asimetria = moments::skewness,
      curtosis  = moments::kurtosis
    ),
    .names = "{.col}__{.fn}"
  )) |>
  pivot_longer(everything(), names_to = c("variable", "medida"), names_sep = "__") |>
  pivot_wider(names_from = medida, values_from = value)


# [fivenum]

# Mínimo, Q1, mediana, Q3, máximo por variable
sapply(Datos, fivenum) |>
  `rownames<-`(c("min", "Q1", "mediana", "Q3", "max"))


### Comparación entre especies ---------------------------------------------

# [descriptivas-por-especie]

datos_completos |>
  pivot_longer(cols = -species, names_to = "variable", values_to = "valor") |>
  group_by(species, variable) |>
  summarise(
    media = mean(valor),
    de    = sd(valor),
    cv    = de / media,
    .groups = "drop"
  )


# [boxplot-especies]
#| fig-width: 8
#| fig-height: 6

datos_completos |>
  pivot_longer(cols = -species, names_to = "variable", values_to = "valor") |>
  ggplot(aes(x = species, y = valor, fill = species)) +
  geom_boxplot(alpha = 0.7) +
  facet_wrap(~variable, scales = "free_y") +
  labs(x = NULL, y = NULL, fill = "Especie")


## Matrices de covarianza y correlación ------------------------------------

# [matrices-cov-cor]

S <- var(Datos)
R <- cor(Datos)

S
R


# [R-desde-S]

D <- diag(1 / sqrt(diag(S)))
R_reconstruida <- D %*% S %*% D
dimnames(R_reconstruida) <- dimnames(S)

R_reconstruida

all.equal(unname(R), unname(R_reconstruida))


## Medidas de variabilidad multivariada ------------------------------------

# [variabilidad-multivariada]

varianza_total   <- sum(diag(S))
varianza_media   <- varianza_total / p
varianza_gen     <- det(S)
desv_gen         <- sqrt(varianza_gen)
varianza_efec    <- varianza_gen^(1 / p)
desv_efec        <- sqrt(varianza_efec)
dependencia_efec <- 1 - det(R)^(1 / (p - 1))

tibble(
  medida = c("Varianza total", "Varianza media", "Varianza generalizada",
             "Desv. generalizada", "Varianza efectiva", "Desv. efectiva",
             "Dependencia efectiva"),
  valor  = c(varianza_total, varianza_media, varianza_gen,
             desv_gen, varianza_efec, desv_efec, dependencia_efec)
)


### Comparación entre especies ---------------------------------------------

# [variabilidad-por-especie]

datos_completos |>
  group_by(species) |>
  group_modify(~{
    Sg <- var(.x)
    Rg <- cor(.x)
    pg <- ncol(.x)
    tibble(
      n                 = nrow(.x),
      varianza_total    = sum(diag(Sg)),
      varianza_gen      = det(Sg),
      desv_gen          = sqrt(det(Sg)),
      varianza_efec     = det(Sg)^(1 / pg),
      desv_efec         = sqrt(det(Sg)^(1 / pg)),
      dependencia_efec  = 1 - det(Rg)^(1 / (pg - 1))
    )
  })


# [correlacion-parcial]

Sinv <- solve(S)                        # matriz de precisión
Dp   <- diag(diag(Sinv)^(-1 / 2))       # diag(S^-1)^(-1/2)

Rp <- -(Dp %*% Sinv %*% Dp)             # el signo cambiado va fuera de la diagonal
diag(Rp) <- 1
dimnames(Rp) <- dimnames(S)

round(Rp, 4)


# [verificacion-correlacion-parcial]

vars <- names(Datos)

expand.grid(i = 1:p, j = 1:p) |>
  filter(i < j) |>
  rowwise() |>
  mutate(
    par     = paste(vars[i], "vs", vars[j]),
    simple  = cor(Datos[[i]], Datos[[j]]),
    parcial_formula = Rp[i, j],
    parcial_residuos = {
      otras <- vars[-c(i, j)]
      f_i <- as.formula(paste(vars[i], "~", paste(otras, collapse = " + ")))
      f_j <- as.formula(paste(vars[j], "~", paste(otras, collapse = " + ")))
      cor(residuals(lm(f_i, data = Datos)), residuals(lm(f_j, data = Datos)))
    }
  ) |>
  ungroup() |>
  select(par, simple, parcial_formula, parcial_residuos)


# [equivalencia-pcor]

Rinv <- solve(R)

Rp_desde_R <- -(diag(diag(Rinv)^(-1 / 2)) %*% Rinv %*% diag(diag(Rinv)^(-1 / 2)))
diag(Rp_desde_R) <- 1
dimnames(Rp_desde_R) <- dimnames(R)

all.equal(Rp, Rp_desde_R)


# [condicionamiento]

c(condicion_S = kappa(S), condicion_R = kappa(R))


## Descomposición espectral ------------------------------------------------

# [descomposicion-espectral]

eS <- eigen(S)
eS

# La varianza total y generalizada también se recuperan a partir de los autovalores
sum(eS$values)   # = varianza total = tr(S)
prod(eS$values)  # = varianza generalizada = det(S)


## Estandarización multivariante -------------------------------------------

# [estandarizacion-multivariante]

U    <- eS$vectors
D1_2 <- diag(1 / sqrt(eS$values))

Xc  <- scale(Datos, center = TRUE, scale = FALSE)         # datos centrados
YM  <- Xc %*% U %*% D1_2 %*% t(U)                          # estandarización multivariante

round(colMeans(YM), 2)  # medias ~ 0
round(var(YM), 2)       # matriz de covarianzas ~ identidad


# Distancias multivariadas -------------------------------------------------


## Entre observaciones -----------------------------------------------------

# [distancia-euclidea-manual]

Datos_ej <- head(Datos, 6)

# Distancia euclídea "clásica"
dist(Datos_ej, method = "euclidean") |> as.matrix() |> round(2)


# [distancia-mahalanobis-manual]

S1 <- solve(S)
n_ej <- nrow(Datos_ej)
dih <- matrix(0, n_ej, n_ej)

for (i in 1:n_ej) {
  for (j in 1:n_ej) {
    diff <- as.matrix(Datos_ej[i, ] - Datos_ej[j, ])
    dih[i, j] <- sqrt(diff %*% S1 %*% t(diff))
  }
}
round(dih, 2)


# [distancia-mahalanobis-completa]

# Distancia euclídea completa (matriz n x n)
dist_eucl <- as.matrix(dist(Datos))

# La distancia máxima de cada observación al resto, y a qué observación corresponde
dist_max <- apply(dist_eucl, 2, max)
obs_mas_lejana <- apply(dist_eucl, 2, which.max)

tibble(observacion = 1:n, dist_max, obs_mas_lejana) |> head(10)


## Distancia de cada observación a la media: detección de atípicos ---------

# [mahalanobis-a-la-media]

media <- colMeans(Datos)
mh2   <- mahalanobis(Datos, media, S)  # distancia al cuadrado

corte <- qchisq(0.95, df = p)

resultado_atipicos <- tibble(
  observacion = 1:n,
  mh2 = mh2,
  p_valor = pchisq(mh2, df = p, lower.tail = FALSE),
  atipico = mh2 > corte
)

resultado_atipicos |> filter(atipico)


# [plot-mahalanobis]
#| fig-cap: "Distancia de Mahalanobis al cuadrado de cada observación respecto de la media. La línea horizontal marca el percentil 95 de una χ²(4)."

ggplot(resultado_atipicos, aes(x = observacion, y = mh2, color = atipico)) +
  geom_point() +
  geom_hline(yintercept = corte, linetype = "dashed") +
  scale_color_manual(values = c("FALSE" = "steelblue", "TRUE" = "firebrick")) +
  labs(x = "Observación", y = expression(D[M]^2), color = "¿Atípico?")


# Asimetría y curtosis multivariada ----------------------------------------

# [mardia-manual]

S1 <- solve(S)
Xc_mat <- as.matrix(scale(Datos, scale = FALSE))

# Matriz de productos g_ij = (x_i - media)' S^-1 (x_j - media)
G <- Xc_mat %*% S1 %*% t(Xc_mat)

asimetria_mardia <- sum(G^3) / n^2
curtosis_mardia  <- sum(diag(G)^2) / n

tibble(coeficiente = c("Asimetría de Mardia", "Curtosis de Mardia"),
       valor = c(asimetria_mardia, curtosis_mardia))


# Gráficos multivariados ---------------------------------------------------


## Matriz de dispersión ----------------------------------------------------

# [ggpairs]
#| fig-width: 8
#| fig-height: 7

ggpairs(Datos, progress = FALSE) +
  labs(title = "Matriz de dispersión — Pingüinos Adelie")


# [ggpairs-especies]
#| fig-width: 8
#| fig-height: 7

ggpairs(datos_completos, columns = 2:5, aes(color = species, alpha = 0.6),
        progress = FALSE) +
  labs(title = "Matriz de dispersión por especie")


## Gráficos de glifos ------------------------------------------------------

# [glifos]
#| fig-width: 7
#| fig-height: 7

set.seed(123)
muestra_glifos <- Datos[sample(1:n, 15), ]

stars(muestra_glifos, key.loc = c(4, 2),
      main = "Gráfico de estrellas - muestra de 15 pingüinos Adelie")


# [caras-chernoff]
#| fig-width: 7
#| fig-height: 7

TeachingDemos::faces(muestra_glifos,
                      main = "Caras de Chernoff - muestra de 15 pingüinos Adelie")


## Gráfico 3D --------------------------------------------------------------

# [scatter3d]
#| fig-width: 7
#| fig-height: 6

scatterplot3d::scatterplot3d(
  Datos$bill_length_mm, Datos$flipper_length_mm, Datos$body_mass_g,
  type = "h", angle = 55,
  xlab = "Largo de pico (mm)", ylab = "Largo de aleta (mm)", zlab = "Masa (g)",
  main = "Pingüinos Adelie"
)


# Distribución Normal Multivariada -----------------------------------------


## Visualización de la densidad normal bivariada ---------------------------

# [normal-bivariada]
#| fig-width: 8
#| fig-height: 7

normal_bivariada <- function(x, y, rho, mu1, sigma1, mu2, sigma2) {
  1 / (2 * pi * sigma1 * sigma2 * sqrt(1 - rho^2)) *
    exp(-1 / (2 * (1 - rho^2)) *
          (((x - mu1) / sigma1)^2 - 2 * rho * ((x - mu1) / sigma1) * ((y - mu2) / sigma2) +
             ((y - mu2) / sigma2)^2))
}

x <- seq(-4, 4, length.out = 60)
y <- seq(-4, 4, length.out = 60)

par(mfrow = c(2, 2), mar = c(1, 1, 4, 1))
for (rho in c(0.85, 0.5, 0, -0.85)) {
  f <- outer(x, y, normal_bivariada, rho = rho, mu1 = 0, sigma1 = 1, mu2 = 0, sigma2 = 1)
  persp(x, y, f, theta = 30, phi = 30, col = "lightblue",
        xlab = "X", ylab = "Y", zlab = "Z", main = paste("rho =", rho), cex.main = 0.9)
}
mtext("Distribución Normal Bivariada", outer = TRUE, line = -1.5, font = 2)


## Elipses de confianza ----------------------------------------------------

# [elipses]
#| fig-width: 8
#| fig-height: 5

par(mfrow = c(2, 3))
for (rho in c(0.70, -0.70, 0.85, 0, 0.5)) {
  plot(ellipse::ellipse(rho), type = "l", main = paste("rho =", rho))
}


## Evaluación del supuesto de normalidad multivariada ----------------------


### Gráfico QQ multivariado ------------------------------------------------

# [qq-multivariado]
#| fig-cap: "Gráfico QQ multivariado: distancias de Mahalanobis ordenadas vs. cuantiles de una χ²(4)."

qq_datos <- tibble(
  prob = (1:n - 0.5) / n,
  chi_teorico = qchisq(prob, df = p),
  mh2_ordenado = sort(mh2)
)

ggplot(qq_datos, aes(x = chi_teorico, y = mh2_ordenado)) +
  geom_point(color = "steelblue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(x = expression(chi[4]^2), y = "Distancias de Mahalanobis ordenadas",
       title = "Evaluación gráfica de normalidad multivariada")


### Test de Mardia ---------------------------------------------------------

# [mardia-test-manual]

asimetria_ajustada <- asimetria_mardia * n / 6
gl_asimetria       <- p * (p + 1) * (p + 2) / 6
p_valor_asimetria  <- pchisq(asimetria_ajustada, gl_asimetria, lower.tail = FALSE)

curtosis_z      <- (curtosis_mardia - p * (p + 2)) / sqrt(8 * p * (p + 2) / n)
p_valor_curtosis <- 2 * pnorm(abs(curtosis_z), lower.tail = FALSE)

tibble(
  componente = c("Asimetría", "Curtosis"),
  estadistico = c(asimetria_ajustada, curtosis_z),
  p_valor = c(p_valor_asimetria, p_valor_curtosis)
)


# [psych-mardia]

mardia(Datos, plot = FALSE)


# Inferencia multivariada --------------------------------------------------


## T² de Hotelling para una muestra ----------------------------------------

# [t2-una-muestra-manual]

mu0 <- c(bill_length_mm = 38.8, bill_depth_mm = 18.3,
         flipper_length_mm = 190, body_mass_g = 3700)

m <- colMeans(Datos)
T2 <- n * t(m - mu0) %*% solve(S) %*% (m - mu0)

F_crit  <- qf(0.95, p, n - p)
T2_crit <- (n - 1) * p * F_crit / (n - p)

tibble(estadistico = c("T2 observado", "T2 crítico (α=0.05)"),
       valor = c(as.numeric(T2), T2_crit))


# [t2-una-muestra-conclusion]
#| output: asis

if (T2 > T2_crit) {
  cat("Como T² observado supera al valor crítico, **se rechaza H0**.")
} else {
  cat("Como T² observado no supera al valor crítico, **no se rechaza H0**.")
}


# [t2-una-muestra-icsnp]

HotellingsT2(Datos, mu = mu0)


### Intervalos de confianza simultáneos ------------------------------------

# [intervalos-simultaneos]

# Basados en T2
li_t2 <- m - sqrt(T2_crit * diag(S) / n)
ls_t2 <- m + sqrt(T2_crit * diag(S) / n)

# Basados en Bonferroni
alpha_bonf <- 0.05 / (2 * p)
t_crit <- qt(1 - alpha_bonf, n - 1)
li_bonf <- m - t_crit * sqrt(diag(S) / n)
ls_bonf <- m + t_crit * sqrt(diag(S) / n)

tibble(
  variable = names(m),
  li_T2 = li_t2, ls_T2 = ls_t2, amplitud_T2 = ls_t2 - li_t2,
  li_Bonferroni = li_bonf, ls_Bonferroni = ls_bonf, amplitud_Bonferroni = ls_bonf - li_bonf
)


## T² de Hotelling para dos muestras ---------------------------------------

# [dos-grupos]

grupo1 <- datos_completos |> filter(species == "Adelie")    |> select(-species) |> as.data.frame()
grupo2 <- datos_completos |> filter(species == "Chinstrap") |> select(-species) |> as.data.frame()

n1 <- nrow(grupo1)
n2 <- nrow(grupo2)
g  <- 2

m1 <- colMeans(grupo1)
m2 <- colMeans(grupo2)

tibble(variable = names(m1), media_Adelie = m1, media_Chinstrap = m2, diferencia = m1 - m2)


# [t2-dos-muestras-manual]

c1 <- cov(grupo1)
c2 <- cov(grupo2)

# Matriz de covarianzas combinada (pooled)
Sp <- ((n1 - 1) * c1 + (n2 - 1) * c2) / (n1 + n2 - 2)

T2_dif <- t(m1 - m2) %*% ((n1 * n2 / (n1 + n2)) * solve(Sp)) %*% (m1 - m2)

F_crit_2m  <- qf(0.95, p, n1 + n2 - p - 1)
T2_crit_2m <- (n1 + n2 - 2) * p * F_crit_2m / (n1 + n2 - p - 1)

tibble(estadistico = c("T2 observado", "T2 crítico (α=0.05)"),
       valor = c(as.numeric(T2_dif), T2_crit_2m))


# [t2-dos-muestras-conclusion]
#| output: asis

if (T2_dif > T2_crit_2m) {
  cat("T² observado supera ampliamente al T² crítico: **se rechaza H0**. ",
      "Las especies Adelie y Chinstrap tienen vectores de medias morfológicas distintos — ",
      "resultado esperable, dado que son especies diferentes.")
} else {
  cat("No se rechaza H0.")
}


# [t2-dos-muestras-icsnp]

HotellingsT2(grupo1, grupo2)


## Igualdad de matrices de covarianza (test de Box M) ----------------------

# [box-m-manual]

ldc1 <- log(det(c1))
ldc2 <- log(det(c2))
ldsp <- log(det(Sp))

logBoxM <- ldc1 * (n1 - 1) / 2 + ldc2 * (n2 - 1) / 2 - ldsp * (n1 + n2 - g) / 2

gc1 <- (2 * p^2 + 3 * p - 1) * (1 / (n1 - 1) + 1 / (n2 - 1) - 1 / (n1 + n2 - g)) / (6 * (p + 1) * (g - 1))
gc2 <- (p - 1) * (p + 2) * (1 / (n1 - 1)^2 + 1 / (n2 - 1)^2 - 1 / (n1 + n2 - g)^2) / (6 * (g - 1))

v1 <- 0.5 * p * (p + 1) * (g - 1)
v2 <- (v1 + 2) / abs(gc2 - gc1^2)

b  <- (1 - gc1 - v1 / v2) / v1
b1 <- (1 - gc1 - 2 / v2) / v2

condicion <- gc2 - gc1^2
FBoxM <- if (condicion > 0) {
  -2 * b * logBoxM
} else {
  -(2 * b1 * v2 * logBoxM) / (v1 + 2 * b1 * v2 * logBoxM)
}

F_crit_box <- qf(0.95, v1, v2)
p_valor_box <- 1 - pf(FBoxM, v1, v2)

tibble(estadistico = c("F de Box M", "F crítico (α=0.05)", "p-valor"),
       valor = c(FBoxM, F_crit_box, p_valor_box))


# [box-m-biotools]

boxM(datos_completos |> filter(species %in% c("Adelie", "Chinstrap")) |> select(-species),
     datos_completos |> filter(species %in% c("Adelie", "Chinstrap")) |> pull(species))

