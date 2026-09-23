
# Introducción -------------------------------------------------------------


## Los conjuntos de datos --------------------------------------------------


# Carga de paquetes y datos ------------------------------------------------

# [setup]

for (loc in c("es_AR.UTF-8", "es_ES.UTF-8", "en_US.UTF-8", "C.UTF-8", "C.utf8")) {
  if (Sys.setlocale("LC_CTYPE", loc) != "") break
}

library(cluster)      # dist, daisy, silhouette
library(fpc)           # calinhara
library(clValid)       # comparación de agrupamientos
library(mclust)        # adjustedRandIndex
library(factoextra)    # visualización (fviz_dend, fviz_cluster, fviz_nbclust)
library(tidyverse)      # al final, para que dplyr::select() no quede enmascarado

theme_set(theme_minimal(base_size = 11))


# [carga-datos]

data(USArrests)

url_penguins <- "https://raw.githubusercontent.com/allisonhorst/palmerpenguins/main/inst/extdata/penguins.csv"
penguins <- read_csv(url_penguins, show_col_types = FALSE) |>
  drop_na()

dim(USArrests)
head(USArrests)


# Análisis descriptivo por grupos ------------------------------------------


## Un ejemplo con dos grupos conocidos -------------------------------------

# [empresas-datos]

empresas <- tibble(
  empresa   = paste0("E", 1:8),
  inversion = c(16, 12, 10, 12, 45, 50, 45, 50),
  ventas    = c(10, 14, 22, 25, 10, 15, 25, 27)
) |>
  column_to_rownames("empresa")

grupo_empresas <- factor(rep(1:2, each = 4))

empresas


# [empresas-twb]

xbar <- colMeans(empresas)
T_mat <- crossprod(as.matrix(sweep(empresas, 2, xbar)))

W_grupo <- function(g) {
  datos_g <- empresas[grupo_empresas == g, ]
  crossprod(as.matrix(sweep(datos_g, 2, colMeans(datos_g))))
}

W_mat <- W_grupo(1) + W_grupo(2)
B_mat <- T_mat - W_mat

list(T = T_mat, W = W_mat, B = B_mat)


# [empresas-verificacion]

stopifnot(isTRUE(all.equal(T_mat, W_mat + B_mat)))

round(T_mat, 2)
round(W_mat + B_mat, 2)


# [empresas-sp]

n_empresas <- nrow(empresas)
G_empresas <- nlevels(grupo_empresas)

Sp_mat <- W_mat / (n_empresas - G_empresas)
Sp_mat


# Medidas de proximidad ----------------------------------------------------


## Variables métricas ------------------------------------------------------

# [distancias-metricas-datos]

# Trabajamos con las 4 primeras filas de USArrests, alcanza para comparar medidas
sub_arrests <- USArrests[1:4, ]
sub_arrests


# [distancia-euclidea]

dist(sub_arrests, method = "euclidean")
dist(sub_arrests, method = "euclidean")^2   # al cuadrado: evita la raíz, más liviana de calcular


# [distancia-manhattan]

dist(sub_arrests, method = "manhattan")


# [minkowski-convergencia]

d_chebyshev <- dist(sub_arrests, method = "maximum")

# La distancia de Minkowski se acerca a la de Tchebycheff a medida que k crece
# (no llega a converger exactamente por el límite de precisión numérica con exponentes grandes)
tibble(
  k = c(1, 2, 10, 50, 100),
  diferencia_maxima_vs_chebyshev = map_dbl(k, \(k) max(abs(dist(sub_arrests, method = "minkowski", p = k) - d_chebyshev)))
)

d_chebyshev


# [distancias-invariantes]

dist(sub_arrests, method = "canberra")

# Karl Pearson: R no la trae armada. Se ve en la fórmula que equivale a una euclídea
# sobre las variables ponderadas por 1/sqrt(s_jj), así que se arma reescalando primero
s_jj <- diag(var(USArrests))
sub_pearson <- sweep(sub_arrests, 2, sqrt(s_jj), FUN = "/")
dist(sub_pearson, method = "euclidean")

# Mahalanobis: tampoco está en dist(), pero S^-1 = L'L (descomposición de Cholesky),
# así que es la distancia euclídea de siempre sobre los datos transformados por L
S    <- cov(USArrests)
L    <- chol(solve(S))
d_maha <- dist(as.matrix(sub_arrests) %*% t(L))
round(as.matrix(d_maha), 3)


# [verificacion-mahalanobis]

Sinv <- solve(S)
maha_forma_cuadratica <- function(r, s) sqrt(t(r - s) %*% Sinv %*% (r - s))

pares <- combn(1:4, 2)
maha_manual <- apply(pares, 2, \(idx) {
  maha_forma_cuadratica(as.numeric(sub_arrests[idx[1], ]), as.numeric(sub_arrests[idx[2], ]))
})
maha_cholesky <- as.matrix(d_maha)[t(pares)]

stopifnot(isTRUE(all.equal(as.numeric(maha_manual), as.numeric(maha_cholesky))))

round(maha_manual, 3)
round(maha_cholesky, 3)


# [dist-almacenamiento]

d_empresas <- dist(empresas)
d_empresas


# [dist-as-matrix]

round(as.matrix(d_empresas), 1)


# [dist-par-mas-cercano]

m_empresas <- as.matrix(d_empresas)
diag(m_empresas) <- Inf

idx_min <- which(m_empresas == min(m_empresas), arr.ind = TRUE)[1, ]
rownames(m_empresas)[idx_min]


# [verificacion-canberra]

r <- as.numeric(sub_arrests[1, ])
s <- as.numeric(sub_arrests[2, ])

canberra_manual <- sum(abs(r - s) / (abs(r) + abs(s)))
canberra_dist    <- as.matrix(dist(sub_arrests, method = "canberra"))[1, 2]

stopifnot(isTRUE(all.equal(canberra_manual, canberra_dist)))

c(manual = canberra_manual, dist = canberra_dist)


## Variables dicotómicas ---------------------------------------------------

# [dicotomicas-datos]

bin_penguins <- penguins |>
  transmute(es_torgersen = as.integer(island == "Torgersen"),
            es_macho     = as.integer(sex == "male")) |>
  slice(1:6)

bin_penguins


# [jaccard-dist]

dist(bin_penguins, method = "binary")


# [verificacion-jaccard]

jaccard_similaridad <- function(x, y) {
  a  <- sum(x == 1 & y == 1)
  bc <- sum(x != y)
  a / (a + bc)
}

m <- as.matrix(bin_penguins)
n <- nrow(m)
sim <- outer(1:n, 1:n, Vectorize(\(i, j) jaccard_similaridad(m[i, ], m[j, ])))

dist_manual <- as.dist(1 - sim)
dist_R      <- dist(bin_penguins, method = "binary")

stopifnot(isTRUE(all.equal(as.numeric(dist_manual), as.numeric(dist_R))))

round(as.matrix(dist_manual), 2)
round(as.matrix(dist_R), 2)


# [sokal-sorensen]

conteo_abcd <- function(x, y) {
  c(a = sum(x == 1 & y == 1), b = sum(x == 1 & y == 0),
    c = sum(x == 0 & y == 1), d = sum(x == 0 & y == 0))
}

similaridad_matriz <- function(datos, formula) {
  m_bin <- as.matrix(datos)
  n_bin <- nrow(m_bin)
  outer(1:n_bin, 1:n_bin, Vectorize(\(i, j) {
    abcd <- conteo_abcd(m_bin[i, ], m_bin[j, ])
    formula(abcd["a"], abcd["b"], abcd["c"], abcd["d"])
  }))
}

sim_jaccard        <- similaridad_matriz(bin_penguins, \(a, b, c, d) a / (a + b + c))
sim_sokal_michener <- similaridad_matriz(bin_penguins, \(a, b, c, d) (a + d) / (a + b + c + d))
sim_sorensen       <- similaridad_matriz(bin_penguins, \(a, b, c, d) (2 * a) / (2 * a + b + c))

round(sim_sokal_michener, 2)
round(sim_sorensen, 2)


# [verificacion-sokal-michener]

stopifnot(isTRUE(all.equal(sim_sokal_michener, sim_jaccard)))

round(sim_jaccard, 2)


# [transformacion-similaridad-distancia]

D2_transformado <- 2 * (1 - sim)
round(D2_transformado[1:4, 1:4], 3)


## Distancia de Gower para datos mixtos ------------------------------------

# [gower]

mix_penguins <- penguins |>
  transmute(bill_length_mm, bill_depth_mm,
            island = factor(island), sex = factor(sex)) |>
  slice(1:6)

mix_penguins

daisy(mix_penguins, metric = "gower")


## Estandarización ---------------------------------------------------------

# [estandarizacion]

rango1   <- \(x) x / (max(x) - min(x))
rango01  <- \(x) (x - min(x)) / (max(x) - min(x))

USArrests |>
  mutate(across(everything(), list(z = scale, rango1 = rango1, rango01 = rango01))) |>
  select(starts_with("Murder")) |>
  head(4)


# Agrupamiento jerárquico --------------------------------------------------


## Paso a paso: el método del centroide ------------------------------------

# [empresas-distancias]

d2_empresas <- dist(empresas)^2
round(as.matrix(d2_empresas))


# [empresas-centroide-e34]

centroide_e34 <- colMeans(empresas[c("E3", "E4"), ])
centroide_e34


# [empresas-hclust-centroide]

hc_centroide <- hclust(d2_empresas, method = "centroid")

# Historial de conglomeración: en qué distancia (height) se unió cada par
tibble(paso = seq_along(hc_centroide$height),
       height = hc_centroide$height,
       fusiona_1 = hc_centroide$merge[, 1],
       fusiona_2 = hc_centroide$merge[, 2])


# [empresas-dendrograma-centroide]
#| fig-width: 6
#| fig-height: 4.5

fviz_dend(hc_centroide, main = "Método del centroide — 8 empresas")


# [verificacion-centroide]

# Paso 1: fusiona las observaciones -3 y -4 (E3, E4), a distancia 13
stopifnot(
  identical(hc_centroide$merge[1, ], c(-3L, -4L)),
  isTRUE(all.equal(hc_centroide$height[1], 13))
)

# Paso 5: fusiona el grupo del paso 1 (E34) con el grupo del paso 3 (E12) --
# hclust() codifica esa referencia como merge = c(1, 3). La distancia en ese
# paso debe coincidir con la distancia al cuadrado entre los centroides ya
# recalculados de ambos grupos.
centroide_e12 <- colMeans(empresas[c("E1", "E2"), ])
d2_e34_e12    <- sum((centroide_e34 - centroide_e12)^2)

stopifnot(
  identical(hc_centroide$merge[5, ], c(1L, 3L)),
  isTRUE(all.equal(d2_e34_e12, hc_centroide$height[5]))
)

c(paso1_merge = hc_centroide$merge[1, ], paso1_height = hc_centroide$height[1])
c(paso5_distancia_recalculada = d2_e34_e12, paso5_height = hc_centroide$height[5])


## Métodos de enlace -------------------------------------------------------

# [empresas-metodos-comparacion]

metodos <- c("single", "complete", "average", "ward.D")

map(metodos, \(m) hclust(d2_empresas, method = m)$height) |>
  set_names(metodos) |>
  as_tibble() |>
  mutate(paso = row_number(), .before = 1)


# [verificacion-ward]

I_formula <- (2 * 2 * 2 / (2 + 2)) * sum((centroide_e12 - centroide_e34)^2)

hc_ward <- hclust(d2_empresas, method = "ward.D")

stopifnot(isTRUE(all.equal(I_formula, hc_ward$height[5])))

c(I_formula = I_formula, height_hclust = hc_ward$height[5])


# [arrests-dendrogramas]
#| fig-width: 8
#| fig-height: 5

d_arrests <- dist(USArrests)

fviz_dend(hclust(d_arrests, method = "single"), main = "Vecino más cercano", cex = 0.5)
fviz_dend(hclust(d_arrests, method = "complete"), main = "Vecino más lejano", cex = 0.5)
fviz_dend(hclust(d_arrests, method = "average"), main = "Promedio", cex = 0.5)
fviz_dend(hclust(d_arrests, method = "ward.D2"), main = "Ward", cex = 0.5)


## Coeficiente de correlación cofenético -----------------------------------

# [cofenetica]

tibble(
  metodo = c("single", "complete", "average", "centroid", "ward.D2"),
  cor_cofenetica = map_dbl(metodo, \(m) cor(d_arrests, cophenetic(hclust(d_arrests, method = m))))
) |>
  arrange(desc(cor_cofenetica))


## El efecto de estandarizar -----------------------------------------------

# [efecto-estandarizar]

d_sin_est <- dist(USArrests)
d_con_est <- dist(scale(USArrests))

grupos_sin_est <- cutree(hclust(d_sin_est, method = "ward.D2"), k = 4)
grupos_con_est <- cutree(hclust(d_con_est, method = "ward.D2"), k = 4)

table(sin_estandarizar = grupos_sin_est, con_estandarizar = grupos_con_est)


# [efecto-estandarizar-ari]

adjustedRandIndex(grupos_sin_est, grupos_con_est)


# [arrests-dendrograma-final]
#| fig-width: 7
#| fig-height: 5

hc_arrests <- hclust(d_con_est, method = "ward.D2")

fviz_dend(hc_arrests, k = 4, rect = TRUE, cex = 0.55,
          main = "USArrests estandarizado — Ward, 4 grupos")


# Selección del número de conglomerados {#sec-nro-conglomerados} -----------

# [funcion-scd]

scd_total <- function(datos, grupos) {
  datos |>
    as.data.frame() |>
    split(grupos) |>
    map_dbl(\(g) if (nrow(g) > 1) sum(diag(var(g))) * (nrow(g) - 1) else 0) |>
    sum()
}

scd_por_k <- tibble(
  k   = 2:10,
  scd = map_dbl(k, \(k) scd_total(scale(USArrests), cutree(hc_arrests, k = k)))
)

scd_por_k


# [scd-codo]
#| fig-width: 7
#| fig-height: 4.5

scd_por_k |>
  ggplot(aes(k, scd)) +
  geom_line() +
  geom_point(size = 2) +
  scale_x_continuous(breaks = scd_por_k$k) +
  labs(title = "Criterio del codo", x = "Número de grupos", y = "Suma de cuadrados dentro")


## F de Hartigan -----------------------------------------------------------

# [hartigan-f]

n <- nrow(USArrests)

hartigan <- tibble(
  de = scd_por_k$k[-nrow(scd_por_k)],
  a  = scd_por_k$k[-1],
  F  = (scd_por_k$scd[-nrow(scd_por_k)] - scd_por_k$scd[-1]) /
       (scd_por_k$scd[-1] / (n - scd_por_k$k[-nrow(scd_por_k)] - 1))
)

hartigan


# [hartigan-conclusion]
#| output: asis

# La regla es secuencial: se parte de g=2 y se suma un grupo mientras F > 10;
# el número de grupos sugerido es el primer g en el que la condición falla
# (no el mayor g con F > 10, que asume sin comprobarlo que F es monótona)
primera_falla <- hartigan |> filter(F <= 10) |> slice(1)
k_hartigan <- if (nrow(primera_falla) > 0) primera_falla$de else max(hartigan$a)

cat(glue::glue(
  "Siguiendo la regla de Hartigan (agregar un grupo mientras F > 10), ",
  "el criterio sugiere **{k_hartigan} grupos**: el paso de {k_hartigan - 1} a {k_hartigan} ",
  "todavía supera el umbral, mientras que el paso de {k_hartigan} a {k_hartigan + 1} ya no lo hace."
))


## Calinski y Harabasz -----------------------------------------------------

# [calinski-harabasz]
#| fig-width: 7
#| fig-height: 4.5

ch <- tibble(
  k  = 2:10,
  CH = map_dbl(k, \(k) calinhara(scale(USArrests), cutree(hc_arrests, k = k)))
)

ch |>
  ggplot(aes(k, CH)) +
  geom_line() +
  geom_point(size = 2) +
  scale_x_continuous(breaks = 2:10) +
  labs(title = "Índice de Calinski-Harabasz", x = "Número de grupos", y = "CH")


# [calinski-harabasz-maximo]

ch |> slice_max(CH, n = 1)


## Coeficiente de silueta --------------------------------------------------

# [silueta]

silueta <- tibble(
  k  = 2:10,
  SC = map_dbl(k, \(k) mean(silhouette(cutree(hc_arrests, k = k), d_con_est)[, 3]))
)

fviz_nbclust(scale(USArrests), FUNcluster = hcut, method = "silhouette", k.max = 10) +
  labs(title = "Coeficiente de silueta promedio por número de grupos")

silueta |> slice_max(SC, n = 1)


# [silueta-individual]
#| fig-width: 7
#| fig-height: 7

grupos4 <- cutree(hc_arrests, k = 4)
sil4    <- silhouette(grupos4, d_con_est)
rownames(sil4) <- rownames(USArrests)

fviz_silhouette(sil4) +
  labs(title = "Silueta individual por estado — 4 grupos")


# [silueta-negativos]

sil4_df <- as_tibble(sil4[, 1:3], rownames = "estado") |>
  arrange(sil_width)

sil4_df |> filter(sil_width < 0)


# [silueta-negativos-conclusion]
#| output: asis

estados_negativos <- sil4_df |> filter(sil_width < 0) |> pull(estado)

lista_estados <- glue::glue_collapse(glue::glue("**{estados_negativos}**"), sep = ", ", last = " y ")
verbo <- if (length(estados_negativos) == 1) "queda" else "quedan"

if (length(estados_negativos) > 0) {
  cat(glue::glue(
    "Casi todas las siluetas son positivas, pero unos pocos estados quedan con ",
    "$s_i$ cercano a cero o negativo: son los que están en el borde entre dos ",
    "grupos, más parecidos al grupo vecino que al propio. En esta partición ",
    "{lista_estados} {verbo} con $s_i<0$: `hclust()` los asigna a un grupo, pero ",
    "su distancia al grupo vecino es, en promedio, menor que la distancia dentro ",
    "del propio grupo. No es un error del algoritmo —cada observación tiene que ",
    "ir a algún grupo— sino una señal de que, para esos estados en particular, ",
    "la frontera entre grupos no está bien definida."
  ))
} else {
  cat(glue::glue(
    "En esta partición todas las siluetas son positivas: no hay estados en el ",
    "borde entre dos grupos, más parecidos al grupo vecino que al propio."
  ))
}


# Descripción de los grupos ------------------------------------------------

# [descripcion-grupos]

media_general <- colMeans(USArrests)
var_general   <- diag(var(USArrests))

describir_grupo <- function(g) {
  datos_g <- USArrests[grupos4 == g, ]
  n_g     <- nrow(datos_g)
  media_g <- colMeans(datos_g)

  error_est <- sqrt((var_general / n_g) * ((n - n_g) / (n - 1)))
  t_stat    <- (media_g - media_general) / error_est
  p_val     <- 2 * pt(abs(t_stat), df = n_g - 1, lower.tail = FALSE)

  tibble(variable = names(t_stat), media_grupo = media_g, t = t_stat, p_valor = p_val)
}

descripcion <- map(1:4, describir_grupo) |>
  set_names(paste0("Grupo ", 1:4, " (n=", table(grupos4), ")"))

descripcion


# K-medias -----------------------------------------------------------------

# [kmeans-basico]

set.seed(123)
km_arrests <- kmeans(scale(USArrests), centers = 4, nstart = 25)

km_arrests$size
km_arrests$centers


# [kmeans-cluster-plot]
#| fig-width: 7
#| fig-height: 5.5

fviz_cluster(km_arrests, data = scale(USArrests), geom = "point",
             ellipse.type = "norm", main = "K-medias — USArrests estandarizado, k=4")


## Estabilidad: por qué importa `nstart` -----------------------------------

# [kmeans-estabilidad]

comparar_nstart <- function(nstart, semilla) {
  set.seed(semilla)
  kmeans(scale(USArrests), centers = 4, nstart = nstart)$tot.withinss
}

tibble(semilla = 1:10) |>
  mutate(
    nstart_1  = map_dbl(semilla, \(s) comparar_nstart(1, s)),
    nstart_25 = map_dbl(semilla, \(s) comparar_nstart(25, s))
  )


## K-medias contra el agrupamiento jerárquico ------------------------------

# [kmeans-vs-jerarquico]

adjustedRandIndex(km_arrests$cluster, grupos4)


# Validación con grupos conocidos ------------------------------------------

# [penguins-preparacion]

penguins_num <- penguins |>
  select(bill_length_mm, bill_depth_mm, flipper_length_mm, body_mass_g)

penguins_num_std <- scale(penguins_num)


# [penguins-jerarquico]

hc_penguins <- hclust(dist(penguins_num_std), method = "ward.D2")
grupos_penguins <- cutree(hc_penguins, k = 3)

table(especie_real = penguins$species, conglomerado = grupos_penguins)


# [penguins-ari]

ari_jerarquico <- adjustedRandIndex(grupos_penguins, penguins$species)
ari_jerarquico


# [penguins-kmedias]

set.seed(123)
km_penguins <- kmeans(penguins_num_std, centers = 3, nstart = 25)

table(especie_real = penguins$species, conglomerado = km_penguins$cluster)

ari_kmedias <- adjustedRandIndex(km_penguins$cluster, penguins$species)
ari_kmedias


# [penguins-conclusion]
#| output: asis

mejor <- if (ari_jerarquico > ari_kmedias) "el agrupamiento jerárquico" else "k-medias"

cat(glue::glue(
  "El índice de Rand ajustado es **{round(ari_jerarquico, 2)}** para el agrupamiento jerárquico ",
  "y **{round(ari_kmedias, 2)}** para k-medias: ambos recuperan la estructura de especies ",
  "razonablemente bien (contra un valor esperado de 0 si la asignación fuera al azar), ",
  "y en este caso {mejor} se ajusta algo mejor a las especies reales. ",
  "La confusión, cuando ocurre, es sobre todo entre Adelie y Chinstrap: son las dos especies ",
  "más parecidas en las variables morfológicas usadas."
))


# [penguins-cluster-plot]
#| fig-width: 7
#| fig-height: 5.5

fviz_cluster(list(data = penguins_num_std, cluster = grupos_penguins), geom = "point",
             ellipse.type = "norm", main = "Conglomerados jerárquicos — penguins") +
  labs(subtitle = "Comparar con la especie real en la tabla de contingencia de arriba")


# Comparación de métodos de agrupamiento -----------------------------------

# [clvalid]

comparacion <- clValid(
  obj        = scale(USArrests),
  nClust     = 2:6,
  clMethods  = c("hierarchical", "kmeans", "pam"),
  validation = c("internal", "stability")
)

summary(comparacion)


# Otros algoritmos ---------------------------------------------------------


## K-medias recortado ------------------------------------------------------

# [rskc]

rskc <- RSKC::RSKC(scale(USArrests), ncl = 4, alpha = 0.1)

table(rskc$labels)
rownames(USArrests)[rskc$oW]   # observaciones tratadas como atípicas y excluidas


## K-modas -----------------------------------------------------------------

# [kmodas]

set.seed(1)
categoricas_penguins <- penguins |>
  transmute(island = factor(island), sex = factor(sex)) |>
  as.data.frame()   # kmodes() indexa con data[, j]: no funciona sobre un tibble
kmodas <- klaR::kmodes(categoricas_penguins, modes = 3)

table(especie = penguins$species, conglomerado = kmodas$cluster)


## K-prototipos ------------------------------------------------------------

# [kprototipos]

mixtas_penguins <- penguins |>
  transmute(bill_length_mm, bill_depth_mm, flipper_length_mm, body_mass_g,
            island = factor(island), sex = factor(sex))

set.seed(1)
kproto <- clustMixType::kproto(mixtas_penguins, k = 3, verbose = FALSE)

table(especie = penguins$species, conglomerado = kproto$cluster)

