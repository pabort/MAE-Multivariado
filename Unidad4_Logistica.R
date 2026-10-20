
# Introducción -------------------------------------------------------------


# Carga de paquetes --------------------------------------------------------

# [setup]

for (loc in c("es_AR.UTF-8", "es_ES.UTF-8", "en_US.UTF-8", "C.UTF-8", "C.utf8")) {
  if (Sys.setlocale("LC_CTYPE", loc) != "") break
}

library(ISLR2)              # base Default
library(MASS)               # lda(), qda()
library(pROC)               # curva ROC, AUC, punto de corte
library(caret)              # confusionMatrix()
library(ResourceSelection)  # hoslem.test()
library(broom)              # tidy(): coeficientes como tabla
library(knitr)              # kable()
library(tidyverse)          # al final, para que dplyr::select() no quede enmascarado

theme_set(theme_minimal(base_size = 11))

# Formato de números para el texto: coma decimal y punto de miles
num <- function(x, d = 3) formatC(x, format = "f", digits = d, decimal.mark = ",", big.mark = ".")
pval <- function(p) if (p < 1e-4) "menor que 0,0001" else num(p, 4)


# De clasificar a estimar una probabilidad ---------------------------------


## Por qué no alcanza con una recta ----------------------------------------

# [lineal-vs-logistica]
#| fig-width: 7
#| fig-height: 4

default_01 <- Default |> mutate(y = as.numeric(default == "Yes"))

ajuste_lineal    <- lm(y ~ balance, data = default_01)
ajuste_logistico <- glm(y ~ balance, family = binomial, data = default_01)

grilla_balance <- tibble(balance = seq(0, max(Default$balance), length.out = 300))

curvas <- grilla_balance |>
  mutate(Lineal    = predict(ajuste_lineal, grilla_balance),
         Logística = predict(ajuste_logistico, grilla_balance, type = "response")) |>
  pivot_longer(-balance, names_to = "modelo", values_to = "prob")

ggplot(default_01, aes(balance, y)) +
  geom_jitter(height = 0.02, width = 0, alpha = 0.08, size = 0.6) +
  geom_line(data = curvas, aes(y = prob, colour = modelo), linewidth = 1) +
  geom_hline(yintercept = c(0, 1), linetype = "dashed", colour = "grey50") +
  labs(title = "Recta contra curva logística",
       x = "Saldo de la tarjeta (balance)", y = "P(default = Yes)", colour = NULL)

predicciones_negativas <- sum(fitted(ajuste_lineal) < 0)


## El modelo logístico -----------------------------------------------------


## Probabilidad, odds y logit: tres escalas {#sec-escalas} -----------------

# [tres-escalas]

tibble(probabilidad = c(0.01, 0.10, 0.20, 0.50, 0.80, 0.90, 0.99)) |>
  mutate(odds  = probabilidad / (1 - probabilidad),
         logit = log(odds)) |>
  kable(digits = 3)


# Primer ejemplo: `Default` {#sec-default} ---------------------------------

# [default-exploracion]

glimpse(Default)

Default |>
  count(default) |>
  mutate(proporcion = n / sum(n))


## Modelo simple con `balance` ---------------------------------------------

# [default-balance]

modelo_balance <- glm(default ~ balance, family = binomial, data = Default)
summary(modelo_balance)

b_balance <- coef(modelo_balance)["balance"]


# [default-balance-prob]

saldos <- tibble(balance = c(1000, 1500, 2000))

saldos |>
  mutate(logit        = predict(modelo_balance, saldos, type = "link"),
         odds         = exp(logit),
         probabilidad = predict(modelo_balance, saldos, type = "response")) |>
  kable(digits = 4)

prob_saldos <- predict(modelo_balance, saldos, type = "response")


## Modelo simple con `student` ---------------------------------------------

# [default-student]

modelo_student <- glm(default ~ student, family = binomial, data = Default)
summary(modelo_student)

b_student_simple <- coef(modelo_student)["studentYes"]


# [default-student-cociente]

tasa_por_grupo <- tapply(Default$default == "Yes", Default$student, mean)
cociente_tasas_default <- tasa_por_grupo[["Yes"]] / tasa_por_grupo[["No"]]

c(odds_ratio = unname(exp(b_student_simple)), cociente_de_probabilidades = cociente_tasas_default)


## Modelo múltiple: el coeficiente de `student` cambia de signo ------------

# [default-multiple]

modelo_multiple <- glm(default ~ balance + income + student, family = binomial, data = Default)
summary(modelo_multiple)

b_student_multiple <- coef(modelo_multiple)["studentYes"]

stopifnot(sign(b_student_simple) != sign(b_student_multiple))

c(simple = unname(b_student_simple), multiple = unname(b_student_multiple))


# [default-confusion-figura]
#| fig-width: 9
#| fig-height: 4
#| layout-ncol: 2

grilla_student <- expand_grid(balance = seq(0, 2650, by = 10),
                              student = c("No", "Yes"),
                              income  = mean(Default$income)) |>
  mutate(prob = predict(modelo_multiple, pick(everything()), type = "response"))

tasas_globales <- Default |>
  group_by(student) |>
  summarise(tasa = mean(default == "Yes"))

ggplot(grilla_student, aes(balance, prob, colour = student)) +
  geom_line(linewidth = 1) +
  geom_hline(data = tasas_globales, aes(yintercept = tasa, colour = student),
             linetype = "dashed") +
  labs(title = "A igual saldo, los estudiantes fallan menos",
       x = "Saldo (balance)", y = "P(default = Yes)", colour = "Estudiante")

ggplot(Default, aes(student, balance, fill = student)) +
  geom_boxplot(show.legend = FALSE) +
  labs(title = "Los estudiantes tienen saldos más altos",
       x = "Estudiante", y = "Saldo (balance)")

saldo_medio <- tapply(Default$balance, Default$student, mean)


## Lo que deja planteado: el punto de corte --------------------------------

# [default-tabla]

prob_default <- fitted(modelo_multiple)
clase_default <- factor(ifelse(prob_default > 0.5, "Yes", "No"), levels = c("Yes", "No"))

tabla_default <- table(Real = factor(Default$default, levels = c("Yes", "No")),
                       Predicho = clase_default)
tabla_default

sens_default <- tabla_default["Yes", "Yes"] / sum(tabla_default["Yes", ])
exac_default <- sum(diag(tabla_default)) / sum(tabla_default)


# Estimación por máxima verosimilitud {#sec-estimacion} --------------------


## La verosimilitud --------------------------------------------------------


## Newton-Raphson, paso a paso ---------------------------------------------

# [irls-manual]

X <- model.matrix(modelo_multiple)
y <- as.numeric(Default$default == "Yes")

beta <- rep(0, ncol(X))   # valor inicial: todos los coeficientes en cero
historial <- tibble()

for (iteracion in 1:25) {
  eta   <- as.vector(X %*% beta)
  p_hat <- plogis(eta)
  pesos <- p_hat * (1 - p_hat)
  z     <- eta + (y - p_hat) / pesos
  beta_nuevo <- solve(crossprod(X, pesos * X), crossprod(X, pesos * z))

  log_verosimilitud <- sum(y * log(p_hat) + (1 - y) * log(1 - p_hat))
  cambio <- max(abs(beta_nuevo - beta))
  historial <- bind_rows(historial,
                         tibble(iteracion, log_verosimilitud, cambio_max = cambio))
  beta <- as.vector(beta_nuevo)
  if (cambio < 1e-10) break
}

historial


# [irls-verificacion]

stopifnot(isTRUE(all.equal(beta, unname(coef(modelo_multiple)))))

cbind(IRLS = beta, glm = coef(modelo_multiple))


## La evidencia en la salida de `glm()` ------------------------------------

# [irls-convergencia]

c(iteraciones = modelo_multiple$iter, convergio = modelo_multiple$converged)


# La aplicación: deserción de clientela ------------------------------------


## El problema -------------------------------------------------------------


## Diccionario de variables ------------------------------------------------


## Lectura de los datos y niveles de referencia ----------------------------

# [desercion-lectura]

url_datos <- "https://raw.githubusercontent.com/pabort/MAE-Multivariado/refs/heads/main/desercion.csv"
datos <- read.csv(url_datos, stringsAsFactors = FALSE, encoding = "UTF-8")

# El orden de los niveles no es alfabético: fija el nivel de referencia
# contra el que se leen todos los coeficientes.
niveles <- list(
  TPROM4_2 = c("De 0,02 a 1,2487", "Más de 1,2487 a 500", "> 500"),
  MMODELO  = c("Motorola", "Nokia", "Otros"),
  TDEUDA   = c("Sin deuda", "hasta 157,25", "> a 157,25"),
  REGION_3 = c("Bs As, Sta Fe y La Pampa", "Resto", "Patagonia"),
  ADQUISI1 = c("Comodato", "Dist. y Leasing", "Propio"),
  REST_AG1 = c("Activo", "Baja", "Tramita baja y OE"),
  RTIP_CTA = c("Negocios", "Otros", "Personal", "Top"),
  RESTADO  = c("ACTIVO", "CANCELADO")
)
for (v in names(niveles)) datos[[v]] <- factor(datos[[v]], levels = niveles[[v]])

stopifnot(
  ncol(datos) == 14, nrow(datos) == 9702,
  levels(datos$TPROM4_2)[1] == "De 0,02 a 1,2487",
  levels(datos$TDEUDA)[1] == "Sin deuda",
  !anyNA(datos)
)

dim(datos)
str(datos)


## Exploración por grupo ---------------------------------------------------

# [desercion-grupos]

datos |>
  count(RESTADO) |>
  mutate(proporcion = n / sum(n))


# [desercion-numericas]

datos |>
  select(RESTADO, where(is.numeric)) |>
  pivot_longer(-RESTADO, names_to = "variable", values_to = "valor") |>
  group_by(variable, RESTADO) |>
  summarise(media = mean(valor), mediana = median(valor), desvio = sd(valor),
            .groups = "drop") |>
  kable(digits = 2)


# [desercion-categoricas]
#| fig-width: 9
#| fig-height: 7

datos |>
  select(where(is.factor)) |>
  pivot_longer(-RESTADO, names_to = "variable", values_to = "nivel") |>
  group_by(variable, nivel) |>
  summarise(tasa_cancelacion = mean(RESTADO == "CANCELADO"), n = n(), .groups = "drop") |>
  ggplot(aes(tasa_cancelacion, nivel)) +
  geom_col(fill = "steelblue") +
  geom_vline(xintercept = mean(datos$RESTADO == "CANCELADO"), linetype = "dashed") +
  facet_wrap(~ variable, scales = "free_y", ncol = 2) +
  labs(title = "Tasa de cancelación por nivel de cada factor",
       subtitle = "La línea punteada es la tasa global",
       x = "Proporción de clientes que cancelan", y = NULL)


## Entrenamiento y test ----------------------------------------------------

# [particion]

train <- seq_len(5771)
entrenamiento <- datos[train, ]
test          <- datos[-train, ]

composicion <- bind_rows(
  entrenamiento |> count(RESTADO) |> mutate(muestra = "entrenamiento"),
  test          |> count(RESTADO) |> mutate(muestra = "test")
) |>
  group_by(muestra) |>
  mutate(proporcion = n / sum(n)) |>
  ungroup()

composicion

stopifnot(
  nrow(entrenamiento) == 5771, nrow(test) == 3931,
  identical(as.vector(table(entrenamiento$RESTADO)), c(3421L, 2350L)),
  identical(as.vector(table(test$RESTADO)), c(2291L, 1640L))
)

prop_canc <- composicion |> filter(RESTADO == "CANCELADO") |> pull(proporcion, name = muestra)


# [particion-alternativa]
#| eval: false

# set.seed(2026)
# train_aleatorio <- sample(nrow(datos), size = round(0.6 * nrow(datos)))
# 
# set.seed(2026)
# train_estratificado <- createDataPartition(datos$RESTADO, p = 0.6, list = FALSE)


# El modelo y cómo se lee su salida {#sec-modelo} --------------------------


## Ajuste ------------------------------------------------------------------

# [modelo-ajuste]

modelo_desercion <- glm(RESTADO ~ ., family = binomial, data = entrenamiento)
summary(modelo_desercion)

n_ent <- nrow(entrenamiento)
k_coef <- length(coef(modelo_desercion))


# [modelo-x-parcial]

"x" %in% names(modelo_desercion)

stopifnot(identical(modelo_desercion$x, modelo_desercion$xlevels))

modelo_desercion$x


## Los coeficientes --------------------------------------------------------


## El test de Wald {#sec-wald} ---------------------------------------------

# [wald-manual]

coefs <- tidy(modelo_desercion)

wald_ejemplo <- coefs |>
  filter(term %in% c("D_EDAD", "MMODELOOtros")) |>
  mutate(W = estimate / std.error,
         p_normal = 2 * pnorm(-abs(W)),
         p_chi2   = pchisq(W^2, df = 1, lower.tail = FALSE))

stopifnot(
  isTRUE(all.equal(wald_ejemplo$W, wald_ejemplo$statistic)),
  isTRUE(all.equal(wald_ejemplo$p_normal, wald_ejemplo$p.value)),
  isTRUE(all.equal(wald_ejemplo$p_chi2, wald_ejemplo$p.value))
)

wald_ejemplo |> select(term, estimate, std.error, statistic, W, p.value, p_normal, p_chi2)

w_edad <- wald_ejemplo |> filter(term == "D_EDAD")
w_otros <- wald_ejemplo |> filter(term == "MMODELOOtros")


## La devianza {#sec-devianza} ---------------------------------------------

# [devianzas]

D_nula     <- modelo_desercion$null.deviance
D_residual <- deviance(modelo_desercion)
gl_nula     <- modelo_desercion$df.null
gl_residual <- df.residual(modelo_desercion)

modelo_nulo <- glm(RESTADO ~ 1, family = binomial, data = entrenamiento)

stopifnot(
  round(D_nula, 1) == 7800.4, gl_nula == 5770,
  round(D_residual, 1) == 3555.7, gl_residual == 5749,
  isTRUE(all.equal(D_residual, -2 * as.numeric(logLik(modelo_desercion)))),
  isTRUE(all.equal(D_nula, deviance(modelo_nulo)))
)

tibble(modelo = c("nulo (solo constante)", "completo"),
       devianza = c(D_nula, D_residual),
       gl = c(gl_nula, gl_residual),
       menos_2_logL = -2 * c(logLik(modelo_nulo), logLik(modelo_desercion)))


## El contraste global -----------------------------------------------------

# [contraste-global]

G_global  <- D_nula - D_residual
gl_global <- gl_nula - gl_residual
p_global  <- pchisq(G_global, df = gl_global, lower.tail = FALSE)

stopifnot(round(G_global, 1) == 4244.7, gl_global == 21)

c(G = G_global, gl = gl_global, valor_critico_5pc = qchisq(0.95, gl_global), p_valor = p_global)


## El AIC ------------------------------------------------------------------

# [aic]

stopifnot(
  round(AIC(modelo_desercion), 1) == 3599.7,
  isTRUE(all.equal(AIC(modelo_desercion), D_residual + 2 * k_coef))
)

c(AIC = AIC(modelo_desercion), D_mas_2k = D_residual + 2 * k_coef)


# Odds ratios e interpretación {#sec-odds} ---------------------------------


## Odds ratios con intervalos de confianza ---------------------------------

# [odds-ratios]

tabla_or <- coefs |>
  filter(term != "(Intercept)") |>
  mutate(odds_ratio = exp(estimate),
         ic_inf = exp(estimate - qnorm(0.975) * std.error),
         ic_sup = exp(estimate + qnorm(0.975) * std.error),
         contiene_1 = ic_inf < 1 & ic_sup > 1) |>
  select(term, estimate, odds_ratio, ic_inf, ic_sup, p.value, contiene_1) |>
  arrange(desc(odds_ratio))

ic_wald_r <- exp(confint.default(modelo_desercion))[tabla_or$term, ]

stopifnot(
  isTRUE(all.equal(tabla_or$odds_ratio, unname(exp(coef(modelo_desercion))[tabla_or$term]))),
  isTRUE(all.equal(cbind(tabla_or$ic_inf, tabla_or$ic_sup), unname(ic_wald_r))),
  # el intervalo contiene al 1 exactamente cuando Wald no rechaza al 5 %
  identical(tabla_or$contiene_1, tabla_or$p.value > 0.05)
)

tabla_or |> kable(digits = 4)


# [odds-contienen-1]
#| output: asis

contienen <- tabla_or |> filter(contiene_1)
otros <- tabla_or |> filter(term == "MMODELOOtros")

cat(glue::glue(
  "Hay {nrow(contienen)} odds ratios cuyo intervalo contiene al 1: ",
  "{paste0('`', contienen$term, '`', collapse = ', ')}. ",
  "Por ejemplo, `MMODELOOtros` tiene odds ratio {num(otros$odds_ratio, 3)} con intervalo ",
  "({num(otros$ic_inf, 3)}; {num(otros$ic_sup, 3)}) y p-valor {num(otros$p.value, 2)}: el intervalo ",
  "va desde una reducción de los odds de alrededor del {num(100 * (1 - otros$ic_inf), 0)} % hasta un ",
  "aumento del {num(100 * (otros$ic_sup - 1), 0)} %, así que es compatible con que el modelo del equipo ",
  "\"Otros\" no difiera de Motorola."
))


# [odds-grafico]
#| fig-width: 8
#| fig-height: 6

ggplot(tabla_or, aes(odds_ratio, fct_reorder(term, odds_ratio), colour = contiene_1)) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_pointrange(aes(xmin = ic_inf, xmax = ic_sup)) +
  scale_x_log10() +
  scale_colour_manual(values = c(`FALSE` = "steelblue", `TRUE` = "firebrick"),
                      labels = c("No contiene al 1", "Contiene al 1")) +
  labs(title = "Odds ratios e intervalos de confianza del 95 %",
       subtitle = "Eje horizontal en escala logarítmica",
       x = "Odds ratio", y = NULL, colour = NULL)


# [odds-publicados]

or_publicados <- c(
  PROM_ACT = 1.0005, D_EDAD = 0.9933, D_FINCON = 0.9987, PORC_DES = 1.3734,
  HIST_P = 1.9486, D_DURAC = 0.4805,
  "TPROM4_2Más de 1,2487 a 500" = 1.7713, "TPROM4_2> 500" = 6.4174,
  MMODELONokia = 1.4847, MMODELOOtros = 0.9056,
  "TDEUDAhasta 157,25" = 0.5084, "TDEUDA> a 157,25" = 3.0611,
  REGION_3Resto = 1.2222, REGION_3Patagonia = 1.6303,
  "ADQUISI1Dist. y Leasing" = 0.7871, ADQUISI1Propio = 1.6280,
  REST_AG1Baja = 2.1847, "REST_AG1Tramita baja y OE" = 1.4195,
  RTIP_CTAOtros = 0.4950, RTIP_CTAPersonal = 0.7700, RTIP_CTATop = 0.1225
)

comparacion_or <- tibble(term = names(or_publicados),
                         publicado = or_publicados,
                         modelo = exp(coef(modelo_desercion))[names(or_publicados)])

stopifnot(all(abs(comparacion_or$modelo - comparacion_or$publicado) < 1e-3))

comparacion_or |> kable(digits = 5)


## Cómo se lee cada tipo de predictora -------------------------------------

# [or-a-probabilidad]

# Probabilidad -> odds -> odds multiplicados por el odds ratio -> probabilidad
prob_tras_or <- function(p_inicial, odds_ratio) {
  odds_final <- p_inicial / (1 - p_inicial) * odds_ratio
  odds_final / (1 + odds_final)
}

tasa_global <- mean(datos$RESTADO == "CANCELADO")
or_histp <- exp(coef(modelo_desercion)[["HIST_P"]])
or_top   <- exp(coef(modelo_desercion)[["RTIP_CTATop"]])
or_tprom <- exp(coef(modelo_desercion)[["TPROM4_2> 500"]])

tibble(termino = c("HIST_P", "RTIP_CTATop", "TPROM4_2> 500"),
       odds_ratio = c(or_histp, or_top, or_tprom),
       p_inicial = tasa_global,
       p_final = prob_tras_or(tasa_global, odds_ratio),
       cociente_probabilidades = p_final / p_inicial) |>
  kable(digits = 3)


## Un odds ratio no multiplica la probabilidad {#sec-or-prob} --------------


### El recorrido con un cliente concreto -----------------------------------

# [or-recorrido]

prob_ajustadas <- fitted(modelo_desercion)
candidatos <- which(entrenamiento$HIST_P < 7)
fila_cliente <- candidatos[which.min(abs(prob_ajustadas[candidatos] - median(prob_ajustadas)))]
cliente <- entrenamiento[fila_cliente, ]

p_cliente     <- prob_ajustadas[[fila_cliente]]
odds_cliente  <- p_cliente / (1 - p_cliente)
odds_nuevos   <- odds_cliente * or_histp
p_nueva       <- odds_nuevos / (1 + odds_nuevos)

# El mismo resultado sale de predecir con HIST_P un punto más alto
cliente_mas_1 <- cliente |> mutate(HIST_P = HIST_P + 1)
p_predicha_mas_1 <- predict(modelo_desercion, newdata = cliente_mas_1, type = "response")

stopifnot(isTRUE(all.equal(p_nueva, unname(p_predicha_mas_1))))

tibble(paso = c("1. probabilidad estimada", "2. odds = p / (1 - p)",
                "3. odds x odds ratio de HIST_P", "4. probabilidad = odds / (1 + odds)",
                "   control: predict() con HIST_P + 1"),
       valor = c(p_cliente, odds_cliente, odds_nuevos, p_nueva, p_predicha_mas_1)) |>
  kable(digits = 4)


### El mismo odds ratio desde distintos puntos de partida ------------------

# [or-tabla-partidas]

tabla_partidas <- tibble(p_inicial = c(0.02, 0.10, 0.20, tasa_global, 0.80, 0.90)) |>
  mutate(p_final = prob_tras_or(p_inicial, or_histp),
         cambio_puntos = 100 * (p_final - p_inicial),
         cociente_probabilidades = p_final / p_inicial)

tabla_partidas |> kable(digits = c(4, 4, 2, 3))

menor_cambio <- min(tabla_partidas$cambio_puntos)
mayor_cambio_pp <- max(tabla_partidas$cambio_puntos)


### Cuándo el odds ratio se parece al cociente de probabilidades -----------

# [or-vs-cociente]

tabla_partidas |>
  mutate(odds_ratio = or_histp,
         or_sobre_cociente = odds_ratio / cociente_probabilidades) |>
  select(p_inicial, cociente_probabilidades, odds_ratio, or_sobre_cociente) |>
  kable(digits = 3)

exceso_raro   <- or_histp / tabla_partidas$cociente_probabilidades[1] - 1
exceso_global <- or_histp / tabla_partidas$cociente_probabilidades[4] - 1


## La escala importa -------------------------------------------------------

# [odds-escala]

desvios <- entrenamiento |>
  select(where(is.numeric)) |>
  pivot_longer(everything(), names_to = "term", values_to = "valor") |>
  group_by(term) |>
  summarise(desvio = sd(valor))

tabla_escala <- coefs |>
  inner_join(desvios, by = "term") |>
  mutate(or_por_unidad = exp(estimate),
         or_por_desvio = exp(estimate * desvio)) |>
  select(term, estimate, desvio, or_por_unidad, or_por_desvio) |>
  arrange(desc(abs(log(or_por_desvio))))

tabla_escala |> kable(digits = 4)

or_prom_100 <- exp(100 * coef(modelo_desercion)["PROM_ACT"])


# [ic-perfil]

ic_perfil <- exp(confint(modelo_desercion))[tabla_or$term, ]

tibble(term = tabla_or$term,
       wald_inf = tabla_or$ic_inf, perfil_inf = ic_perfil[, 1],
       wald_sup = tabla_or$ic_sup, perfil_sup = ic_perfil[, 2]) |>
  kable(digits = 4)


# Razón de verosimilitudes contra Wald {#sec-lrt} --------------------------


## El estadístico $G$ ------------------------------------------------------


## El caso que lo justifica: `MMODELO` -------------------------------------

# [mmodelo-wald]

coefs |> filter(str_starts(term, "MMODELO"))


# [mmodelo-lrt]

modelo_sin_mmodelo <- update(modelo_desercion, . ~ . - MMODELO)

G_mmodelo  <- deviance(modelo_sin_mmodelo) - deviance(modelo_desercion)
gl_mmodelo <- df.residual(modelo_sin_mmodelo) - df.residual(modelo_desercion)
p_mmodelo  <- pchisq(G_mmodelo, df = gl_mmodelo, lower.tail = FALSE)

c(D_sin = deviance(modelo_sin_mmodelo), D_con = deviance(modelo_desercion),
  G = G_mmodelo, gl = gl_mmodelo, p_valor = p_mmodelo)


## `drop1()`: el test para cada variable -----------------------------------

# [drop1]

tabla_drop1 <- drop1(modelo_desercion, test = "LRT")
tabla_drop1

stopifnot(
  isTRUE(all.equal(tabla_drop1["MMODELO", "LRT"], G_mmodelo)),
  isTRUE(all.equal(tabla_drop1["MMODELO", "Pr(>Chi)"], p_mmodelo))
)


# [drop1-resumen]
#| output: asis

no_signif <- rownames(tabla_drop1)[-1][tabla_drop1$`Pr(>Chi)`[-1] > 0.05]
p_max <- max(tabla_drop1$`Pr(>Chi)`[-1])
var_max <- rownames(tabla_drop1)[-1][which.max(tabla_drop1$`Pr(>Chi)`[-1])]

if (length(no_signif) == 0) {
  cat(glue::glue(
    "Las trece predictoras son significativas al 5 % por razón de verosimilitudes: el p-valor más alto ",
    "es el de `{var_max}`, {num(p_max, 4)}. Ninguna variable puede sacarse del modelo sin un aumento ",
    "significativo de la devianza."
  ))
} else {
  cat(glue::glue(
    "Las variables con p-valor mayor que 0,05 son {paste0('`', no_signif, '`', collapse = ', ')}: ",
    "son candidatas a salir del modelo."
  ))
}


# [anova-secuencial]

tabla_anova <- anova(modelo_desercion, test = "LRT")
tabla_anova

stopifnot(
  isTRUE(all.equal(tabla_anova["RTIP_CTA", "Deviance"], tabla_drop1["RTIP_CTA", "LRT"])),
  isTRUE(all.equal(sum(tabla_anova$Deviance, na.rm = TRUE), G_global))
)


## ¿Cuándo preferir uno u otro? --------------------------------------------

# [wald-vs-lrt]

numericas <- c("PROM_ACT", "D_EDAD", "D_FINCON", "PORC_DES", "HIST_P", "D_DURAC")

tibble(variable = numericas,
       W2_wald  = (coefs$statistic[match(numericas, coefs$term)])^2,
       G_lrt    = tabla_drop1[numericas, "LRT"],
       p_wald   = coefs$p.value[match(numericas, coefs$term)],
       p_lrt    = tabla_drop1[numericas, "Pr(>Chi)"]) |>
  kable(digits = 4)


# [wald-vs-lrt-texto]
#| output: asis

cociente <- (coefs$statistic[match(numericas, coefs$term)])^2 / tabla_drop1[numericas, "LRT"]
names(cociente) <- numericas
parecidos <- names(cociente)[cociente > 0.9]
menores   <- names(cociente)[cociente <= 0.9]

cat(glue::glue(
  "$W^2$ y $G$ tienen la misma distribución de referencia, $\\chi^2_{{(1)}}$. ",
  "Para {paste0('`', parecidos, '`', collapse = ', ')} los dos estadísticos son muy parecidos. ",
  "{if (length(menores) > 0) paste0('Para ', paste0('`', menores, '`', collapse = ', '), ', en cambio, $W^2$ es claramente menor que $G$ (hasta un ', num(100 * (1 - min(cociente)), 0), ' % menor). ', if (max(tabla_drop1[menores, 'Pr(>Chi)'], coefs$p.value[match(menores, coefs$term)]) < 1e-4) 'Las conclusiones no cambian, porque con cualquiera de los dos tests esos p-valores son ínfimos, ' else 'Acá las diferencias sí pueden cambiar conclusiones, ', 'pero se ve que, cuando discrepan, Wald es el más conservador de los dos.') else ''}"
))


## Una decisión de modelado: `D_DURAC` numérica o factor {#sec-ddurac} -----

# [ddurac-factor]

table(entrenamiento$D_DURAC)

modelo_durac_factor <- update(modelo_desercion, . ~ . - D_DURAC + factor(D_DURAC))
lrt_durac <- anova(modelo_desercion, modelo_durac_factor, test = "LRT")
lrt_durac

p_durac <- lrt_durac$`Pr(>Chi)`[2]


# [ddurac-texto]
#| output: asis

cat(if (p_durac > 0.05) {
  "No se rechaza $H_0$: tratar a `D_DURAC` como factor no mejora significativamente el ajuste. La tendencia lineal es una simplificación aceptable y el modelo la mantiene como numérica, con un parámetro menos."
} else {
  "Se rechaza $H_0$: la tendencia lineal no alcanza, y conviene tratar a `D_DURAC` como factor."
})


# [histp-factor]

modelo_histp_factor <- update(modelo_desercion, . ~ . - HIST_P + factor(HIST_P))
lrt_histp <- anova(modelo_desercion, modelo_histp_factor, test = "LRT")
lrt_histp

p_histp <- lrt_histp$`Pr(>Chi)`[2]
aic_histp <- c(lineal = AIC(modelo_desercion), factor = AIC(modelo_histp_factor))
aic_histp


# [histp-texto]
#| output: asis

cat(glue::glue(
  "Para `HIST_P`, $G=$ {num(lrt_histp$Deviance[2], 2)} con {lrt_histp$Df[2]} grados de libertad y ",
  "p-valor {pval(p_histp)}. ",
  if (p_histp < 0.05) {
    paste0("Acá **sí** se rechaza la linealidad: los efectos de los siete valores de la historia de pago ",
           "no caen sobre una recta en la escala del logit. El AIC apunta en la misma dirección (",
           num(aic_histp[["factor"]], 1), " con `HIST_P` como factor, contra ", num(aic_histp[["lineal"]], 1),
           " como numérica). En lo que sigue el modelo mantiene a `HIST_P` como numérica, para conservar ",
           "la especificación de la que salen los odds ratios de la sección anterior; pero es una decisión ",
           "con costo, y un análisis definitivo debería considerar la versión con `HIST_P` como factor.")
  } else {
    "Tampoco se rechaza la linealidad: `HIST_P` puede quedar como numérica."
  }
))


# Bondad de ajuste y diagnóstico -------------------------------------------


## Residuos {#sec-residuos} ------------------------------------------------

# [residuos]

y_ent <- as.numeric(entrenamiento$RESTADO == "CANCELADO")
p_ent <- fitted(modelo_desercion)

res_pearson  <- (y_ent - p_ent) / sqrt(p_ent * (1 - p_ent))
res_devianza <- sign(y_ent - p_ent) * sqrt(-2 * (y_ent * log(p_ent) + (1 - y_ent) * log(1 - p_ent)))

stopifnot(
  isTRUE(all.equal(res_pearson, residuals(modelo_desercion, type = "pearson"))),
  isTRUE(all.equal(res_devianza, residuals(modelo_desercion, type = "deviance"))),
  isTRUE(all.equal(sum(res_devianza^2), D_residual))
)

c(chi2_pearson = sum(res_pearson^2), D = sum(res_devianza^2), gl = gl_residual)


# [residuos-grafico]
#| fig-width: 8
#| fig-height: 4

tibble(p_hat = p_ent, d = res_devianza, RESTADO = entrenamiento$RESTADO) |>
  ggplot(aes(p_hat, d, colour = RESTADO)) +
  geom_point(alpha = 0.3, size = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "Residuos de devianza contra probabilidad estimada",
       x = expression(hat(p)[i]), y = expression(d[i]), colour = NULL)


## Test de Hosmer-Lemeshow -------------------------------------------------

# [hosmer-lemeshow]

hl <- hoslem.test(y_ent, p_ent, g = 10)
hl

tabla_hl <- tibble(grupo = rownames(hl$observed),
                   cancelados_obs = hl$observed[, "y1"], cancelados_esp = hl$expected[, "yhat1"],
                   activos_obs = hl$observed[, "y0"], activos_esp = hl$expected[, "yhat0"])

stopifnot(isTRUE(all.equal(sum((hl$observed - hl$expected)^2 / hl$expected),
                           unname(hl$statistic))))

tabla_hl |> kable(digits = 1)


# [hosmer-lemeshow-texto]
#| output: asis

aportes <- tabla_hl |>
  mutate(aporte = (cancelados_obs - cancelados_esp)^2 / cancelados_esp +
                  (activos_obs - activos_esp)^2 / activos_esp,
         proporcion = aporte / sum(aporte)) |>
  arrange(desc(aporte))

cat(glue::glue(
  "$\\hat C=$ {num(hl$statistic, 2)} con {hl$parameter} grados de libertad, p-valor {num(hl$p.value, 4)}. ",
  "{if (hl$p.value < 0.05) 'Al 5 % se rechaza $H_0$.' else 'Al 5 % no se rechaza $H_0$.'} ",
  "Calculando el aporte de cada grupo al estadístico se ve de dónde sale: los dos grupos que más aportan ",
  "son los de probabilidades {aportes$grupo[1]} ({aportes$cancelados_obs[1]} cancelaciones observadas ",
  "contra {num(aportes$cancelados_esp[1], 1)} esperadas) y {aportes$grupo[2]} ({aportes$cancelados_obs[2]} ",
  "contra {num(aportes$cancelados_esp[2], 1)}), que juntos explican el ",
  "{num(100 * sum(aportes$proporcion[1:2]), 0)} % de $\\hat C$. En los otros {nrow(aportes) - 2} grupos ",
  "las diferencias entre observado y esperado son de a lo sumo ",
  "{num(max(abs(aportes$cancelados_obs[-(1:2)] - aportes$cancelados_esp[-(1:2)])), 1)} cancelaciones."
))


## Pseudo $R^2$ ------------------------------------------------------------

# [pseudo-r2]

r2_mcfadden   <- 1 - as.numeric(logLik(modelo_desercion)) / as.numeric(logLik(modelo_nulo))
r2_cox_snell  <- 1 - exp((D_residual - D_nula) / n_ent)
r2_nagelkerke <- r2_cox_snell / (1 - exp(-D_nula / n_ent))

stopifnot(
  isTRUE(all.equal(r2_mcfadden, 1 - D_residual / D_nula)),
  round(r2_mcfadden, 3) == 0.544
)

c(McFadden = r2_mcfadden, Cox_Snell = r2_cox_snell, Nagelkerke = r2_nagelkerke)


## Multicolinealidad -------------------------------------------------------

# [gvif]

R_coef <- cov2cor(vcov(modelo_desercion))[-1, -1]
variable_de_coef <- attr(model.matrix(modelo_desercion), "assign")[-1]
nombres_variables <- attr(terms(modelo_desercion), "term.labels")

tabla_gvif <- map_dfr(seq_along(nombres_variables), function(j) {
  bloque <- variable_de_coef == j
  gvif <- det(R_coef[bloque, bloque, drop = FALSE]) * det(R_coef[!bloque, !bloque]) / det(R_coef)
  tibble(variable = nombres_variables[j], gl = sum(bloque), GVIF = gvif,
         GVIF_ajustado = gvif^(1 / (2 * sum(bloque))))
}) |>
  arrange(desc(GVIF_ajustado))

tabla_gvif |> kable(digits = 3)


## Observaciones influyentes -----------------------------------------------

# [influencia]

leverage <- hatvalues(modelo_desercion)
cook     <- cooks.distance(modelo_desercion)

tibble(medida = c("leverage", "distancia de Cook"),
       umbral = c(2 * k_coef / n_ent, 4 / n_ent),
       maximo = c(max(leverage), max(cook)),
       casos_sobre_umbral = c(sum(leverage > 2 * k_coef / n_ent), sum(cook > 4 / n_ent)))


# [influencia-grafico]
#| fig-width: 8
#| fig-height: 4

tibble(leverage, cook, RESTADO = entrenamiento$RESTADO) |>
  ggplot(aes(leverage, cook, colour = RESTADO)) +
  geom_point(alpha = 0.4, size = 0.8) +
  geom_hline(yintercept = 4 / n_ent, linetype = "dashed") +
  geom_vline(xintercept = 2 * k_coef / n_ent, linetype = "dashed") +
  labs(title = "Leverage y distancia de Cook",
       subtitle = "Líneas punteadas: 4/n y 2(p+1)/n",
       x = "Leverage", y = "Distancia de Cook", colour = NULL)


# [influencia-reajuste]

mas_influyentes <- order(cook, decreasing = TRUE)[1:20]
modelo_sin_influyentes <- update(modelo_desercion, data = entrenamiento[-mas_influyentes, ])

cambio_coef <- tibble(term = names(coef(modelo_desercion)),
                      completo = coef(modelo_desercion),
                      sin_20 = coef(modelo_sin_influyentes),
                      ee = sqrt(diag(vcov(modelo_desercion)))) |>
  mutate(cambio_en_ee = (sin_20 - completo) / ee) |>
  arrange(desc(abs(cambio_en_ee)))

cambio_coef |> mutate(across(where(is.numeric), \(x) signif(x, 3))) |> kable()

mayor_cambio <- cambio_coef |> filter(term != "(Intercept)") |> slice(1)


# [influencia-texto]
#| output: asis

p_completo <- coef(summary(modelo_desercion))[, 4]
p_sin_20   <- coef(summary(modelo_sin_influyentes))[, 4]
cambian_signif <- setdiff(names(p_completo)[(p_completo < 0.05) != (p_sin_20 < 0.05)], "(Intercept)")
predictoras <- cambio_coef |> filter(term != "(Intercept)")
n_mas_1ee <- sum(abs(predictoras$cambio_en_ee) > 1)
cambian_signo <- predictoras$term[sign(predictoras$completo) != sign(predictoras$sin_20)]

cat(glue::glue(
  "Sacar 20 de {num(n_ent, 0)} clientes (el {num(100 * 20 / n_ent, 1)} %) mueve {n_mas_1ee} coeficientes de ",
  "predictoras más de un error estándar; el que más se mueve es `{mayor_cambio$term}`, ",
  "{num(abs(mayor_cambio$cambio_en_ee), 1)} errores estándar. ",
  "{if (length(cambian_signo) == 0) 'Ningún coeficiente cambia de signo' else paste0('Cambian de signo ', paste(cambian_signo, collapse = ', '))}",
  "{if (length(cambian_signif) == 0) ' y ningún coeficiente pasa de significativo a no significativo al 5 %, ni al revés.' else paste0('; cambia de significación al 5 % ', paste0('`', cambian_signif, '`', collapse = ', '), '.')} ",
  "{if (length(cambian_signif) == 0) 'Es decir, un puñado de clientes tiene una influencia apreciable sobre la magnitud de algunos coeficientes, pero no sobre las conclusiones cualitativas del modelo.' else paste0('El cambio de significación es de un coeficiente que ya estaba cerca del límite: su p-valor pasa de ', paste(num(p_completo[cambian_signif], 3), collapse = ', '), ' a ', paste(num(p_sin_20[cambian_signif], 3), collapse = ', '), '. Un puñado de clientes tiene una influencia apreciable sobre la magnitud de algunos coeficientes, y alcanza para mover un p-valor que estaba en el borde, pero no cambia signos ni la lectura general del modelo.')} ",
  "Eso no justifica eliminarlos: ",
  "son clientes reales, mal ajustados por el modelo, y el reajuste sirve para saber cuánto dependen ",
  "los resultados de ellos, no para descartarlos."
))


## Separación perfecta {#sec-separacion} -----------------------------------

# [separacion-ejemplo]

ejemplo_separado <- tibble(x = 1:10, y = as.numeric(x > 5))
modelo_separado <- suppressWarnings(glm(y ~ x, family = binomial, data = ejemplo_separado))

coef(summary(modelo_separado))
c(iteraciones = modelo_separado$iter, convergio = modelo_separado$converged)


# [separacion-desercion]

celdas_minimas <- entrenamiento |>
  select(where(is.factor)) |>
  pivot_longer(-RESTADO, names_to = "variable", values_to = "nivel") |>
  count(variable, nivel, RESTADO) |>
  group_by(variable) |>
  summarise(celda_minima = min(n))

celdas_minimas

c(ee_maximo = max(coefs$std.error), p_min = min(p_ent), p_max = max(p_ent),
  iteraciones = modelo_desercion$iter)

stopifnot(all(celdas_minimas$celda_minima > 0), modelo_desercion$converged)


# Selección de variables {#sec-seleccion} ----------------------------------


## Eliminación hacia atrás por AIC -----------------------------------------

# [step-aic]

modelo_step_aic <- step(modelo_desercion, direction = "backward", trace = 1)

vars_aic <- attr(terms(modelo_step_aic), "term.labels")


# [step-aic-texto]
#| output: asis

aum_aic <- tabla_drop1$AIC[-1] - AIC(modelo_desercion)
names(aum_aic) <- rownames(tabla_drop1)[-1]

if (length(vars_aic) == length(nombres_variables)) {
  cat(glue::glue(
    "`step()` no saca ninguna variable: el modelo completo ya es el de menor AIC. Es lo mismo que se ",
    "lee en la columna `AIC` de `drop1()`: sacar cualquiera de las trece predictoras aumenta el AIC. ",
    "El menor aumento es el de `{names(which.min(aum_aic))}`, que lleva el AIC de ",
    "{num(AIC(modelo_desercion), 1)} a {num(AIC(modelo_desercion) + min(aum_aic), 1)}. ",
    "El modelo completo es, por este criterio, el modelo final."
  ))
} else {
  cat(glue::glue(
    "`step()` saca {paste0('`', setdiff(nombres_variables, vars_aic), '`', collapse = ', ')}."
  ))
}


## Un criterio más exigente: BIC -------------------------------------------

# [step-bic]

modelo_step_bic <- step(modelo_desercion, direction = "backward", k = log(n_ent), trace = 0)
formula(modelo_step_bic)

vars_bic <- attr(terms(modelo_step_bic), "term.labels")
vars_fuera_bic <- setdiff(nombres_variables, vars_bic)


# [bic-vs-completo]

lrt_bic <- anova(modelo_step_bic, modelo_desercion, test = "LRT")
lrt_bic

p_bic <- lrt_bic$`Pr(>Chi)`[2]


# Evaluación del modelo ----------------------------------------------------


## Predicciones: `type = "link"` y `type = "response"` ---------------------

# [predict-tipos]

logit_test <- predict(modelo_desercion, newdata = test, type = "link")
prob_test  <- predict(modelo_desercion, newdata = test, type = "response")

stopifnot(isTRUE(all.equal(prob_test, plogis(logit_test))))

tibble(logit = logit_test, plogis_logit = plogis(logit_test), probabilidad = prob_test) |>
  head()

range(logit_test)
range(prob_test)


## Regla de clasificación y matriz de confusión ----------------------------

# [clases-predichas]

clase_train <- factor(ifelse(p_ent > 0.5, "CANCELADO", "ACTIVO"),
                      levels = levels(datos$RESTADO))
clase_test  <- factor(ifelse(prob_test > 0.5, "CANCELADO", "ACTIVO"),
                      levels = levels(datos$RESTADO))

stopifnot(!anyNA(clase_test), length(clase_test) == nrow(test))


# [tablas-confusion]

# Filas: condición real. Columnas: grupo predicho. Clase positiva primero.
tabla_confusion <- function(real, predicho) {
  orden <- c("CANCELADO", "ACTIVO")
  table(Real = factor(real, levels = orden), Predicho = factor(predicho, levels = orden))
}

tabla_train <- tabla_confusion(entrenamiento$RESTADO, clase_train)
tabla_test  <- tabla_confusion(test$RESTADO, clase_test)

esperada_train <- matrix(c(1850, 500, 187, 3234), 2, byrow = TRUE)
esperada_test  <- matrix(c(1303, 337, 138, 2153), 2, byrow = TRUE)

stopifnot(
  all(unclass(tabla_train) == esperada_train),
  all(unclass(tabla_test) == esperada_test),
  sum(tabla_test) == 3931
)

tabla_train
tabla_test


## Medidas de desempeño ----------------------------------------------------

# [medidas]

medidas_desempenio <- function(tabla) {
  VP <- tabla["CANCELADO", "CANCELADO"]; FN <- tabla["CANCELADO", "ACTIVO"]
  FP <- tabla["ACTIVO", "CANCELADO"];    VN <- tabla["ACTIVO", "ACTIVO"]
  n  <- sum(tabla)
  exactitud <- (VP + VN) / n
  sensib    <- VP / (VP + FN)
  precision <- VP / (VP + FP)
  p_e       <- sum(rowSums(tabla) * colSums(tabla)) / n^2
  c(exactitud = exactitud, tasa_error = 1 - exactitud,
    sensibilidad = sensib, especificidad = VN / (VN + FP),
    precision_VPP = precision, VPN = VN / (VN + FN),
    F1 = 2 * precision * sensib / (precision + sensib),
    kappa = (exactitud - p_e) / (1 - p_e))
}

medidas <- tibble(medida = names(medidas_desempenio(tabla_train)),
                  entrenamiento = medidas_desempenio(tabla_train),
                  test = medidas_desempenio(tabla_test))

m_test  <- medidas_desempenio(tabla_test)
m_train <- medidas_desempenio(tabla_train)

stopifnot(
  round(m_train["tasa_error"], 4) == 0.1190,
  round(m_test["tasa_error"], 4) == 0.1208,
  round(m_test["exactitud"], 4) == 0.8792
)

medidas |> kable(digits = 4)


# [caret]

cm_test <- confusionMatrix(data = clase_test, reference = test$RESTADO, positive = "CANCELADO")
cm_test

stopifnot(
  isTRUE(all.equal(unname(cm_test$overall["Accuracy"]), unname(m_test["exactitud"]))),
  isTRUE(all.equal(unname(cm_test$overall["Kappa"]), unname(m_test["kappa"]))),
  isTRUE(all.equal(unname(cm_test$byClass["Sensitivity"]), unname(m_test["sensibilidad"]))),
  isTRUE(all.equal(unname(cm_test$byClass["Specificity"]), unname(m_test["especificidad"]))),
  isTRUE(all.equal(unname(cm_test$byClass["F1"]), unname(m_test["F1"])))
)


## La tasa de no información -----------------------------------------------

# [nir]

nir <- max(prop.table(table(test$RESTADO)))

stopifnot(isTRUE(all.equal(nir, unname(cm_test$overall["AccuracyNull"]))))

c(NIR = nir, exactitud = unname(m_test["exactitud"]),
  p_valor_exactitud_mayor_NIR = unname(cm_test$overall["AccuracyPValue"]))


## Entrenamiento contra test -----------------------------------------------

# [train-vs-test]

c(error_entrenamiento = unname(m_train["tasa_error"]), error_test = unname(m_test["tasa_error"]),
  diferencia = unname(m_test["tasa_error"] - m_train["tasa_error"]))


## Curva ROC y AUC ---------------------------------------------------------

# [roc]

roc_test <- roc(test$RESTADO, prob_test, levels = c("ACTIVO", "CANCELADO"), direction = "<",
                quiet = TRUE)
auc_test <- as.numeric(auc(roc_test))

stopifnot(round(auc_test, 4) == 0.9182)

puntos_roc <- coords(roc_test, "all", transpose = FALSE)

ggplot(puntos_roc, aes(1 - specificity, sensitivity)) +
  geom_path(colour = "firebrick", linewidth = 1) +
  geom_abline(linetype = "dashed", colour = "grey50") +
  coord_equal() +
  labs(title = "Curva ROC en la muestra de test",
       subtitle = paste0("AUC = ", num(auc_test, 4)),
       x = "1 - especificidad (tasa de falsos positivos)", y = "Sensibilidad")


# [auc-pares]

es_cancelado <- test$RESTADO == "CANCELADO"
n1 <- sum(es_cancelado); n0 <- sum(!es_cancelado)
auc_pares <- (sum(rank(prob_test)[es_cancelado]) - n1 * (n1 + 1) / 2) / (n1 * n0)

stopifnot(isTRUE(all.equal(auc_pares, auc_test)))

c(AUC_pROC = auc_test, proporcion_de_pares = auc_pares)


## El punto de corte {#sec-corte} ------------------------------------------

# [youden]

youden <- coords(roc_test, "best", best.method = "youden", transpose = FALSE)
youden

stopifnot(
  abs(youden$threshold - 0.5117) < 1e-3,
  round(youden$sensitivity, 3) == 0.795,
  round(youden$specificity, 3) == 0.943
)


# [tabla-umbrales]

umbrales <- c(0.2, 0.3, 0.4, 0.5, youden$threshold, 0.6, 0.7, 0.8)

coords(roc_test, x = umbrales, input = "threshold",
       ret = c("threshold", "sensitivity", "specificity", "ppv", "accuracy"),
       transpose = FALSE) |>
  kable(digits = 3)


# [youden-entrenamiento]

roc_train <- roc(entrenamiento$RESTADO, p_ent, levels = c("ACTIVO", "CANCELADO"),
                 direction = "<", quiet = TRUE)
youden_train <- coords(roc_train, "best", best.method = "youden", transpose = FALSE)
youden_train


# [youden-entrenamiento-en-test]

clase_test_youden <- ifelse(prob_test > youden_train$threshold, "CANCELADO", "ACTIVO")
m_youden_train <- medidas_desempenio(tabla_confusion(test$RESTADO, clase_test_youden))
m_youden_train[c("sensibilidad", "especificidad", "tasa_error")]


# [youden-default]

roc_default <- roc(Default$default, prob_default, levels = c("No", "Yes"), direction = "<",
                   quiet = TRUE)
youden_default <- coords(roc_default, "best", best.method = "youden", transpose = FALSE)

tibble(umbral = c(0.5, youden_default$threshold),
       sensibilidad = c(sens_default, youden_default$sensitivity),
       especificidad = c(tabla_default["No", "No"] / sum(tabla_default["No", ]),
                         youden_default$specificity)) |>
  kable(digits = 3)


# Logística contra discriminante lineal {#sec-lda} -------------------------


## Las diferencias de fondo ------------------------------------------------


## La comparación con datos ------------------------------------------------

# [lda-qda]

modelo_lda <- lda(RESTADO ~ ., data = entrenamiento)
modelo_qda <- qda(RESTADO ~ ., data = entrenamiento)

pred_lda <- predict(modelo_lda, newdata = test)
pred_qda <- predict(modelo_qda, newdata = test)

resumen_metodo <- function(metodo, clase, prob) {
  tabla <- tabla_confusion(test$RESTADO, clase)
  m <- medidas_desempenio(tabla)
  tibble(metodo,
         error_test = m["tasa_error"], sensibilidad = m["sensibilidad"],
         especificidad = m["especificidad"],
         AUC = as.numeric(auc(roc(test$RESTADO, prob, levels = c("ACTIVO", "CANCELADO"),
                                  direction = "<", quiet = TRUE))))
}

comparacion <- bind_rows(
  resumen_metodo("Logística", clase_test, prob_test),
  resumen_metodo("Discriminante lineal", pred_lda$class, pred_lda$posterior[, "CANCELADO"]),
  resumen_metodo("Discriminante cuadrático", pred_qda$class, pred_qda$posterior[, "CANCELADO"])
) |>
  arrange(error_test)

comparacion |> kable(digits = 4)


# [lda-qda-texto]
#| output: asis

mejor_error <- comparacion$metodo[which.min(comparacion$error_test)]
mejor_auc   <- comparacion$metodo[which.max(comparacion$AUC)]
fila <- function(m) comparacion |> filter(metodo == m)

cat(glue::glue(
  "La logística tiene {if (mejor_error == 'Logística') 'el menor' else 'un mayor'} error de test ",
  "({num(fila('Logística')$error_test, 4)}) y {if (mejor_auc == 'Logística') 'el mayor' else 'un menor'} ",
  "AUC ({num(fila('Logística')$AUC, 4)}). ",
  "El discriminante lineal queda cerca (error {num(fila('Discriminante lineal')$error_test, 4)}, ",
  "AUC {num(fila('Discriminante lineal')$AUC, 4)}): ",
  "{num(100 * (fila('Discriminante lineal')$error_test - fila('Logística')$error_test), 1)} puntos más de error. ",
  "El cuadrático es el peor de los tres (error {num(fila('Discriminante cuadrático')$error_test, 4)}, ",
  "AUC {num(fila('Discriminante cuadrático')$AUC, 4)})."
))

