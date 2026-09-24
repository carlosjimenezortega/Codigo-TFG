#-----------------------OBTENCIÓN COMPONENTES PRINCIPALES-----------------------

# Cargamos la librería necesaria
library(xtable)

# Cargamos la matriz de datos en bruto X tilde
X_tilde <- read.csv("maraton.csv", header=TRUE)

# Centramos la matriz (restando la media de la columna a cada observación)
X <- scale(X_tilde, center = TRUE, scale = FALSE)

# Aplicamos la Descomposición en Valores Singulares (SVD)
descomposicion <- svd(X)

# Extraemos los tres argumentos
U <- descomposicion$u
d <- descomposicion$d
V <- descomposicion$v

# Construimos la matriz diagonal D a partir del vector d
D_mat <- diag(d)

# Calculamos la matriz de puntuaciones Z multiplicando U por la matriz diagonal D
Z <- U %*% D_mat

# Calculamos los valores propios (varianzas) a partir de los valores singulares
n <- nrow(X)
valores_propios <- (d^2) / (n - 1)

# Calculamos el porcentaje de varianza explicada por cada componente
varianza_explicada <- (valores_propios / sum(valores_propios)) * 100

# Calculamos el porcentaje acumulado (t_m)
varianza_acumulada <- cumsum(varianza_explicada)

# Mostramos la tabla con los resultados
tabla_varianzas <- data.frame(
  Componente = 1:10,
  Valor_Propio = round(valores_propios, 4),
  Varianza_Porcentaje = round(varianza_explicada, 2),
  Porcentaje_Acumulada = round(varianza_acumulada, 2)
)
print(tabla_varianzas)

# Generamos el Scree Graph 
pdf("scree_graph_maraton.pdf", width = 7, height = 5)
plot(1:10, valores_propios, type="b", pch=19, col="blue",
     main="Scree Graph",
     xlab="Componente Principal (k)",
     ylab="Varianza (Valor Propio l_k)")
grid() # Cuadrícula
dev.off()

# Calculamos e imprimimos la media de los valores propios
media_valores_propios <- mean(valores_propios)
cat("\nLa media de los valores propios es:", media_valores_propios, "\n")

0.7*media_valores_propios

valores_propios

#--------------INTERPRETACIÓN DE LAS COMPONENTES PRINCIPALES--------------------

# Extraemos los coeficientes para las dos primeras componentes
coeficientes <- V[, 1:2]

# Asignamos los nombres de los 10 puntos de control
puntos_control <- c("5K", "10K", "15K", "20K", "Media", 
                    "25K", "30K", "35K", "40K", "Meta")
rownames(coeficientes) <- puntos_control
colnames(coeficientes) <- c("CP1", "CP2")

# Mostramos la tabla exacta 
cat("\n--- TABLA DE COEFICIENTES EXACTOS ---\n")
print(round(coeficientes, 2))

# Creamos la tabla simplificada 
tabla_simplificada <- ifelse(coeficientes > 0.15, "+", 
                             ifelse(coeficientes < -0.15, "-", ""))

cat("\n--- TABLA SIMPLIFICADA DE SIGNOS ---\n")
print(tabla_simplificada, quote = FALSE)

#------------------------------------BIPLOT-------------------------------------


# Definimos alfa 
alfa <- 0 

# Construimos D^alfa y D^(1-alfa)
D_alfa <- diag(d^alfa)
D_1_menos_alfa <- diag(d^(1 - alfa))

# Construimos las matrices G y H
G <- U %*% D_alfa
H <- V %*% D_1_menos_alfa

# Truncamos
G_estrella <- G[, 1:2]
H_estrella <- H[, 1:2]

# Factor de reescalado
factor_escala <- sqrt(n - 1)
G_estrella_reescalado <- G_estrella * factor_escala
H_estrella_reescalado <- H_estrella / factor_escala

# Dibujar el Biplot:
# Límites
xlim <- range(c(G_estrella_reescalado[,1], H_estrella_reescalado[,1]))
ylim <- range(c(G_estrella_reescalado[,2], H_estrella_reescalado[,2]))

# Trazamos las observaciones
pdf("biplot.pdf", width = 7, height = 5)
plot(G_estrella_reescalado, pch = 19, col = "lightblue", 
     xlab = "CP1", ylab = "CP2",
     xlim = xlim, ylim = ylim, main = expression(paste("Biplot  (", alpha, " = 0)")))
grid()
abline(h=0, v=0, lty=2, col="gray")

# Trazamos las variables 
arrows(x0 = 0, y0 = 0, x1 = H_estrella_reescalado[,1], y1 = H_estrella_reescalado[,2], 
       col = "darkred", length = 0.1, lwd = 1.5)

# Añadimos etiquetas ajustando 'pos' y 'cex'
text(H_estrella_reescalado[,1], H_estrella_reescalado[,2], labels = puntos_control, 
     col = "black", pos = 2, offset = 0.3, cex = 0.7, font = 2)
dev.off()

xtable(coeficientes)
xtable(tabla_simplificada)
xtable(tabla_varianzas)

#-------------------------- BIPLOT DEGRADADO (CP1) -----------------------------

# Extraemos las puntuaciones de la primera componente principal
cp1_puntuaciones <- G_estrella_reescalado[, 1]

# Generador del degradado
degradado <- colorRampPalette(c("#440154", "#21908C", "#FDE725"))
paleta <- degradado(100)

# Normalizamos los valores de CP2
cp1_normal <- (cp1_puntuaciones - min(cp1_puntuaciones)) / (max(cp1_puntuaciones) - min(cp1_puntuaciones))
indices_color <- round(cp1_normal * 99) + 1

# Asignamos el color a cada corredor
colores <- paleta[indices_color]

# Generamos el Biplot con degradado
pdf("biplot_degradado_cp1.pdf", width = 8, height = 6)

# Dibujamos el Biplot
plot(G_estrella_reescalado, pch = 19, col = colores, 
     xlab = "CP1", ylab = "CP2",
     xlim = xlim, ylim = ylim, 
     main = expression(paste("Biplot con gradiente CP1 (", alpha, " = 0)")),
     panel.first = list(grid(), abline(h = 0, v = 0, lty = 2, col = "gray")))

# Trazamos las variables 
arrows(x0 = 0, y0 = 0, x1 = H_estrella_reescalado[, 1], y1 = H_estrella_reescalado[, 2], 
       col = "darkred", length = 0.1, lwd = 1.5)

# Añadimos etiquetas 
text(H_estrella_reescalado[, 1], H_estrella_reescalado[, 2], labels = puntos_control, 
     col = "black", pos = 2, offset = 0.3, cex = 0.7, font = 2)

# Leyenda
legend("bottomleft", 
       legend = c("Mejor nivel", "Nivel medio", "Peor nivel"),
       col = c("#FDE725", paleta[50], "#440154"), 
       pch = 19, cex = 0.75, bty = "n")

dev.off()

#---------------------- BIPLOT CON COLORES (CP2) -------------------------------

# Extraemos las puntuaciones de la segunda componente principal
cp2_puntuaciones <- G_estrella_reescalado[, 2]

# Umbrales para separar las tres vertientes
grupo <- ifelse(cp2_puntuaciones > 0.75, "De menos a más",
                           ifelse(cp2_puntuaciones < -0.75, "El muro", "Constante"))

# Asignamos los colores:
colores <- ifelse(grupo == "De menos a más", "dodgerblue3",
                             ifelse(grupo == "El muro", "firebrick3", "goldenrod2"))

# Generamos el Biplot con colores
pdf("biplot_colores_cp2.pdf", width = 8, height = 6)
plot(G_estrella_reescalado, pch = 19, col = colores, 
     xlab = "CP1", ylab = "CP2",
     xlim = xlim, ylim = ylim, 
     main = expression(paste("Biplot estratificado CP2 (", alpha, " = 0)")),
     panel.first = list(grid(), abline(h = 0, v = 0, lty = 2, col = "gray")))

# Trazamos las variables
arrows(x0 = 0, y0 = 0, x1 = H_estrella_reescalado[, 1], y1 = H_estrella_reescalado[, 2], 
       col = "darkred", length = 0.1, lwd = 1.5)

# Añadimos etiquetas
text(H_estrella_reescalado[, 1], H_estrella_reescalado[, 2], labels = puntos_control, 
     col = "black", pos = 2, offset = 0.3, cex = 0.7, font = 2)

# Leyenda
legend("bottomleft", 
       legend = c("De menos a más (CP2 > 0.75)", "Constante (-0.75 <= CP2 <= 0.75)", "El muro (CP2 < -0.75)"),
       col = c("dodgerblue3", "goldenrod2", "firebrick3"), 
       pch = 19, cex = 0.75, bty = "n")
dev.off()

# #-------------------------- BIPLOT DEGRADADO (CP2) -----------------------------
# 
# # Extraemos las puntuaciones de la segunda componente principal
# cp2_puntuaciones <- G_estrella_reescalado[, 2]
# 
# # Generador del degradado
# degradado <- colorRampPalette(c("firebrick3", "gray90", "dodgerblue3"))
# paleta <- degradado(100)
# 
# # Normalizamos los valores de CP2
# cp2_normal <- (cp2_puntuaciones - min(cp2_puntuaciones)) / (max(cp2_puntuaciones) - min(cp2_puntuaciones))
# indices_color <- round(cp2_normal * 99) + 1
# 
# # Asignamos el color a cada corredor
# colores <- paleta[indices_color]
# 
# # Generamos el Biplot con degradado
# pdf("biplot_degradado.pdf", width = 8, height = 6)
# 
# # Dibujamos el Biplot
# plot(G_estrella_reescalado, pch = 19, col = colores,
#      xlab = "CP1", ylab = "CP2",
#      xlim = xlim, ylim = ylim,
#      main = expression(paste("Biplot con gradiente CP2 (", alpha, " = 0)")),
#      panel.first = list(grid(), abline(h = 0, v = 0, lty = 2, col = "gray")))
# 
# # Trazamos las variables
# arrows(x0 = 0, y0 = 0, x1 = H_estrella_reescalado[, 1], y1 = H_estrella_reescalado[, 2],
#        col = "darkred", length = 0.1, lwd = 1.5)
# 
# # Añadimos etiquetas
# text(H_estrella_reescalado[, 1], H_estrella_reescalado[, 2], labels = puntos_control,
#      col = "black", pos = 2, offset = 0.3, cex = 0.7, font = 2)
# 
# # Leyenda
# legend("bottomleft",
#        legend = c("De menos a más", "Constante", "El muro"),
#        col = c("dodgerblue3", "gray90", "firebrick3"),
#        pch = 19, cex = 0.75, bty = "n")
# 
# dev.off()

# -------------------------- BIPLOT INTRODUCTORIO ------------------------------

# Usamos la base de datos de R mtcars
# Tomamos 15 coches y 4 variables (mpg, disp, hp, wt)
X_tilde <- mtcars[1:15, c("mpg", "disp", "hp", "wt")]

# Guardamos los nombres
variables <- colnames(X_tilde)

# Número de observaciones
n <- nrow(X_tilde)

# Centramos y escalamos la matriz de datos
X <- scale(X_tilde, center = TRUE, scale = TRUE)

# Realizamos la Descomposición en Valores Singulares (SVD)
descomposicion <- svd(X)
U <- descomposicion$u
d <- descomposicion$d
V <- descomposicion$v

# Definimos alfa 
alfa <- 0 

# Construimos D^alfa y D^(1-alfa)
D_alfa <- diag(d^alfa)
D_1_menos_alfa <- diag(d^(1 - alfa))

# Construimos las matrices G y H
G <- U %*% D_alfa
H <- V %*% D_1_menos_alfa

# Truncamos
G_estrella <- G[, 1:2]
H_estrella <- H[, 1:2]

# Factor de reescalado
factor_escala <- sqrt(n - 1)
G_estrella_reescalado <- G_estrella * factor_escala
H_estrella_reescalado <- H_estrella / factor_escala

# Dibujar el Biplot:
# Límites
xlim <- range(c(G_estrella_reescalado[,1], H_estrella_reescalado[,1]))
ylim <- range(c(G_estrella_reescalado[,2], H_estrella_reescalado[,2]))

# Trazamos las observaciones
pdf("biplot_teorico.pdf", width = 7, height = 5)
plot(G_estrella_reescalado, pch = 19, col = "lightblue", 
     xlab = "CP1", ylab = "CP2",
     xlim = xlim, ylim = ylim)
grid()
abline(h=0, v=0, lty=2, col="gray")

# Trazamos las variables 
arrows(x0 = 0, y0 = 0, x1 = H_estrella_reescalado[,1], y1 = H_estrella_reescalado[,2], 
       col = "darkred", length = 0.1, lwd = 1.5)

# Añadimos etiquetas ajustando 'pos' y 'cex'
text(H_estrella_reescalado[,1] * 1.15, H_estrella_reescalado[,2] * 1.15, 
     labels = variables, col = "black", cex = 0.7, font = 2)
dev.off()
