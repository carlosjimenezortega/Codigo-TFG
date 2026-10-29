#------------------------CARGAMOS LA IMAGEN---------------------------------

# Cargamos las librerías necesarias
install.packages(c("magick", "pdftools"))
library(imager)
library(rrcov)
library(HDclassif)
library(magick)
library(pdftools)

# Cargamos la imagen original
img_original <- load.image("i4edit.jpg")

# Mostramos la imagen original
par(mar = c(2, 2, 3, 2))
plot(img_original, main = "Imagen Original")

# Ajustamos el tamaño total al múltiplo de 8 más cercano
ancho_real <- width(img_original)
alto_real <- height(img_original)

ancho_ajustado <- floor(ancho_real / 8) * 8
alto_ajustado <- floor(alto_real / 8) * 8

# Ajustamos la imagen
img_trabajo <- imsub(img_original, x <= ancho_ajustado, y <= alto_ajustado)

#---------------------CONSTRUCCIÓN DE LA MATRIZ---------------------------

# Definimos el tamaño del bloque (b x b)
b <- 8 

# Convertimos la imagen en un array
img_array <- as.array(img_trabajo)

# Calculamos las dimensiones de la nueva matriz X para toda la imagen
num_bloques_ancho <- ancho_ajustado / b
num_bloques_alto <- alto_ajustado / b
n <- num_bloques_ancho * num_bloques_alto # Total de bloques 
p <- b * b * 3 # Dimensión de cada vector (192)

cat("La matriz X tendrá", n, "observaciones y", p, "variables.\n")

# Inicializamos la matriz X con ceros
X <- matrix(0, nrow = n, ncol = p)

# Extraemos cada bloque de toda la imagen, lo aplanamos y rellenamos X
contador <- 1
for (i in 1:num_bloques_ancho) {
  for (j in 1:num_bloques_alto) {
    x_inicio <- (i - 1) * b + 1
    x_fin <- i * b
    y_inicio <- (j - 1) * b + 1
    y_fin <- j * b
    
    bloque <- img_array[x_inicio:x_fin, y_inicio:y_fin, 1, ]
    X[contador, ] <- as.numeric(bloque)
    contador <- contador + 1
  }
}

# Análisis dimensiones y canales
cat("Dimensiones originales:", ancho_real, "x", alto_real, "\n")
cat("Dimensiones ajustadas:", ancho_ajustado, "x", alto_ajustado, "\n")
cat("Bloques ancho x alto:", num_bloques_ancho, "x", num_bloques_alto, "\n")

media_R <- mean(X[, 1:64])
media_G <- mean(X[, 65:128])
media_B <- mean(X[, 129:192])
sd_R <- sd(as.vector(X[, 1:64]))
sd_G <- sd(as.vector(X[, 65:128]))
sd_B <- sd(as.vector(X[, 129:192]))

cat("Canal Rojo - Media:", round(media_R, 4), "SD:", round(sd_R, 4), "\n")
cat("Canal Verde - Media:", round(media_G, 4), "SD:", round(sd_G, 4), "\n")
cat("Canal Azul - Media:", round(media_B, 4), "SD:", round(sd_B, 4), "\n")

#---------------------------LIMPIEZA DE RUIDO----------------------------------

# Extraemos una muestra aleatoria
set.seed(42) 
idx_muestra <- sample(1:nrow(X), size = floor(0.15 * nrow(X)))
X_muestra <- X[idx_muestra, ]

# Aplicamos técnicas ACP Robusto 
modelo_robusto <- PcaProj(X_muestra, k = 3)

# Centro robusto y pesos
centro_robusto <- modelo_robusto@center
pesos <- modelo_robusto@loadings

# Proyectamos la matriz X
X_centrada <- scale(X, center = centro_robusto, scale = FALSE)
puntuaciones <- X_centrada %*% pesos 

# Reconstruimos la matriz y residuos
X_reconstruida <- puntuaciones %*% t(pesos)
residuos <- X_centrada - X_reconstruida
distancias_ortogonales <- sqrt(rowSums(residuos^2))

# Umbral teórico Wilson-Hilferty
dist_ort_norm <- distancias_ortogonales^(2/3)
media_estim <- median(dist_ort_norm)
desv_estim <- mad(dist_ort_norm)
cuantil_z <- qnorm(0.975)

umbral_teorico <- (media_estim + desv_estim * cuantil_z)^(3/2)
outliers <- distancias_ortogonales > umbral_teorico
cat("Se han detectado", sum(outliers), "bloques atípicos.\n")

# Análisis limpieza
cat("Tamaño de la submuestra:", nrow(X_muestra), "\n")
cat("Autovalores robustos:", round(modelo_robusto@eigenvalues, 4), "\n")
cat("Mediana (d^2/3):", round(media_estim, 4), "MAD (d^2/3):", round(desv_estim, 4), "\n")
cat("Umbral Wilson-Hilferty:", round(umbral_teorico, 4), "\n")
cat("Porcentaje de outliers:", round(mean(outliers) * 100, 2), "%\n")

cat("Media R/G/B en bloques normales:",
    round(mean(X[!outliers, 1:64]), 4),
    round(mean(X[!outliers, 65:128]), 4),
    round(mean(X[!outliers, 129:192]), 4), "\n")

cat("Media R/G/B en bloques outliers:",
    round(mean(X[outliers, 1:64]), 4),
    round(mean(X[outliers, 65:128]), 4),
    round(mean(X[outliers, 129:192]), 4), "\n")

cat("Media R/G/B del centro robusto:",
    round(mean(centro_robusto[1:64]), 4),
    round(mean(centro_robusto[65:128]), 4),
    round(mean(centro_robusto[129:192]), 4), "\n")

# Gráfico outliers

pdf("diagnostico_outliers.pdf", width = 7, height = 5)
par(mar = c(4, 4, 3, 2))
colores_pts <- ifelse(outliers, "red", "darkolivegreen4")
pch_pts <- ifelse(outliers, 4, 16)
plot(1:n, distancias_ortogonales,
     col = colores_pts, pch = pch_pts, cex = 0.5,
     xlab = "Índice del bloque (i)",
     ylab = "Distancia ortogonal",
     main = "Diagnóstico de Outliers")
abline(h = umbral_teorico, col = "darkred", lwd = 2, lty = 2)
legend("topright",
       legend = c("Bloques normales", "Bloques atípicos", "Umbral crítico"),
       col = c("darkolivegreen4", "red", "darkred"),
       pch = c(16, 4, NA), lty = c(NA, NA, 2), lwd = c(NA, NA, 2),
       cex = 0.85, bg = "white")
dev.off()

# Biplot robusto

pdf("biplot_robusto.pdf", width = 7, height = 5)
par(mar = c(4, 4, 3, 2))
plot(puntuaciones[, 1], puntuaciones[, 2],
     col = colores_pts, pch = pch_pts, cex = 0.45,
     xlab = "CP1 Robusta", ylab = "CP2 Robusta",
     main = "Biplot Robusto")
abline(h = 0, v = 0, lty = 3, col = "gray60")

flecha_R <- colMeans(pesos[1:64, 1:2])
flecha_G <- colMeans(pesos[65:128, 1:2])
flecha_B <- colMeans(pesos[129:192, 1:2])
escala_flechas <- max(abs(puntuaciones[, 1:2])) * 9.5

arrows(0, 0, flecha_R[1] * escala_flechas, flecha_R[2] * escala_flechas,
       col = "darkred", lwd = 2, length = 0.1)
arrows(0, 0, flecha_G[1] * escala_flechas, flecha_G[2] * escala_flechas,
       col = "darkgreen", lwd = 2, length = 0.1)
arrows(0, 0, flecha_B[1] * escala_flechas, flecha_B[2] * escala_flechas,
       col = "darkblue", lwd = 2, length = 0.1)

text(flecha_R[1] * escala_flechas * 1.12, flecha_R[2] * escala_flechas * 1.12,
     "Canal R", col = "darkred", cex = 0.85, font = 2)
text(flecha_G[1] * escala_flechas * 1.12, flecha_G[2] * escala_flechas * 1.12,
     "Canal G", col = "darkgreen", cex = 0.85, font = 2)
text(flecha_B[1] * escala_flechas * 1.12, flecha_B[2] * escala_flechas * 1.12,
     "Canal B", col = "darkblue", cex = 0.85, font = 2)

legend("topright",
       legend = c("Bloques normales", "Bloques atípicos"),
       col = c("darkolivegreen4", "red"),
       pch = c(16, 4), cex = 0.85, bg = "white")
dev.off()

# Limpiamos tales bloques
X_limpio <- X
centro_robusto <- modelo_robusto@center

for (i in 1:nrow(X)) {
  if (outliers[i]) {
    X_limpio[i, ] <- centro_robusto
  }
}

# Reconstruimos la imagen limpia
img_limpia_array <- as.array(img_trabajo) 
contador <- 1

for (i in 1:num_bloques_ancho) {
  for (j in 1:num_bloques_alto) {
    x_inicio <- (i - 1) * b + 1
    x_fin <- i * b
    y_inicio <- (j - 1) * b + 1
    y_fin <- j * b
    
    bloque_vector <- X_limpio[contador, ]
    bloque_array <- array(bloque_vector, dim = c(b, b, 1, 3))
    img_limpia_array[x_inicio:x_fin, y_inicio:y_fin, 1, ] <- bloque_array
    
    contador <- contador + 1
  }
}

img_limpia <- as.cimg(img_limpia_array)
plot(img_limpia, main = "Imagen Limpia")

#-------------------------------SEGMENTACIÓN------------------------------------

# Ejecutamos el modelo de mezclas con K=3
modelo_mppca <- hddc(X_limpio, K = 3) 

# Extraemos la clasificación final
clases_bloques <- modelo_mppca$class

# Análisis clases
tabla_clases <- table(clases_bloques)
prop_clases <- prop.table(tabla_clases) * 100

cat("Bloques por clase (1, 2, 3):", as.numeric(tabla_clases), "\n")
cat("Porcentaje por clase (1, 2, 3):", round(as.numeric(prop_clases), 2), "\n")
cat("Total deforestado (Clases 1 + 2):", round(sum(prop_clases[1:2]), 2), "%\n")
cat("Total selva densa (Clase 3):", round(prop_clases[3], 2), "%\n")

# Reconstruimos con la segmentación 
# Tenemos 3 clases, pintamos una de verde pino y dos de verde césped
colores_clases <- list(
  c(0.55, 0.71, 0.45),  # Verde césped
  c(0.55, 0.71, 0.45),  # Verde césped
  c(0.12, 0.28, 0.12)   # Verde pino
)

img_segmentada_array <- array(0, dim = c(ancho_ajustado, alto_ajustado, 1, 3))
contador <- 1

for (i in 1:num_bloques_ancho) {
  for (j in 1:num_bloques_alto) {
    x_inicio <- (i - 1) * b + 1
    x_fin <- i * b
    y_inicio <- (j - 1) * b + 1
    y_fin <- j * b
    
    clase_asignada <- clases_bloques[contador]
    color_bloque <- colores_clases[[clase_asignada]]
    
    img_segmentada_array[x_inicio:x_fin, y_inicio:y_fin, 1, 1] <- color_bloque[1]
    img_segmentada_array[x_inicio:x_fin, y_inicio:y_fin, 1, 2] <- color_bloque[2]
    img_segmentada_array[x_inicio:x_fin, y_inicio:y_fin, 1, 3] <- color_bloque[3]
    
    contador <- contador + 1
  }
}

img_segmentada <- as.cimg(img_segmentada_array)

#-----------------------------VISUALIZACIÓN FINAL-------------------------------
par(mfrow = c(1, 3), mar = c(1, 1, 2, 1))

plot(img_trabajo, main = "1. Imagen Original con Marca de Agua")
plot(img_limpia, main = "2. Imagen Limpia")
plot(img_segmentada, main = "3. Imagen Segmentada Final")

par(mfrow = c(1, 1))

# Guardamos la imagen limpia
save.image(img_limpia, "imagen_limpia.jpg")

# Guardamos la imagen segmentada
save.image(img_segmentada, "imagen_segmentada.jpg")