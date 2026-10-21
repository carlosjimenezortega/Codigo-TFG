# --------------------- ILUSTRACIÓN ELIPSOIDE MEJOR AJUSTE ---------------------

library(rrcov)

# Generamos los datos simulados
set.seed(42) 

# Subpoblación densa 
x1_clean <- rnorm(80, mean = 0, sd = 1)
x2_clean <- 0.8 * x1_clean + rnorm(80, mean = 0, sd = 0.5)
datos_limpios <- cbind(x1_clean, x2_clean)

# Valores atípicos
x1_out <- rnorm(6, mean = 6, sd = 0.5)
x2_out <- rnorm(6, mean = -4, sd = 0.5)
outliers <- cbind(x1_out, x2_out)

datos <- rbind(datos_limpios, outliers)

# Estimadores
cov_clasica <- cov(datos)
centro_clasico <- colMeans(datos)

modelo_mcd <- CovMcd(datos)
cov_mcd <- getCov(modelo_mcd)
centro_mcd <- getCenter(modelo_mcd)

# Función matemática auxiliar
dibujar_elipse <- function(centro, covar, nivel = 0.975, color, tipo_linea, grosor) {
  radio <- sqrt(qchisq(nivel, df = 2))
  theta <- seq(0, 2 * pi, length.out = 100)
  caja <- cbind(cos(theta), sin(theta))
  
  eig <- eigen(covar)
  T_mat <- eig$vectors %*% diag(sqrt(eig$values))
  
  puntos <- caja %*% t(T_mat) * radio
  puntos[,1] <- puntos[,1] + centro[1]
  puntos[,2] <- puntos[,2] + centro[2]
  
  lines(puntos, col = color, lwd = grosor, lty = tipo_linea)
}

# Guardamos la gráfica y ajustes varios
pdf("grafica_elipsoide.pdf", width = 8, height = 5)

par(mar = c(1.5, 1.5, 0.5, 0.5))

plot(datos, type = "n", xlab = "", ylab = "", main = "",
     xlim = c(-4, 9), ylim = c(-6, 5), xaxt = "n", yaxt = "n")

axis(1, labels = FALSE, tick = TRUE)
axis(2, labels = FALSE, tick = TRUE)

# Dibujamos los puntos
points(datos_limpios, pch = 16, col = "darkgrey")
points(outliers, pch = 4, col = "red", lwd = 2)

# Dibujamos los elipsoides
dibujar_elipse(centro_clasico, cov_clasica, nivel = 0.975,
               color = "red", tipo_linea = 2, grosor = 2)

dibujar_elipse(centro_mcd, cov_mcd, nivel = 0.975,
               color = "blue", tipo_linea = 1, grosor = 2)

dev.off()
