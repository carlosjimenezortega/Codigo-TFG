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
