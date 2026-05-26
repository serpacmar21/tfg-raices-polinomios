
## Regla de Descartes

# Cuenta cuantas veces cambia de signo el polinomio
cambios_signo <- function(coefs) {
  # Quitamos los coeficientes que tengan ceros 
  coefs_sin_ceros <- coefs[coefs != 0]
  
  # Sacamos los signos y sumamos cada vez que son distintos
  signos <- sign(coefs_sin_ceros)
  cambios <- sum(diff(signos) != 0)
  return(cambios)
}

# Coeficientes del polinomio
coeficientes <- c(1, -1, -3, 0, 2, 4)

# Raices positivas en P(x)
cambios_pos <- cambios_signo(coeficientes)

# Raices negativas en P(-x)
grados <- (length(coeficientes) - 1):0
coeficientes_neg <- coeficientes * (-1)^grados
cambios_neg <- cambios_signo(coeficientes_neg)


## Metodo de McLaurin

# Coeficientes del polinomio 
coefs <- c(1, -1, -3, 0, 2, 4)
n <- length(coefs) - 1

# Cota superior (calculamos lambda)
a_n <- abs(coefs[1])
# Buscamos el maximo 
lambda <- max(abs(coefs[2:(n + 1)]) / a_n)
cota_sup <- 1 + lambda

# Cota inferior (calculamos mu)
a_0 <- abs(coefs[n + 1])
# Buscamos el maximo 
mu <- max(abs(coefs[1:n]) / a_0)
cota_inf <- 1 / (1 + mu)


## Grafica P(x) con limites de McLaurin

# Definimos el polinomio 
p <- function(x) {
  x^5 - x^4 - 3*x^3 + 2*x + 4
}

# Dibujamos la funcion 
curve(p, from = -3, to = 3, col = "red", lwd = 2, ylim = c(-5, 5),
      main = "Grafica de P(x)",
      xlab = "x", ylab = "P(x)")

# Pintamos los ejes coordenados
abline(h = 0, v = 0, col = "black")


## Algoritmo de Horner y Deflacion

# Evaluamos P(x) y sacamos el cociente por Ruffini
horner_deflacion <- function(coefs, x0) {
  n <- length(coefs)
  cociente <- numeric(n - 1) 
  cociente[1] <- coefs[1]
  
  # Aplicamos Ruffini
  for (i in 2:(n - 1)) {
    cociente[i] <- coefs[i] + x0 * cociente[i - 1]
  }
  
  # El resto es P(x0)
  resto <- coefs[n] + x0 * cociente[n - 1]
  
  return(list(valor = resto, cociente = cociente))
}

# Coeficientes de P(x)
coeficientes <- c(1, -1, -3, 0, 2, 4)

# Posibles raices (divisores del termino independiente)
candidatos <- c(1, -1, 2, -2, 4, -4)

# Probamos los candidatos
for (cand in candidatos) {
  resultado <- horner_deflacion(coeficientes, cand)
  
  cat("Probando candidato x =", cand, ", Resto =", resultado$valor, "\n")
  
  # Si da 0, es raiz y hacemos la deflacion
  if (resultado$valor == 0) {
    raiz_encontrada <- cand
    coeficientes_Q <- resultado$cociente
    
    cat("Raiz encontrada: x =", raiz_encontrada, "\n")
    cat("Coeficientes del nuevo polinomio Q(x):", coeficientes_Q, "\n")
    
    break 
  }
}


## Metodo de Sturm

# Cargamos la libreria 
# install.packages("polynom")
library(polynom)

# Polinomio Q(x) (de menor a mayor grado)
coefs_Q <- c(-2, -2, -1, 1, 1)
P0 <- polynomial(coefs_Q)
P1 <- deriv(P0)

# Generamos la secuencia de Sturm
generar_sturm <- function(p0, p1) {
  secuencia <- list(p0, p1)
  i <- 2
  
  # Repetimos hasta llegar a una constante 
  while (length(secuencia[[i]]) > 1) {
    # Calculamos el resto de la division
    division <- secuencia[[i-1]] %% secuencia[[i]]
    # Cambiamos de signo el resto
    nuevo_P <- -division
    secuencia[[i+1]] <- nuevo_P
    i <- i + 1
  }
  return(secuencia)
}

secuencia_sturm <- generar_sturm(P0, P1)

# Calculamos N(x) evaluando la secuencia en un punto
contar_Nx <- function(secuencia, valor_x) {
  evaluaciones <- sapply(secuencia, predict, valor_x)
  # Quitamos los ceros para que no sumen cambios falsos
  signos <- sign(evaluaciones[evaluaciones != 0])
  cambios <- sum(diff(signos) != 0)
  return(cambios)
}

cat("N(-2) - N(-1) = ", contar_Nx(secuencia_sturm, -2) - contar_Nx(secuencia_sturm, -1), "\n")
cat("N(1) - N(2) = ", contar_Nx(secuencia_sturm, 1) - contar_Nx(secuencia_sturm, 2), "\n")

## Metodo de Newton

# Definimos las funciones para Newton
f_Q <- function(x) { x^4 + x^3 - x^2 - 2*x - 2 }
df_Q <- function(x) { 4*x^3 + 3*x^2 - 2*x - 2 }

metodo_newton <- function(f, df, x0, tol = 1e-5, max_iter = 100) {
  x_n <- x0
  iter <- 0
  while (iter < max_iter) {
    x_n1 <- x_n - f(x_n) / df(x_n)
    iter <- iter + 1
    if (abs(x_n1 - x_n) < tol) {
      cat("Convergencia en iteracion", iter, ": x =", x_n1, "\n")
      return(x_n1)
    }
    x_n <- x_n1
  }
  return(NA)
}

raiz_negativa <- metodo_newton(f_Q, df_Q, x0 = -2)
raiz_positiva <- metodo_newton(f_Q, df_Q, x0 = 2)


## Metodo de Bairstow

# Implementacion del metodo de Bairstow
metodo_bairstow <- function(a, u0, v0, tol = 1e-5, max_iter = 100) {
  n <- length(a) - 1
  u <- u0
  v <- v0
  
  for (iter in 1:max_iter) {
    b <- numeric(n + 1); c <- numeric(n + 1)
    
    # Calculo de la secuencia b_k (coeficientes del cociente y resto)
    b[n+1] <- a[n+1]; b[n] <- a[n] + u * b[n+1]
    for (k in (n - 2):0) b[k+1] <- a[k+1] + u * b[k+2] + v * b[k+3]
    
    # Calculo de la secuencia c_k (derivadas parciales)
    c[n+1] <- 0; c[n] <- b[n+1]
    for (k in (n - 2):0) c[k+1] <- b[k+2] + u * c[k+2] + v * c[k+3]
    
    # Sistema lineal para actualizar u y v (Jacobiano)
    detJ <- c[1] * c[3] - c[2]^2
    if (detJ == 0) return(NA) 
    
    du <- (b[2] * c[2] - b[1] * c[3]) / detJ
    dv <- (b[1] * c[2] - b[2] * c[1]) / detJ
    
    u <- u + du
    v <- v + dv
    
    # Criterio de parada
    if (max(abs(du), abs(dv)) < tol) {
      # Mostramos la convergencia de u y v
      cat("Convergencia en iteracion", iter, ": u =", u, ", v =", v, "\n")
      
      # Calculo y extraccion de las raices complejas
      disc <- -u^2 - 4 * v
      if (disc > 0) {
        parte_real <- u / 2
        parte_imag <- sqrt(disc) / 2
        raiz1 <- complex(real = parte_real, imaginary = parte_imag)
        raiz2 <- complex(real = parte_real, imaginary = -parte_imag)
        
        cat("Raices complejas encontradas:\n")
        print(raiz1)
        print(raiz2)
      }
      return(list(u = u, v = v, n = iter))
    }
  }
  return(NA)
}

# Coeficientes de Q(x) 
coefs_Q <- c(-2, -2, -1, 1, 1)

# Ejecucion con punto inicial (-0.5, -0.5)
resultado <- metodo_bairstow(coefs_Q, u0 = -0.5, v0 = -0.5)










