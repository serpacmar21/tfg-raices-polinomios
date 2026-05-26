# Métodos para aproximar raíces de polinomios - Implementación en R

Este repositorio contiene el código fuente y los algoritmos desarrollados para el **Trabajo Fin de Grado (TFG) en Estadística** de la Universidad de Sevilla (Curso 2025-2026).

**Autor:** Sergio Pacheco Márquez  
**Directora:** María Anguiano Moreno  
**Departamento:** Análisis Matemático (Facultad de Matemáticas, Universidad de Sevilla)

---

## 📋 Descripción del Proyecto

El objetivo principal de este trabajo es el estudio, desarrollo e implementación computacional de diversos métodos numéricos y algebraicos para la separación, acotación y aproximación de raíces (tanto reales como complejas) en ecuaciones polinómicas. 

Toda la base teórica se consolida en este repositorio mediante scripts optimizados en el lenguaje de programación **R**, permitiendo el análisis práctico completo de polinomios de forma rigurosa y eficiente.

---

## 🛠️ Métodos Implementados

El script principal incluye las funciones programadas desde cero para los siguientes métodos estudiados en la memoria:

1. **Regla de los signos de Descartes:** Estimación del número máximo de raíces reales positivas y negativas analizando las variaciones de signo de los coeficientes.
2. **Método de McLaurin:** Cálculo sistemático de cotas superiores e inferiores para el módulo de todas las raíces (reales y complejas).
3. **Algoritmo de Horner y Deflación:** Evaluación eficiente de polinomios y reducción de grado (regla de Ruffini) tras localizar una raíz exacta.
4. **Método de Sturm:** Construcción de la secuencia de Sturm mediante divisiones euclidianas sucesivas y separación exacta de raíces reales irracionales en intervalos disjuntos.
5. **Método de Newton-Raphson:** Aproximación numérica de alta precisión para raíces reales a partir de los intervalos obtenidos por Sturm.
6. **Método de Bairstow:** Resolución de sistemas no lineales mediante un esquema iterativo de Newton 2D para extraer factores cuadráticos de la forma $x^2 - ux - v$, aislando con éxito raíces complejas conjugadas.

---

## 📂 Estructura del Repositorio

* `codigo_TFG_SergioPacheco.R`: Script principal en R que contiene todas las funciones de los métodos descritos, así como la resolución y análisis práctico del polinomio de estudio de grado 5:  
  $$P(x) = x^5 - x^4 - 3x^3 + 2x + 4$$
* `README.md`: Este archivo con las instrucciones de uso y guía de ejecución.

---

## 🚀 Requisitos e Instalación

Para ejecutar este código de forma local, es necesario disponer de:
1. **R** (versión 4.0 o superior recomendada).
2. **RStudio** como entorno de desarrollo integrado (IDE).

### Dependencias
El script utiliza el paquete oficial `polynom` de CRAN para la gestión avanzada y manipulación de estructuras polinómicas en R (fundamental para el cálculo algebraico de la secuencia de Sturm). Puedes instalarlo ejecutando en la consola de RStudio:

```R
install.packages("polynom")
