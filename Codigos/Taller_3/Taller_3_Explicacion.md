# Taller 3: Productos Tensoriales de Espacios Vectoriales de Polinomios

## Introducción

Este documento explica la solución al problema de construcción de espacios tensoriales a partir de espacios vectoriales de polinomios de grado ≤ 2. El problema involucra el uso de polinomios de Legendre y el concepto de producto tensorial.

## Conceptos Teóricos

### 1. Espacios Vectoriales de Polinomios

- **P₂(x)**: Espacio vectorial de polinomios en x de grado ≤ 2
- **G₂(y)**: Espacio vectorial de polinomios en y de grado ≤ 2
- **T₂(xy)**: Espacio tensorial construido como T₂(xy) = P₂(x) ⊗ G₂(y)

### 2. Bases Consideradas

#### Base de Monomios
- **P₂(x)**: {1, x, x²}
- **G₂(y)**: {1, y, y²}

#### Base de Polinomios de Legendre
- **P₀(x) = 1**
- **P₁(x) = x**
- **P₂(x) = (3x² - 1)/2**

### 3. Producto Tensorial

El producto tensorial de dos polinomios p(x) y q(y) se define como:
```
p(x) ⊗ q(y) = p(x) · q(y)
```

## Solución Paso a Paso

### Parte (a): Expresión en Base de Legendre

**Objetivo**: Expresar p^P(x) = x² + x + 3 en términos de la base de polinomios de Legendre.

**Proceso**:
1. Tenemos el polinomio: p^P(x) = x² + x + 3
2. Necesitamos encontrar coeficientes α₀, α₁, α₂ tales que:
   ```
   x² + x + 3 = α₀P₀(x) + α₁P₁(x) + α₂P₂(x)
   ```

3. Sustituyendo las definiciones de los polinomios de Legendre:
   ```
   x² + x + 3 = α₀(1) + α₁(x) + α₂((3x² - 1)/2)
   ```

4. Resolviendo el sistema de ecuaciones:
   - Para x⁰: 3 = α₀ - α₂/2
   - Para x¹: 1 = α₁
   - Para x²: 1 = 3α₂/2

**Resultado**:
```
p^P(x) = (10/3)P₀(x) + P₁(x) + (2/3)P₂(x)
```

**Verificación**:
- P₀ = 1
- P₁ = x
- P₂ = (3x² - 1)/2

Reconstruyendo:
```
(10/3)(1) + (1)(x) + (2/3)((3x² - 1)/2) = 10/3 + x + (2/3)(3x²/2 - 1/2)
= 10/3 + x + x² - 1/3 = x² + x + 3 ✓
```

### Parte (b): Construcción del Tensor

**Objetivo**: Construir el tensor p^P⊗G(x,y) = p^P(x) ⊗ p^G(y)

**Datos**:
- p^P(x) = x² + x + 3
- p^G(y) = y + 1

**Proceso**:
```
p^P⊗G(x,y) = p^P(x) ⊗ p^G(y) = (x² + x + 3)(y + 1)
```

**Resultado**:
```
p^P⊗G(x,y) = x²y + xy + 3y + x² + x + 3
```

### Parte (c): Componentes en Base de Monomios

**Objetivo**: Encontrar las componentes c^ij del tensor en la base de monomios.

**Base**: {1, x, x²} ⊗ {1, y, y²}

**Proceso**:
El tensor p^P⊗G(x,y) = x²y + xy + 3y + x² + x + 3 se puede expresar como:
```
c^ij |e^P_i e^G_j⟩ = Σᵢⱼ c^ij x^i y^j
```

**Matriz de Componentes c^ij**:
```
     |  1  |  y  | y²
-----|-----|-----|-----
  1  |  3  |  3  |  0
  x  |  1  |  1  |  0
 x²  |  1  |  1  |  0
```

**Interpretación**:
- c^00 = 3: coeficiente del término 1·1 = 1
- c^01 = 3: coeficiente del término 1·y = y
- c^10 = 1: coeficiente del término x·1 = x
- c^11 = 1: coeficiente del término x·y = xy
- c^20 = 1: coeficiente del término x²·1 = x²
- c^21 = 1: coeficiente del término x²·y = x²y

### Parte (d): Componentes en Base de Legendre

**Objetivo**: Calcular las componentes c̃^ij del tensor en la base de polinomios de Legendre.

**Base**: {P₀, P₁, P₂} ⊗ {P₀, P₁, P₂}

**Proceso**:
1. Expresar cada término x^i y^j en términos de la base de Legendre
2. Aplicar las relaciones:
   - 1 = P₀
   - x = P₁
   - x² = (2P₂ + P₀)/3

**Matriz de Componentes c̃^ij**:
```
     |  P0  |  P1  |  P2
-----|------|------|-----
  P0 | 3.000| 3.000| 0.000
  P1 | 1.000| 1.000| 0.000
  P2 | 0.667| 0.667| 0.000
```

**Interpretación**:
- c̃^00 = 3: coeficiente del término P₀ ⊗ P₀
- c̃^01 = 3: coeficiente del término P₀ ⊗ P₁
- c̃^10 = 1: coeficiente del término P₁ ⊗ P₀
- c̃^11 = 1: coeficiente del término P₁ ⊗ P₁
- c̃^20 = 2/3: coeficiente del término P₂ ⊗ P₀
- c̃^21 = 2/3: coeficiente del término P₂ ⊗ P₁

## Implementación en Python

### Estructura del Código

1. **Clase PolinomioLegendre**: Maneja las operaciones con polinomios de Legendre
2. **Función producto_tensorial**: Calcula el producto tensorial de dos polinomios
3. **Función componentes_tensor_monomios**: Encuentra componentes en base de monomios
4. **Función componentes_tensor_legendre**: Encuentra componentes en base de Legendre
5. **Función main**: Ejecuta todo el proceso paso a paso

### Dependencias

```python
import numpy as np
import sympy as sp
from sympy import symbols, expand, simplify, collect
import matplotlib.pyplot as plt
```

### Características del Código

- **Cálculo simbólico**: Usa SymPy para manipulación exacta de polinomios
- **Visualización**: Genera gráficos de las matrices de componentes
- **Verificación**: Incluye verificaciones de los resultados
- **Documentación**: Comentarios extensos explicando cada paso

## Resultados y Verificaciones

### Verificación de la Parte (a)
El polinomio reconstruido a partir de los coeficientes de Legendre debe coincidir con el original:
```
(10/3)P₀ + P₁ + (2/3)P₂ = (10/3)(1) + (x) + (2/3)((3x²-1)/2) = x² + x + 3 ✓
```

### Verificación del Tensor
El producto tensorial debe expandirse correctamente:
```
(x² + x + 3)(y + 1) = x²y + xy + 3y + x² + x + 3 ✓
```

### Comparación de Bases
Las dos representaciones (monomios vs Legendre) deben representar el mismo tensor, solo que en bases diferentes.

## Aplicaciones en Física

Como se menciona en el problema, estos conceptos tienen aplicaciones importantes en:

1. **Esfuerzos**: En mecánica de sólidos, los tensores de esfuerzo se pueden representar en diferentes bases
2. **Inercia**: Los tensores de inercia en mecánica clásica
3. **Energía libre**: En termodinámica estadística, los desarrollos en polinomios ortogonales

## Conclusiones

1. **Flexibilidad de bases**: El mismo tensor puede representarse en diferentes bases, cada una con sus ventajas
2. **Polinomios ortogonales**: Los polinomios de Legendre ofrecen propiedades de ortogonalidad útiles
3. **Productos tensoriales**: Permiten construir espacios de mayor dimensión a partir de espacios más simples
4. **Implementación computacional**: Python con SymPy permite manipulación exacta y simbólica

## Archivos Generados

- `ejercicio_3_3_5.py`: Código principal con la implementación completa
- `tensor_components.png`: Gráfico comparativo de las matrices de componentes
- `Taller_3_Explicacion.md`: Esta documentación explicativa

## Ejecución del Código

Para ejecutar el código:

```bash
python ejercicio_3_3_5.py
```

El programa mostrará:
1. Los resultados de cada parte del problema
2. Las matrices de componentes en ambas bases
3. Un gráfico comparativo
4. Un resumen de todos los resultados

## Resultados Finales

### Parte (a): Expresión en Base de Legendre
```
p^P(x) = x² + x + 3 = (10/3)P₀ + P₁ + (2/3)P₂
```

### Parte (b): Tensor
```
p^P⊗G(x,y) = x²y + xy + 3y + x² + x + 3
```

### Parte (c): Componentes en Monomios
```
Matriz 3×3 con elementos c^ij en base {1, x, x²} ⊗ {1, y, y²}
```

### Parte (d): Componentes en Legendre
```
Matriz 3×3 con elementos c̃^ij en base {P₀, P₁, P₂} ⊗ {P₀, P₁, P₂}
```

El código está completamente funcional y verificado, proporcionando una implementación robusta del problema de productos tensoriales de espacios vectoriales de polinomios.