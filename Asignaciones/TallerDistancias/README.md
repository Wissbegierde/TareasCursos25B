# Taller de Distancias - Calibración de Sensores IoT

Este proyecto implementa métodos matemáticos para la calibración de sensores IoT de bajo costo que miden concentración de PM2.5, utilizando distancias euclidianas y promedios móviles.

## Descripción del Problema

Los sensores IoT de bajo costo son parte de la revolución del "Internet de las Cosas" pero a menudo no son lo suficientemente precisos y necesitan ser calibrados usando un estándar de referencia. Este taller demuestra que esta calibración está íntimamente ligada al concepto de "métrica".

### Objetivos

1. **Cuantificar el error de medición** del sensor de bajo costo
2. **Calibrar el sensor** para establecer lecturas más precisas
3. **Determinar el alcance de validez** del modelo de calibración
4. **Optimizar parámetros** como el ancho de ventana para promedios móviles

## Metodología

### 1. Cálculo de Distancia Euclidiana

Se calcula la distancia euclidiana entre mediciones de estaciones de referencia y sensores IoT:

```
D(D_i, D̂_i) = √(∑_{i,î} (D_i - D̂_i)²)
```

Donde:
- `D_i = {(x_1, y_1), (x_2, y_2) ... (x_n, y_n)}` es el conjunto de datos de referencia
- `D̂_i = {(x̂_1, ŷ_1), (x̂_2, ŷ_2) ... (x̂_m, ŷ_m)}` es el conjunto de datos a calibrar
- `f(x_i) = y_i` representa la variable dependiente medida por la estación de referencia
- `f(x̂_i) = ŷ_i` representa las mediciones de la estación a calibrar

### 2. Promedios Móviles

Se utiliza el criterio de promedio móvil para identificar datos "más cercanos":
- Comparar promedios locales de ambos conjuntos: `f(ξ_j)` y `f(ξ̂_j)`
- Calculados para una ventana común: `a_j ≤ x_i, x̂_i ≤ b_j`
- Elección posible: `ξ_j = a_j + (b_j - a_j)/2`

### 3. Ajuste Lineal

Se realiza un ajuste por mínimos cuadrados para determinar el modelo:
```
f(ξ̂_j) = αf(ξ_j)
```

Si ambas estaciones midieran lo mismo, `α` sería 1.

### 4. Alcance de Validez

Se determina el alcance de validez del modelo lineal definiendo una tolerancia y encontrando el rango en `x̂_i` para la validez del modelo.

## Instalación

1. Instalar las dependencias:
```bash
pip install -r requirements.txt
```

2. Ejecutar el taller:
```bash
python taller_distancias.py
```

## Uso

### Ejemplo Básico

```python
from taller_distancias import CalibradorSensores, generar_datos_ejemplo
import numpy as np

# Generar datos de ejemplo
datos_ref, datos_sensor = generar_datos_ejemplo(n_puntos=200, ruido=0.15)

# Crear calibrador
calibrador = CalibradorSensores(datos_ref, datos_sensor)

# Encontrar mejor ventana
ventanas = np.linspace(0.5, 5.0, 20)
mejor_ventana, menor_distancia, distancias = calibrador.encontrar_mejor_ventana(ventanas)

# Ajuste lineal
modelo = calibrador.ajuste_lineal(mejor_ventana)

# Alcance de validez
alcance = calibrador.determinar_alcance_validez(tolerancia=0.15)

# Visualizaciones
calibrador.visualizar_resultados(mejor_ventana)

# Reporte
reporte = calibrador.generar_reporte(mejor_ventana)
print(reporte)
```

### Con Datos Reales

```python
import pandas as pd
import numpy as np
from taller_distancias import CalibradorSensores

# Cargar datos de referencia
df_ref = pd.read_csv('Datos Estaciones AMB.csv')
tiempos_ref = df_ref['tiempo'].values
concentraciones_ref = df_ref['PM2.5'].values
datos_ref = np.column_stack([tiempos_ref, concentraciones_ref])

# Cargar datos del sensor
df_sensor = pd.read_csv('mediciones_sensor.csv')
tiempos_sensor = df_sensor['tiempo'].values
concentraciones_sensor = df_sensor['PM2.5'].values
datos_sensor = np.column_stack([tiempos_sensor, concentraciones_sensor])

# Crear calibrador y realizar análisis
calibrador = CalibradorSensores(datos_ref, datos_sensor)
# ... resto del análisis
```

## Archivos del Proyecto

- `taller_distancias.py`: Implementación principal del calibrador
- `ejemplo_uso.py`: Ejemplos de uso y análisis de sensibilidad
- `requirements.txt`: Dependencias del proyecto
- `README.md`: Este archivo de documentación

## Resultados

El análisis genera:

1. **Gráficos de visualización**:
   - Datos originales vs promedios móviles
   - Ajuste lineal entre estación de referencia y sensor
   - Optimización de ventana móvil

2. **Métricas de calibración**:
   - Distancia euclidiana óptima
   - Parámetros del modelo lineal (pendiente, intercepto)
   - R² y RMSE del ajuste
   - Porcentaje de puntos válidos para una tolerancia dada

3. **Reporte completo** en formato texto

## Datos de Referencia

Los datos de referencia se pueden obtener del repositorio:
https://github.com/nunezluis/MisCursos/tree/main/MisMateriales/Asignaciones/TallerDistancias/DatosDistancias

- `Datos Estaciones AMB`: Mediciones de referencia de concentración PM2.5
- `mediciones...`: Registros de las estaciones de bajo costo

## Referencias

- [Promedio Móvil](https://en.wikipedia.org/wiki/Moving_average)
- [PM2.5](https://blissair.com/what-is-pm-2-5.htm)
- Taller de Luis A. Núñez - Escuela de Física

## Autor

Implementación basada en el taller de **Luis A. Núñez** - Escuela de Física
