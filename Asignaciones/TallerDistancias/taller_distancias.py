"""
Taller de Distancias - Calibración de Sensores IoT para PM2.5
Métodos Matemáticos - Luis A. Núñez

Este módulo implementa:
1. Cálculo de distancia euclidiana entre conjuntos de datos
2. Promedios móviles con ventanas variables
3. Ajuste lineal por mínimos cuadrados
4. Determinación del alcance de validez del modelo
5. Visualizaciones con Matplotlib

Autor: Implementación basada en el taller de Luis A. Núñez
"""

import numpy as np
import matplotlib.pyplot as plt
import sympy as sp
from scipy import stats
from scipy.optimize import minimize
import pandas as pd
from typing import Tuple, List, Dict, Optional
import warnings
import os
warnings.filterwarnings('ignore')

# Configuración de matplotlib para mejor visualización
plt.rcParams['figure.figsize'] = (12, 8)
plt.rcParams['font.size'] = 10
plt.rcParams['axes.grid'] = True
plt.rcParams['grid.alpha'] = 0.3

class CalibradorSensores:
    """
    Clase principal para la calibración de sensores IoT usando distancias euclidianas
    y promedios móviles.
    """
    
    def __init__(self, datos_referencia: np.ndarray, datos_sensor: np.ndarray):
        """
        Inicializa el calibrador con los datos de referencia y del sensor.
        
        Args:
            datos_referencia: Array de forma (n, 2) con [tiempo, concentracion_PM25]
            datos_sensor: Array de forma (m, 2) con [tiempo, concentracion_PM25]
        """
        self.datos_referencia = np.array(datos_referencia)
        self.datos_sensor = np.array(datos_sensor)
        self.modelo_lineal = None
        self.parametros_modelo = None
        self.tolerancia = None
        
    def distancia_euclidiana(self, ventana: float = 1.0) -> float:
        """
        Calcula la distancia euclidiana entre los conjuntos de datos usando promedios móviles.
        
        Args:
            ventana: Ancho de la ventana para el promedio móvil (en horas)
            
        Returns:
            Distancia euclidiana D(D_i, D̂_i)
        """
        # Calcular promedios móviles
        promedios_ref = self._calcular_promedio_movil(
            self.datos_referencia, ventana
        )
        promedios_sensor = self._calcular_promedio_movil(
            self.datos_sensor, ventana
        )
        
        # Encontrar puntos más cercanos en el tiempo
        distancias = []
        for t_ref, y_ref in promedios_ref:
            # Buscar el punto más cercano en el tiempo
            tiempos_sensor = promedios_sensor[:, 0]
            idx_cercano = np.argmin(np.abs(tiempos_sensor - t_ref))
            y_sensor = promedios_sensor[idx_cercano, 1]
            
            # Calcular distancia euclidiana
            distancia = np.sqrt((y_ref - y_sensor)**2)
            distancias.append(distancia)
        
        return np.mean(distancias)
    
    def _calcular_promedio_movil(self, datos: np.ndarray, ventana: float) -> np.ndarray:
        """
        Calcula el promedio móvil de los datos.
        
        Args:
            datos: Array de forma (n, 2) con [tiempo, valor]
            ventana: Ancho de la ventana en las mismas unidades que el tiempo
            
        Returns:
            Array con [tiempo_centro, promedio]
        """
        tiempos, valores = datos[:, 0], datos[:, 1]
        promedios = []
        
        for i, t in enumerate(tiempos):
            # Definir ventana alrededor del punto actual
            t_inicio = t - ventana/2
            t_fin = t + ventana/2
            
            # Encontrar puntos dentro de la ventana
            mascara = (tiempos >= t_inicio) & (tiempos <= t_fin)
            valores_ventana = valores[mascara]
            
            if len(valores_ventana) > 0:
                promedio = np.mean(valores_ventana)
                promedios.append([t, promedio])
        
        return np.array(promedios)
    
    def encontrar_mejor_ventana(self, ventanas: List[float]) -> Tuple[float, float, List[float]]:
        """
        Encuentra el mejor ancho de ventana para minimizar la distancia euclidiana.
        
        Args:
            ventanas: Lista de anchos de ventana a probar
            
        Returns:
            Tupla con (mejor_ventana, menor_distancia, distancias)
        """
        distancias = []
        
        for ventana in ventanas:
            try:
                dist = self.distancia_euclidiana(ventana)
                distancias.append(dist)
                print(f"Ventana: {ventana:.2f}h, Distancia: {dist:.4f}")
            except Exception as e:
                print(f"Error con ventana {ventana}: {e}")
                distancias.append(np.inf)
        
        # Encontrar la mejor ventana
        idx_mejor = np.argmin(distancias)
        mejor_ventana = ventanas[idx_mejor]
        menor_distancia = distancias[idx_mejor]
        
        return mejor_ventana, menor_distancia, distancias
    
    def ajuste_lineal(self, ventana: float = 1.0) -> Dict:
        """
        Realiza ajuste lineal por mínimos cuadrados entre los promedios móviles.
        
        Args:
            ventana: Ancho de la ventana para el promedio móvil
            
        Returns:
            Diccionario con parámetros del modelo y estadísticas
        """
        # Calcular promedios móviles
        promedios_ref = self._calcular_promedio_movil(self.datos_referencia, ventana)
        promedios_sensor = self._calcular_promedio_movil(self.datos_sensor, ventana)
        
        # Preparar datos para ajuste
        x_vals, y_vals = [], []
        for t_ref, y_ref in promedios_ref:
            tiempos_sensor = promedios_sensor[:, 0]
            idx_cercano = np.argmin(np.abs(tiempos_sensor - t_ref))
            y_sensor = promedios_sensor[idx_cercano, 1]
            
            x_vals.append(y_ref)
            y_vals.append(y_sensor)
        
        x_vals = np.array(x_vals)
        y_vals = np.array(y_vals)
        
        # Ajuste lineal
        slope, intercept, r_value, p_value, std_err = stats.linregress(x_vals, y_vals)
        
        # Calcular predicciones
        y_pred = slope * x_vals + intercept #y = mx + b
        residuos = y_vals - y_pred 
        mse = np.mean(residuos**2)
        rmse = np.sqrt(mse)
        
        self.modelo_lineal = {
            'pendiente': slope,
            'intercepto': intercept,
            'r_cuadrado': r_value**2,
            'p_valor': p_value,
            'error_estandar': std_err,
            'rmse': rmse,
            'mse': mse
        }
        
        return self.modelo_lineal
    
    def determinar_alcance_validez(self, tolerancia: float = 0.1) -> Dict:
        """
        Determina el alcance de validez del modelo lineal para una tolerancia dada.
        
        Args:
            tolerancia: Error máximo aceptable (como fracción)
            
        Returns:
            Diccionario con información del alcance de validez
        """
        if self.modelo_lineal is None:
            raise ValueError("Debe ejecutar ajuste_lineal() primero")
        
        # Usar la mitad de los datos para el modelo
        n_ref = len(self.datos_referencia)
        n_sensor = len(self.datos_sensor)
        
        mitad_ref = n_ref // 2
        mitad_sensor = n_sensor // 2
        
        # Datos para entrenamiento
        datos_ref_train = self.datos_referencia[:mitad_ref]
        datos_sensor_train = self.datos_sensor[:mitad_sensor]
        
        # Datos para validación
        datos_ref_val = self.datos_referencia[mitad_ref:]
        datos_sensor_val = self.datos_sensor[mitad_sensor:]
        
        # Crear calibrador temporal para entrenamiento
        calibrador_temp = CalibradorSensores(datos_ref_train, datos_sensor_train)
        modelo_temp = calibrador_temp.ajuste_lineal()
        
        # Validar con datos de prueba
        promedios_ref_val = calibrador_temp._calcular_promedio_movil(datos_ref_val, 1.0)
        promedios_sensor_val = calibrador_temp._calcular_promedio_movil(datos_sensor_val, 1.0)
        
        errores = []
        for t_ref, y_ref in promedios_ref_val:
            tiempos_sensor = promedios_sensor_val[:, 0]
            if len(tiempos_sensor) > 0:
                idx_cercano = np.argmin(np.abs(tiempos_sensor - t_ref))
                y_sensor = promedios_sensor_val[idx_cercano, 1]
                
                # Predicción del modelo
                y_pred = modelo_temp['pendiente'] * y_ref + modelo_temp['intercepto']
                error_relativo = abs(y_sensor - y_pred) / y_ref
                errores.append(error_relativo)
        
        errores = np.array(errores)
        puntos_validos = np.sum(errores <= tolerancia)
        porcentaje_valido = (puntos_validos / len(errores)) * 100
        
        self.tolerancia = tolerancia
        
        return {
            'tolerancia': tolerancia,
            'puntos_validos': puntos_validos,
            'total_puntos': len(errores),
            'porcentaje_valido': porcentaje_valido,
            'error_promedio': np.mean(errores),
            'error_maximo': np.max(errores),
            'errores': errores
        }
    
    def visualizar_resultados(self, ventana: float = 1.0, guardar: bool = True):
        """
        Crea visualizaciones de los resultados del análisis.
        
        Args:
            ventana: Ancho de la ventana para visualización
            guardar: Si guardar las figuras como archivos
        """
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        fig.suptitle('Análisis de Calibración de Sensores IoT - PM2.5', fontsize=16, fontweight='bold')
        
        # 1. Datos originales
        ax1 = axes[0, 0]
        ax1.plot(self.datos_referencia[:, 0], self.datos_referencia[:, 1], 
                'b-', label='Estación Referencia', alpha=0.7, linewidth=2)
        ax1.plot(self.datos_sensor[:, 0], self.datos_sensor[:, 1], 
                'r-', label='Sensor IoT', alpha=0.7, linewidth=2)
        ax1.set_xlabel('Tiempo (horas)')
        ax1.set_ylabel('Concentración PM2.5 (μg/m³)')
        ax1.set_title('Datos Originales')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        
        # 2. Promedios móviles
        ax2 = axes[0, 1]
        promedios_ref = self._calcular_promedio_movil(self.datos_referencia, ventana)
        promedios_sensor = self._calcular_promedio_movil(self.datos_sensor, ventana)
        
        ax2.plot(promedios_ref[:, 0], promedios_ref[:, 1], 
                'b-o', label='Ref. Promedio Móvil', markersize=4)
        ax2.plot(promedios_sensor[:, 0], promedios_sensor[:, 1], 
                'r-s', label='Sensor Promedio Móvil', markersize=4)
        ax2.set_xlabel('Tiempo (horas)')
        ax2.set_ylabel('Concentración PM2.5 (μg/m³)')
        ax2.set_title(f'Promedios Móviles (Ventana: {ventana}h)')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        
        # 3. Ajuste lineal
        ax3 = axes[1, 0]
        if self.modelo_lineal is not None:
            # Preparar datos para el scatter plot
            x_vals, y_vals = [], []
            for t_ref, y_ref in promedios_ref:
                tiempos_sensor = promedios_sensor[:, 0]
                if len(tiempos_sensor) > 0:
                    idx_cercano = np.argmin(np.abs(tiempos_sensor - t_ref))
                    y_sensor = promedios_sensor[idx_cercano, 1]
                    x_vals.append(y_ref)
                    y_vals.append(y_sensor)
            
            x_vals = np.array(x_vals)
            y_vals = np.array(y_vals)
            
            # Scatter plot
            ax3.scatter(x_vals, y_vals, alpha=0.6, s=30)
            
            # Línea de ajuste
            x_line = np.linspace(x_vals.min(), x_vals.max(), 100)
            y_line = self.modelo_lineal['pendiente'] * x_line + self.modelo_lineal['intercepto']
            ax3.plot(x_line, y_line, 'r-', linewidth=2, 
                    label=f'y = {self.modelo_lineal["pendiente"]:.3f}x + {self.modelo_lineal["intercepto"]:.3f}')
            
            # Línea de referencia (y = x)
            ax3.plot(x_line, x_line, 'k--', alpha=0.5, label='y = x (referencia)')
            
            ax3.set_xlabel('Concentración Referencia (μg/m³)')
            ax3.set_ylabel('Concentración Sensor (μg/m³)')
            ax3.set_title(f'Ajuste Lineal (R² = {self.modelo_lineal["r_cuadrado"]:.3f})')
            ax3.legend()
            ax3.grid(True, alpha=0.3)
        
        # 4. Análisis de errores
        ax4 = axes[1, 1]
        if self.tolerancia is not None:
            # Simular análisis de errores para visualización
            ventanas = np.linspace(0.5, 5.0, 20)
            distancias = [self.distancia_euclidiana(w) for w in ventanas]
            
            ax4.plot(ventanas, distancias, 'g-o', linewidth=2, markersize=6)
            ax4.set_xlabel('Ancho de Ventana (horas)')
            ax4.set_ylabel('Distancia Euclidiana')
            ax4.set_title('Optimización de Ventana Móvil')
            ax4.grid(True, alpha=0.3)
            
            # Marcar el mínimo
            idx_min = np.argmin(distancias)
            ax4.plot(ventanas[idx_min], distancias[idx_min], 'ro', markersize=10, 
                    label=f'Mínimo: {ventanas[idx_min]:.2f}h')
            ax4.legend()
        
        plt.tight_layout()
        
        if guardar:
            # Crear directorio si no existe
            os.makedirs('Asignaciones/TallerDistancias', exist_ok=True)
            plt.savefig('Asignaciones/TallerDistancias/resultados_calibracion.png', 
                       dpi=300, bbox_inches='tight')
            print("Gráfico guardado como 'resultados_calibracion.png'")
        
        plt.show()
    
    def generar_reporte(self, ventana: float = 1.0) -> str:
        """
        Genera un reporte completo del análisis de calibración.
        
        Args:
            ventana: Ancho de la ventana utilizada
            
        Returns:
            String con el reporte formateado
        """
        reporte = []
        reporte.append("="*60)
        reporte.append("REPORTE DE CALIBRACIÓN DE SENSORES IoT - PM2.5")
        reporte.append("="*60)
        reporte.append("")
        
        # Información de los datos
        reporte.append("INFORMACIÓN DE DATOS:")
        reporte.append(f"- Puntos de referencia: {len(self.datos_referencia)}")
        reporte.append(f"- Puntos del sensor: {len(self.datos_sensor)}")
        reporte.append(f"- Rango temporal referencia: {self.datos_referencia[:, 0].min():.2f} - {self.datos_referencia[:, 0].max():.2f} horas")
        reporte.append(f"- Rango temporal sensor: {self.datos_sensor[:, 0].min():.2f} - {self.datos_sensor[:, 0].max():.2f} horas")
        reporte.append("")
        
        # Análisis de distancia euclidiana
        distancia = self.distancia_euclidiana(ventana)
        reporte.append("ANÁLISIS DE DISTANCIA EUCLIDIANA:")
        reporte.append(f"- Distancia euclidiana (ventana {ventana}h): {distancia:.4f}")
        reporte.append("")
        
        # Análisis de ventana óptima
        ventanas = np.linspace(0.5, 5.0, 20)
        mejor_ventana, menor_distancia, distancias = self.encontrar_mejor_ventana(ventanas)
        reporte.append("OPTIMIZACIÓN DE VENTANA MÓVIL:")
        reporte.append(f"- Mejor ventana: {mejor_ventana:.2f} horas")
        reporte.append(f"- Menor distancia: {menor_distancia:.4f}")
        reporte.append("")
        
        # Ajuste lineal
        if self.modelo_lineal is not None:
            reporte.append("AJUSTE LINEAL:")
            reporte.append(f"- Ecuación: y = {self.modelo_lineal['pendiente']:.4f}x + {self.modelo_lineal['intercepto']:.4f}")
            reporte.append(f"- R²: {self.modelo_lineal['r_cuadrado']:.4f}")
            reporte.append(f"- RMSE: {self.modelo_lineal['rmse']:.4f}")
            reporte.append(f"- P-valor: {self.modelo_lineal['p_valor']:.6f}")
            reporte.append("")
        
        # Alcance de validez
        if self.tolerancia is not None:
            alcance = self.determinar_alcance_validez(self.tolerancia)
            reporte.append("ALCANCE DE VALIDEZ:")
            reporte.append(f"- Tolerancia: {self.tolerancia:.1%}")
            reporte.append(f"- Puntos válidos: {alcance['puntos_validos']}/{alcance['total_puntos']}")
            reporte.append(f"- Porcentaje válido: {alcance['porcentaje_valido']:.1f}%")
            reporte.append(f"- Error promedio: {alcance['error_promedio']:.4f}")
            reporte.append(f"- Error máximo: {alcance['error_maximo']:.4f}")
            reporte.append("")
        
        reporte.append("="*60)
        
        return "\n".join(reporte)


def generar_datos_ejemplo(n_puntos: int = 100, ruido: float = 0.1) -> Tuple[np.ndarray, np.ndarray]:
    """
    Genera datos de ejemplo para probar el calibrador.
    
    Args:
        n_puntos: Número de puntos a generar
        ruido: Nivel de ruido (desviación estándar)
        
    Returns:
        Tupla con (datos_referencia, datos_sensor)
    """
    # Generar tiempos
    tiempos = np.linspace(0, 24, n_puntos)  # 24 horas
    
    # Generar datos de referencia (tendencia + estacionalidad + ruido)
    tendencia = 20 + 5 * np.sin(2 * np.pi * tiempos / 24)  # Ciclo diario
    estacionalidad = 10 * np.sin(2 * np.pi * tiempos / 12)  # Ciclo de 12 horas
    ruido_ref = np.random.normal(0, ruido * 5, n_puntos)
    datos_ref = tendencia + estacionalidad + ruido_ref
    
    # Generar datos del sensor (con calibración incorrecta)
    factor_calibracion = 0.8  # Sensor subestima
    offset = 2.0  # Offset constante
    ruido_sensor = np.random.normal(0, ruido * 8, n_puntos)
    datos_sensor = factor_calibracion * datos_ref + offset + ruido_sensor
    
    # Crear arrays con tiempo y concentración
    datos_referencia = np.column_stack([tiempos, datos_ref])
    datos_sensor = np.column_stack([tiempos, datos_sensor])
    
    return datos_referencia, datos_sensor


def main():
    """
    Función principal para ejecutar el taller de distancias.
    """
    print("Taller de Distancias - Calibración de Sensores IoT")
    print("=" * 50)
    
    # Generar datos de ejemplo
    print("Generando datos de ejemplo...")
    datos_ref, datos_sensor = generar_datos_ejemplo(n_puntos=200, ruido=0.15)
    
    # Crear calibrador
    calibrador = CalibradorSensores(datos_ref, datos_sensor)
    
    # Análisis de distancia euclidiana
    print("\n1. Análisis de Distancia Euclidiana:")
    ventanas = np.linspace(0.5, 5.0, 15)
    mejor_ventana, menor_distancia, distancias = calibrador.encontrar_mejor_ventana(ventanas)
    
    # Ajuste lineal
    print(f"\n2. Ajuste Lineal (ventana óptima: {mejor_ventana:.2f}h):")
    modelo = calibrador.ajuste_lineal(mejor_ventana)
    print(f"   Ecuación: y = {modelo['pendiente']:.4f}x + {modelo['intercepto']:.4f}")
    print(f"   R²: {modelo['r_cuadrado']:.4f}")
    print(f"   RMSE: {modelo['rmse']:.4f}")
    
    # Alcance de validez
    print("\n3. Alcance de Validez:")
    tolerancia = 0.15  # 15% de error máximo
    alcance = calibrador.determinar_alcance_validez(tolerancia)
    print(f"   Tolerancia: {tolerancia:.1%}")
    print(f"   Puntos válidos: {alcance['puntos_validos']}/{alcance['total_puntos']}")
    print(f"   Porcentaje válido: {alcance['porcentaje_valido']:.1f}%")
    
    # Generar visualizaciones
    print("\n4. Generando visualizaciones...")
    calibrador.visualizar_resultados(mejor_ventana)
    
    # Generar reporte
    print("\n5. Generando reporte...")
    reporte = calibrador.generar_reporte(mejor_ventana)
    print(reporte)
    
    # Guardar reporte
    os.makedirs('Asignaciones/TallerDistancias', exist_ok=True)
    with open('Asignaciones/TallerDistancias/reporte_calibracion.txt', 'w', encoding='utf-8') as f:
        f.write(reporte)
    print("Reporte guardado como 'reporte_calibracion.txt'")


if __name__ == "__main__":
    main()
