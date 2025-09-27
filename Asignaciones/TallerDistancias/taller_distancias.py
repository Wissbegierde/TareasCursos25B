"""
Taller de Distancias - Análisis con Datos Reales
Calibración de Sensores IoT para PM2.5 usando datos reales

Este script:
1. Carga datos reales del sensor IoT desde la carpeta DatosDistancias
2. Genera datos de referencia simulados
3. Realiza análisis de calibración completo
4. Genera visualizaciones y reportes

Autor: Basado en el taller de Luis A. Núñez
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy import stats
from datetime import datetime
import os
import sys
import atexit
import warnings
warnings.filterwarnings('ignore')

# Configuración de matplotlib
plt.rcParams['figure.figsize'] = (15, 12)
plt.rcParams['font.size'] = 10
plt.rcParams['axes.grid'] = True
plt.rcParams['grid.alpha'] = 0.3

# Utilidad para duplicar stdout/stderr a un archivo
class TeeStream:
    def __init__(self, stream_a, stream_b):
        self.stream_a = stream_a
        self.stream_b = stream_b
        self.encoding = getattr(stream_a, "encoding", "utf-8")

    def write(self, data):
        self.stream_a.write(data)
        self.stream_b.write(data)
        try:
            self.stream_a.flush()
            self.stream_b.flush()
        except Exception:
            pass

    def flush(self):
        try:
            self.stream_a.flush()
            self.stream_b.flush()
        except Exception:
            pass

    def isatty(self):
        try:
            return self.stream_a.isatty()
        except Exception:
            return False

    def fileno(self):
        try:
            return self.stream_a.fileno()
        except Exception:
            return -1

# Variables y limpieza de logging a archivo
_LOG_FILE = None
_OLD_STDOUT = None
_OLD_STDERR = None

def _cleanup_logging():
    global _LOG_FILE, _OLD_STDOUT, _OLD_STDERR
    try:
        if _OLD_STDOUT is not None:
            sys.stdout = _OLD_STDOUT
        if _OLD_STDERR is not None:
            sys.stderr = _OLD_STDERR
        if _LOG_FILE is not None:
            try:
                _LOG_FILE.flush()
            except Exception:
                pass
            try:
                _LOG_FILE.close()
            except Exception:
                pass
    finally:
        _LOG_FILE = None
        _OLD_STDOUT = None
        _OLD_STDERR = None

class CalibradorSensores:
    """
    Clase para la calibración de sensores IoT usando distancias euclidianas.
    """
    
    def __init__(self, datos_referencia, datos_sensor):
        """
        Inicializa el calibrador con los datos de referencia y del sensor.
        
        Args:
            datos_referencia: Array de forma (n, 2) con [tiempo, concentracion_PM25]
            datos_sensor: Array de forma (m, 2) con [tiempo, concentracion_PM25]
        """
        self.datos_referencia = np.array(datos_referencia)
        self.datos_sensor = np.array(datos_sensor)
        self.modelo_lineal = None
        self.tolerancia = None
        
    def distancia_euclidiana(self, ventana=1.0):
        """
        Calcula la distancia euclidiana entre los conjuntos de datos usando promedios móviles.
        
        Args:
            ventana: Ancho de la ventana para el promedio móvil (en horas)
            
        Returns:
            Distancia euclidiana D(D_i, D̂_i)
        """
        # Calcular promedios móviles
        promedios_ref = self._calcular_promedio_movil(self.datos_referencia, ventana)
        promedios_sensor = self._calcular_promedio_movil(self.datos_sensor, ventana)
        
        # Usar pares sincronizados (por tiempos del sensor) con tolerancia estricta
        x_vals, y_vals = self._pares_sincronizados(promedios_ref, promedios_sensor, tolerancia_horas=0.25)
        if x_vals.size == 0:
            return float('inf')
        # Distancia euclidiana promedio (RMSE)
        return float(np.sqrt(np.mean((x_vals - y_vals) ** 2)))
    
    def _calcular_promedio_movil(self, datos, ventana):
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
    
    def encontrar_mejor_ventana(self, ventanas):
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
                print(f"  Ventana: {ventana:.2f}h, Distancia: {dist:.4f}")
            except Exception as e:
                print(f"  Error con ventana {ventana}: {e}")
                distancias.append(np.inf)
        
        # Encontrar la mejor ventana
        idx_mejor = np.argmin(distancias)
        mejor_ventana = ventanas[idx_mejor]
        menor_distancia = distancias[idx_mejor]
        
        return mejor_ventana, menor_distancia, distancias
    
    def ajuste_lineal(self, ventana=1.0):
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
        
        # Preparar datos sincronizados para ajuste (emparejar por tiempo)
        # LaTeX: f(ξ_j) = α \hat{f}(ξ_j) + β ⇒ referencia = f(sensor)
        ref_vals, sen_vals = self._pares_sincronizados(promedios_ref, promedios_sensor, tolerancia_horas=1.0)
        # x = sensor, y = referencia
        x_vals, y_vals = sen_vals, ref_vals
        
        x_vals = np.array(x_vals)
        y_vals = np.array(y_vals)
        
        # Ajuste lineal
        slope, intercept, r_value, p_value, std_err = stats.linregress(x_vals, y_vals)
        
        # Calcular predicciones
        y_pred = slope * x_vals + intercept
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

    def _pares_sincronizados(self, promedios_ref: np.ndarray, promedios_sensor: np.ndarray, tolerancia_horas: float = 1.0):
        """
        Construye pares (x, y) sincronizados por tiempo dentro de una tolerancia.
        x: referencia AMB, y: sensor IoT
        """
        x_vals: list[float] = []
        y_vals: list[float] = []
        if promedios_ref.size == 0 or promedios_sensor.size == 0:
            return np.array([]), np.array([])
        tiempos_sensor = promedios_sensor[:, 0]
        for t_ref, y_ref in promedios_ref:
            difs = np.abs(tiempos_sensor - t_ref)
            idxs = np.where(difs <= tolerancia_horas)[0]
            if idxs.size == 0:
                continue
            idx_cercano = idxs[np.argmin(difs[idxs])]
            y_sensor = promedios_sensor[idx_cercano, 1]
            x_vals.append(y_ref)
            y_vals.append(y_sensor)
        return np.array(x_vals), np.array(y_vals)

    def ajuste_mejor_modelo(self, ventana: float = 1.0):
        """
        Ajusta múltiples modelos (lineal, polinomial grado 2/3, log, exponencial, potencia)
        y selecciona el mejor por R² (y como desempate, menor RMSE).
        Devuelve un diccionario con el mejor modelo y guarda en self.modelo.
        """
        promedios_ref = self._calcular_promedio_movil(self.datos_referencia, ventana)
        promedios_sensor = self._calcular_promedio_movil(self.datos_sensor, ventana)
        x, y = self._pares_sincronizados(promedios_ref, promedios_sensor, tolerancia_horas=1.0)
        if x.size < 10 or y.size < 10:
            # Fallback a lineal si no hay suficientes puntos
            self.modelo = {
                'tipo': 'lineal',
                **self.ajuste_lineal(ventana)
            }
            return self.modelo

        modelos = []

        def r2_rmse(y_true, y_pred):
            ss_res = np.sum((y_true - y_pred) ** 2)
            ss_tot = np.sum((y_true - np.mean(y_true)) ** 2)
            r2 = 1.0 - (ss_res / ss_tot if ss_tot > 0 else np.inf)
            rmse = np.sqrt(np.mean((y_true - y_pred) ** 2))
            return r2, rmse

        # 1) Lineal: y = a x + b
        slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
        y_pred = slope * x + intercept
        r2, rmse = r2_rmse(y, y_pred)
        modelos.append({
            'tipo': 'lineal',
            'pendiente': slope,
            'intercepto': intercept,
            'r_cuadrado': r2,
            'rmse': rmse,
            'predict': lambda xq, a=slope, b=intercept: a * xq + b
        })

        # 2) Polinomial grado 2
        try:
            coef2 = np.polyfit(x, y, deg=2)
            y_pred2 = np.polyval(coef2, x)
            r2_2, rmse_2 = r2_rmse(y, y_pred2)
            modelos.append({
                'tipo': 'polinomio_2',
                'coeficientes': coef2,
                'r_cuadrado': r2_2,
                'rmse': rmse_2,
                'predict': lambda xq, c=coef2: np.polyval(c, xq)
            })
        except Exception:
            pass

        # 3) Polinomial grado 3
        try:
            coef3 = np.polyfit(x, y, deg=3)
            y_pred3 = np.polyval(coef3, x)
            r2_3, rmse_3 = r2_rmse(y, y_pred3)
            modelos.append({
                'tipo': 'polinomio_3',
                'coeficientes': coef3,
                'r_cuadrado': r2_3,
                'rmse': rmse_3,
                'predict': lambda xq, c=coef3: np.polyval(c, xq)
            })
        except Exception:
            pass

        # 4) Logarítmico: y = a ln(x) + b (x>0)
        try:
            mask = x > 0
            if np.sum(mask) > 10:
                Xl = np.log(x[mask])
                yl = y[mask]
                a, b, r, p, se = stats.linregress(Xl, yl)
                y_predl = a * Xl + b
                r2_l, rmse_l = r2_rmse(yl, y_predl)
                modelos.append({
                    'tipo': 'logaritmico',
                    'a': a,
                    'b': b,
                    'r_cuadrado': r2_l,
                    'rmse': rmse_l,
                    'predict': lambda xq, aa=a, bb=b: aa * np.log(xq) + bb
                })
        except Exception:
            pass

        # 5) Exponencial: y = A e^{B x} -> ln y = ln A + B x (y>0)
        try:
            mask = y > 0
            if np.sum(mask) > 10:
                Xe = x[mask]
                Ye = np.log(y[mask])
                B, lnA, r, p, se = stats.linregress(Xe, Ye)
                A = np.exp(lnA)
                y_prede = A * np.exp(B * Xe)
                r2_e, rmse_e = r2_rmse(np.exp(Ye), y_prede)
                modelos.append({
                    'tipo': 'exponencial',
                    'A': A,
                    'B': B,
                    'r_cuadrado': r2_e,
                    'rmse': rmse_e,
                    'predict': lambda xq, AA=A, BB=B: AA * np.exp(BB * xq)
                })
        except Exception:
            pass

        # 6) Potencia: y = A x^B -> ln y = ln A + B ln x (x>0, y>0)
        try:
            mask = (x > 0) & (y > 0)
            if np.sum(mask) > 10:
                Xp = np.log(x[mask])
                Yp = np.log(y[mask])
                B, lnA, r, p, se = stats.linregress(Xp, Yp)
                A = np.exp(lnA)
                y_predp = A * np.exp(B * Xp)
                r2_p, rmse_p = r2_rmse(np.exp(Yp), y_predp)
                modelos.append({
                    'tipo': 'potencia',
                    'A': A,
                    'B': B,
                    'r_cuadrado': r2_p,
                    'rmse': rmse_p,
                    'predict': lambda xq, AA=A, BB=B: AA * (xq ** BB)
                })
        except Exception:
            pass

        # Selección del mejor
        if not modelos:
            self.modelo = {'tipo': 'lineal', **self.ajuste_lineal(ventana)}
            return self.modelo

        # Mejor por R² y luego por menor RMSE
        modelos_ordenados = sorted(modelos, key=lambda m: (m['r_cuadrado'], -m['rmse']))
        mejor = modelos_ordenados[-1]
        self.modelo = mejor
        return mejor

    def _predecir_modelo(self, x_vals: np.ndarray) -> np.ndarray:
        """Devuelve y_pred para el modelo elegido si existe, o el lineal si no."""
        if hasattr(self, 'modelo') and self.modelo is not None and 'predict' in self.modelo and self.modelo.get('tipo') != 'lineal_forzado':
            # Evitar dominios inválidos para log/potencia
            try:
                return self.modelo['predict'](x_vals)
            except Exception:
                pass
        if self.modelo_lineal is not None:
            a = self.modelo_lineal['pendiente']
            b = self.modelo_lineal['intercepto']
            return a * x_vals + b
        return np.zeros_like(x_vals)

    def determinar_intervalo_validez(self, tolerancia: float = 0.25, ventana: float = 1.0):
        """
        Devuelve el intervalo en el eje del sensor (\hat{x}) donde el modelo cumple error relativo <= tolerancia.
        Usa pares sincronizados y el modelo actual (mejor o lineal).
        """
        promedios_ref = self._calcular_promedio_movil(self.datos_referencia, ventana)
        promedios_sensor = self._calcular_promedio_movil(self.datos_sensor, ventana)
        ref_vals, sen_vals = self._pares_sincronizados(promedios_ref, promedios_sensor, tolerancia_horas=1.0)
        if ref_vals.size == 0:
            return None
        # Predicción referencia a partir de sensor: y_pred = f(sensor)
        y_pred = self._predecir_modelo(sen_vals)
        err_rel = np.abs(ref_vals - y_pred) / np.maximum(ref_vals, 1e-6)
        mask = err_rel <= tolerancia
        if not np.any(mask):
            return None
        x_valid = sen_vals[mask]
        return float(np.min(x_valid)), float(np.max(x_valid)), float(np.mean(err_rel[mask])), float(np.max(err_rel[mask]))

    def buscar_tamano_minimo_datos(self, tolerancia: float = 0.25, ventana: float = 1.0, paso: int = 200):
        """
        Busca el tamaño mínimo de entrenamiento que logra intervalo de validez no vacío.
        Incrementa el tamaño de entrenamiento por pasos y evalúa en el resto.
        """
        promedios_ref = self._calcular_promedio_movil(self.datos_referencia, ventana)
        promedios_sensor = self._calcular_promedio_movil(self.datos_sensor, ventana)
        ref_vals, sen_vals = self._pares_sincronizados(promedios_ref, promedios_sensor, tolerancia_horas=1.0)
        n = ref_vals.size
        if n < 2 * paso:
            return None
        for n_train in range(paso, n - paso + 1, paso):
            x_train = sen_vals[:n_train]
            y_train = ref_vals[:n_train]
            x_val = sen_vals[n_train:]
            y_val = ref_vals[n_train:]
            # Ajuste lineal rápido sobre train para evaluar criterio
            try:
                a, b, r, p, se = stats.linregress(x_train, y_train)
                y_pred = a * x_val + b
                err_rel = np.abs(y_val - y_pred) / np.maximum(y_val, 1e-6)
                if np.any(err_rel <= tolerancia):
                    return {
                        'n_train': n_train,
                        'r2_train': float(r**2),
                        'intervalo_valido_val': (
                            float(np.min(x_val[err_rel <= tolerancia])),
                            float(np.max(x_val[err_rel <= tolerancia]))
                        )
                    }
            except Exception:
                continue
        return None
    
    def determinar_alcance_validez(self, tolerancia=0.1):
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
    
    def visualizar_resultados(self, ventana=1.0, guardar=True):
        """
        Crea visualizaciones de los resultados del análisis.
        
        Args:
            ventana: Ancho de la ventana para visualización
            guardar: Si guardar las figuras como archivos
        """
        fig, axes = plt.subplots(2, 2, figsize=(15, 12))
        fig.suptitle('Análisis de Calibración de Sensores IoT - PM2.5 (Datos Reales)', 
                     fontsize=16, fontweight='bold')
        
        # 1. Datos originales
        ax1 = axes[0, 0]
        ax1.plot(self.datos_referencia[:, 0], self.datos_referencia[:, 1], 
                'b-', label='Estaciones AMB (Referencia)', alpha=0.7, linewidth=2)
        ax1.plot(self.datos_sensor[:, 0], self.datos_sensor[:, 1], 
                'r-', label='Sensor IoT (A calibrar)', alpha=0.7, linewidth=2)
        ax1.set_xlabel('Tiempo (horas)')
        ax1.set_ylabel('Concentración PM2.5 (μg/m³)')
        ax1.set_title('Datos Originales - Comparación')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        
        # 2. Promedios móviles
        ax2 = axes[0, 1]
        promedios_ref = self._calcular_promedio_movil(self.datos_referencia, ventana)
        promedios_sensor = self._calcular_promedio_movil(self.datos_sensor, ventana)
        
        ax2.plot(promedios_ref[:, 0], promedios_ref[:, 1], 
                'b-o', label='AMB Promedio Móvil', markersize=4)
        ax2.plot(promedios_sensor[:, 0], promedios_sensor[:, 1], 
                'r-s', label='IoT Promedio Móvil', markersize=4)
        ax2.set_xlabel('Tiempo (horas)')
        ax2.set_ylabel('Concentración PM2.5 (μg/m³)')
        ax2.set_title(f'Promedios Móviles (Ventana: {ventana}h)')
        ax2.legend()
        ax2.grid(True, alpha=0.3)
        
        # 3. Ajuste (mejor modelo)
        ax3 = axes[1, 0]
        # Preparar datos sincronizados
        ref_scatter, sen_scatter = self._pares_sincronizados(promedios_ref, promedios_sensor, tolerancia_horas=1.0)
        if ref_scatter.size > 0:
            # Eje X: sensor, Eje Y: referencia
            ax3.scatter(sen_scatter, ref_scatter, alpha=0.6, s=30)
            x_line = np.linspace(sen_scatter.min(), sen_scatter.max(), 200)
            y_line = self._predecir_modelo(x_line)
            ax3.plot(x_line, y_line, 'r-', linewidth=2, label='Modelo ajustado')
            ax3.plot(x_line, x_line, 'k--', alpha=0.5, label='y = x (línea ideal)')
        ax3.set_xlabel('Concentración IoT (μg/m³)')
        ax3.set_ylabel('Concentración AMB (μg/m³)')
        ax3.set_title('Ajuste del Modelo (pares sincronizados)')
        ax3.legend()
        ax3.grid(True, alpha=0.3)
        
        # 4. Análisis de errores
        ax4 = axes[1, 1]
        # Simular análisis de errores para visualización
        ventanas = np.linspace(0.5, 6.0, 20)
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
            plt.savefig('Asignaciones/TallerDistancias/resultados_datos_reales.png', 
                       dpi=300, bbox_inches='tight')
            print("Gráfico guardado como 'resultados_datos_reales.png'")
        
        plt.show()
    
    def generar_reporte(self, ventana=1.0):
        """
        Genera un reporte completo del análisis de calibración.
        
        Args:
            ventana: Ancho de la ventana utilizada
            
        Returns:
            String con el reporte formateado
        """
        reporte = []
        reporte.append("="*70)
        reporte.append("REPORTE DE CALIBRACIÓN DE SENSORES IoT - PM2.5 (DATOS REALES)")
        reporte.append("="*70)
        reporte.append("")
        
        # Información de los datos
        reporte.append("INFORMACIÓN DE DATOS:")
        reporte.append(f"- Puntos de referencia: {len(self.datos_referencia)}")
        reporte.append(f"- Puntos del sensor: {len(self.datos_sensor)}")
        reporte.append(f"- Rango temporal referencia: {self.datos_referencia[:, 0].min():.2f} - {self.datos_referencia[:, 0].max():.2f} horas")
        reporte.append(f"- Rango temporal sensor: {self.datos_sensor[:, 0].min():.2f} - {self.datos_sensor[:, 0].max():.2f} horas")
        reporte.append(f"- Rango PM2.5 referencia: {self.datos_referencia[:, 1].min():.2f} - {self.datos_referencia[:, 1].max():.2f} μg/m³")
        reporte.append(f"- Rango PM2.5 sensor: {self.datos_sensor[:, 1].min():.2f} - {self.datos_sensor[:, 1].max():.2f} μg/m³")
        reporte.append("")
        
        # Análisis de distancia euclidiana
        distancia = self.distancia_euclidiana(ventana)
        reporte.append("ANÁLISIS DE DISTANCIA EUCLIDIANA:")
        reporte.append(f"- Distancia euclidiana (ventana {ventana}h): {distancia:.4f}")
        reporte.append("")
        
        # Análisis de ventana óptima
        ventanas = np.linspace(0.5, 6.0, 20)
        mejor_ventana, menor_distancia, distancias = self.encontrar_mejor_ventana(ventanas)
        reporte.append("OPTIMIZACIÓN DE VENTANA MÓVIL:")
        reporte.append(f"- Mejor ventana: {mejor_ventana:.2f} horas")
        reporte.append(f"- Menor distancia: {menor_distancia:.4f}")
        reporte.append("")
        
        # Mejor modelo (si existe)
        if hasattr(self, 'modelo') and self.modelo is not None:
            reporte.append("MEJOR MODELO:")
            reporte.append(f"- Tipo: {self.modelo.get('tipo','lineal')}")
            if self.modelo.get('tipo') == 'lineal' and self.modelo_lineal is not None:
                reporte.append(f"- Ecuación: y = {self.modelo_lineal['pendiente']:.4f}x + {self.modelo_lineal['intercepto']:.4f}")
            reporte.append(f"- R²: {self.modelo.get('r_cuadrado', np.nan):.4f}")
            reporte.append(f"- RMSE: {self.modelo.get('rmse', np.nan):.4f}")
            # Intervalo de validez (si puede calcularse)
            inter = self.determinar_intervalo_validez(self.tolerancia if self.tolerancia else 0.25, ventana)
            if inter:
                reporte.append(f"- Intervalo válido IoT: [{inter[0]:.2f}, {inter[1]:.2f}] μg/m³")
                reporte.append(f"- Error promedio (válidos): {inter[2]:.4f}")
                reporte.append(f"- Error máximo (válidos): {inter[3]:.4f}")
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
        
        reporte.append("="*70)
        
        return "\n".join(reporte)

def cargar_datos_referencia():
    """
    Carga datos de referencia desde el archivo Excel de estaciones AMB.
    
    Returns:
        DataFrame con datos de referencia procesados
    """
    print("Cargando datos de referencia (Estaciones AMB)...")
    
    archivo_excel = "DatosDistancias/Datos Estaciones AMB.xlsx"
    
    try:
        # Leer archivo Excel (saltar la primera fila que contiene nombres de columnas)
        df = pd.read_excel(archivo_excel, skiprows=1)
        print(f"Archivo Excel cargado: {df.shape}")
        print(f"Columnas disponibles: {list(df.columns)}")
        
        # Mostrar primeras filas para entender la estructura
        print("\nPrimeras 5 filas del archivo Excel:")
        print(df.head())
        
        # Procesar datos según la estructura real del Excel
        # Después de skiprows=1, las columnas son: Date&Time, ug/m3, ug/m3.1, etc.
        # La columna PM2.5 es la segunda columna (ug/m3.1)
        if 'Date&Time' in df.columns and len(df.columns) >= 3:
            # Usar las columnas reales del Excel
            df_procesado = df[['Date&Time', 'ug/m3.1']].copy()  # PM2.5 es la tercera columna
            df_procesado.columns = ['fecha', 'pm25']
            
            # Convertir fecha a datetime con zona UTC y filtrar 2019+
            df_procesado['fecha'] = pd.to_datetime(df_procesado['fecha'], utc=True)
            df_procesado = df_procesado[df_procesado['fecha'] >= pd.Timestamp('2019-01-01', tz='UTC')].copy()
            
            # Filtrar datos válidos (eliminar "NoData" y valores nulos)
            df_procesado = df_procesado[df_procesado['pm25'] != 'NoData'].copy()
            df_procesado['pm25'] = pd.to_numeric(df_procesado['pm25'], errors='coerce')
            df_procesado = df_procesado.dropna()
            
            # Seleccionar columnas necesarias (dejar fecha cruda, sin convertir a horas aún)
            df_procesado = df_procesado[['fecha', 'pm25']].copy()
        else:
            # Estructura desconocida, usar las primeras dos columnas
            print("Estructura desconocida, usando primeras dos columnas")
            df_procesado = df.iloc[:, [0, 1]].copy()
            df_procesado.columns = ['tiempo_horas', 'pm25']
        
        # Eliminar valores nulos
        df_procesado = df_procesado.dropna()
        
        print(f"\nDatos de referencia procesados: {len(df_procesado)} puntos")
        print(f"Rango temporal: {df_procesado['fecha'].min()} - {df_procesado['fecha'].max()}")
        print(f"Rango de PM2.5: {df_procesado['pm25'].min():.2f} - {df_procesado['pm25'].max():.2f} μg/m³")
        print(f"Promedio PM2.5: {df_procesado['pm25'].mean():.2f} μg/m³")
        
        return df_procesado
        
    except Exception as e:
        print(f"Error cargando datos de referencia: {e}")
        return None

def cargar_datos_sensor_iot():
    """
    Carga datos del sensor IoT desde los archivos CSV.
    
    Returns:
        DataFrame con datos del sensor IoT procesados
    """
    print("\nCargando datos del sensor IoT (archivos CSV)...")
    
    ruta_base = "DatosDistancias"
    archivos_csv = []
    
    # Encontrar todos los archivos CSV
    for archivo in os.listdir(ruta_base):
        if archivo.endswith('.csv') and 'mediciones_clg_normalsup_pm25' in archivo:
            archivos_csv.append(os.path.join(ruta_base, archivo))
    
    print(f"Encontrados {len(archivos_csv)} archivos CSV del sensor IoT")
    
    # Combinar datos de todos los archivos
    datos_combinados = []
    tiempo_inicio_global = None
    
    for archivo in archivos_csv:
        print(f"  Procesando: {os.path.basename(archivo)}")
        
        try:
            df = pd.read_csv(archivo)
            # Asegurar timezone UTC
            df['fecha_hora_med'] = pd.to_datetime(df['fecha_hora_med'], utc=True)
            df_pm25 = df[df['id_parametro'] == 'pm25_a'].copy()
            
            if len(df_pm25) > 0:
                # FILTRAR SOLO DATOS DE 2019 EN ADELANTE
                df_pm25 = df_pm25[df_pm25['fecha_hora_med'] >= pd.Timestamp('2019-01-01', tz='UTC')].copy()
                
                if len(df_pm25) > 0:
                    # Mantener fecha para sincronizar por tiempo real
                    df_pm25 = df_pm25[['fecha_hora_med', 'valor']].copy()
                    df_pm25 = df_pm25.rename(columns={'valor': 'pm25', 'fecha_hora_med': 'fecha'})
                    
                    datos_combinados.append(df_pm25)
                    print(f"    - {len(df_pm25)} puntos de PM2.5 (2019+)")
            
        except Exception as e:
            print(f"    - Error: {e}")
    
    if not datos_combinados:
        print("Error: No se pudieron cargar datos del sensor IoT")
        return None
    
    # Combinar todos los datos
    df_sensor = pd.concat(datos_combinados, ignore_index=True)
    df_sensor = df_sensor.sort_values('fecha').reset_index(drop=True)
    
    print(f"\nTotal de datos del sensor IoT: {len(df_sensor)} puntos")
    print(f"Rango temporal: {df_sensor['fecha'].min()} - {df_sensor['fecha'].max()}")
    print(f"Rango de PM2.5: {df_sensor['pm25'].min():.2f} - {df_sensor['pm25'].max():.2f} μg/m³")
    print(f"Promedio PM2.5: {df_sensor['pm25'].mean():.2f} μg/m³")
    
    return df_sensor

def generar_datos_referencia(tiempo_max, n_puntos=1000):
    """
    Genera datos de referencia simulados basados en patrones reales de PM2.5.
    
    Args:
        tiempo_max: Tiempo máximo en horas
        n_puntos: Número de puntos a generar
        
    Returns:
        Array con datos de referencia
    """
    print("Generando datos de referencia simulados...")
    
    # Crear tiempos
    tiempos = np.linspace(0, tiempo_max, n_puntos)
    
    # Simular patrones reales de PM2.5
    # Tendencia base
    tendencia_base = 15 + 5 * np.sin(2 * np.pi * tiempos / (24*7))  # Ciclo semanal
    
    # Estacionalidad diaria
    estacionalidad_diaria = 8 * np.sin(2 * np.pi * tiempos / 24)  # Ciclo diario
    
    # Picos ocasionales (simulando eventos de contaminación)
    picos = np.zeros_like(tiempos)
    for i in range(5):  # 5 picos aleatorios
        inicio_pico = np.random.uniform(0, tiempo_max-24)
        duracion = np.random.uniform(6, 24)
        intensidad = np.random.uniform(15, 30)
        
        mascara_pico = (tiempos >= inicio_pico) & (tiempos <= inicio_pico + duracion)
        picos[mascara_pico] += intensidad * np.exp(-(tiempos[mascara_pico] - inicio_pico) / 6)
    
    # Ruido
    ruido = np.random.normal(0, 2, len(tiempos))
    
    # Combinar componentes
    pm25_ref = tendencia_base + estacionalidad_diaria + picos + ruido
    
    # Asegurar valores positivos
    pm25_ref = np.maximum(pm25_ref, 1.0)
    
    datos_referencia = np.column_stack([tiempos, pm25_ref])
    
    print(f"Datos de referencia generados: {len(datos_referencia)} puntos")
    print(f"Rango de PM2.5: {pm25_ref.min():.2f} - {pm25_ref.max():.2f} μg/m³")
    print(f"Promedio: {pm25_ref.mean():.2f} μg/m³")
    
    return datos_referencia

def main():
    """
    Función principal para ejecutar el análisis completo con datos reales.
    """
    # Configurar duplicado de salida a archivo
    os.makedirs('Asignaciones/TallerDistancias', exist_ok=True)
    log_path = 'Asignaciones/TallerDistancias/reporte_datos_reales.txt'
    global _LOG_FILE, _OLD_STDOUT, _OLD_STDERR
    _LOG_FILE = open(log_path, 'w', encoding='utf-8')
    _OLD_STDOUT, _OLD_STDERR = sys.stdout, sys.stderr
    sys.stdout = TeeStream(_OLD_STDOUT, _LOG_FILE)
    sys.stderr = TeeStream(_OLD_STDERR, _LOG_FILE)
    atexit.register(_cleanup_logging)

    print("="*70)
    print("TALLER DE DISTANCIAS - ANÁLISIS CON DATOS REALES")
    print("="*70)
    
    # 1. Cargar datos de referencia (Estaciones AMB)
    df_referencia = cargar_datos_referencia()
    
    if df_referencia is None:
        print("Error: No se pudieron cargar datos de referencia")
        return
    
    # 2. Cargar datos del sensor IoT
    df_sensor = cargar_datos_sensor_iot()
    
    if df_sensor is None:
        print("Error: No se pudieron cargar datos del sensor IoT")
        return
    
    # 3. Preparar datos para calibración (sincronizados por tiempo real)
    print("\nPreparando datos para calibración (sincronizados por tiempo real)...")
    # Emparejar por timestamp cercano (tolerancia 1h)
    ref = df_referencia.copy()
    sen = df_sensor.copy()
    ref = ref.sort_values('fecha').reset_index(drop=True)
    sen = sen.sort_values('fecha').reset_index(drop=True)
    
    tolerancia = pd.Timedelta(hours=1)
    y_ref_sync = []
    y_sen_sync = []
    t0 = ref['fecha'].min()
    
    j = 0
    for i in range(len(ref)):
        tr = ref.loc[i, 'fecha']
        # avanzar j hasta estar cerca
        while j + 1 < len(sen) and abs(sen.loc[j+1, 'fecha'] - tr) < abs(sen.loc[j, 'fecha'] - tr):
            j += 1
        if abs(sen.loc[j, 'fecha'] - tr) <= tolerancia:
            y_ref_sync.append(ref.loc[i, 'pm25'])
            y_sen_sync.append(sen.loc[j, 'pm25'])
    
    # convertir a arreglo [tiempo_horas, pm25]
    tiempos_horas = ((ref['fecha'].loc[:len(y_ref_sync)-1] - t0).dt.total_seconds() / 3600).to_numpy()
    datos_referencia = np.column_stack([tiempos_horas, np.array(y_ref_sync)])
    datos_sensor = np.column_stack([tiempos_horas, np.array(y_sen_sync)])
    
    print(f"Datos de referencia (Estaciones AMB): {len(datos_referencia)} puntos")
    print(f"Datos del sensor IoT (a calibrar): {len(datos_sensor)} puntos")
    
    # 4. Crear calibrador y realizar análisis
    print("\n" + "="*50)
    print("REALIZANDO ANÁLISIS DE CALIBRACIÓN")
    print("="*50)
    
    calibrador = CalibradorSensores(datos_referencia, datos_sensor)
    
    # Análisis de distancia euclidiana
    print("\n1. Análisis de distancia euclidiana...")
    ventanas = np.linspace(0.5, 6.0, 20)
    mejor_ventana, menor_distancia, distancias = calibrador.encontrar_mejor_ventana(ventanas)
    
    print(f"\nMejor ventana: {mejor_ventana:.2f} horas")
    print(f"Menor distancia: {menor_distancia:.4f}")
    
    # Ajuste lineal forzado (ax + b) según lo observado en el gráfico
    print(f"\n2. Ajuste lineal (ventana óptima: {mejor_ventana:.2f}h)...")
    modelo = calibrador.ajuste_lineal(mejor_ventana)
    print(f"Ecuación: referencia = {modelo['pendiente']:.4f}·IoT + {modelo['intercepto']:.4f}")
    print(f"R²: {modelo['r_cuadrado']:.4f}")
    print(f"RMSE: {modelo['rmse']:.4f}")
    
    # Alcance de validez
    print("\n3. Alcance de validez...")
    tolerancia = 0.25  # 25% de error máximo
    intervalo = calibrador.determinar_intervalo_validez(tolerancia, mejor_ventana)
    if intervalo is not None:
        x_min, x_max, err_prom, err_max = intervalo
        print(f"Tolerancia: {tolerancia:.1%}")
        print(f"Intervalo válido en IoT: [{x_min:.2f}, {x_max:.2f}] μg/m³")
        print(f"Error promedio (válidos): {err_prom:.4f}")
        print(f"Error máximo (válidos): {err_max:.4f}")
    else:
        print("No se encontró intervalo válido con la tolerancia dada.")
    
    # Visualizaciones
    print("\n4. Generando visualizaciones...")
    calibrador.visualizar_resultados(mejor_ventana)
    
    # Tamaño mínimo de datos
    print("\n5. Tamaño mínimo de datos...")
    tam = calibrador.buscar_tamano_minimo_datos(tolerancia, mejor_ventana, paso=200)
    if tam:
        print(f"Tamaño mínimo de entrenamiento: {tam['n_train']} pares")
        iv = tam['intervalo_valido_val']
        print(f"Intervalo válido (validación) IoT: [{iv[0]:.2f}, {iv[1]:.2f}] μg/m³")
    else:
        print("No se encontró tamaño mínimo con la tolerancia dada.")

    # Reporte
    print("\n6. Generando reporte...")
    reporte = calibrador.generar_reporte(mejor_ventana)
    
    # Guardar reporte: imprimir para que quede en el archivo mediante el tee
    print("\n===== REPORTE COMPLETO =====")
    print(reporte)
    print("===== FIN REPORTE =====\n")
    print("Reporte guardado como 'reporte_datos_reales.txt'")
    
    # Guardar datos procesados
    df_sensor.to_csv('Asignaciones/TallerDistancias/datos_sensor_procesados.csv', index=False)
    np.savetxt('Asignaciones/TallerDistancias/datos_referencia_procesados.csv', 
               datos_referencia, delimiter=',', header='tiempo_horas,pm25', comments='')
    print("Datos procesados guardados")
    
    print("\n" + "="*70)
    print("¡ANÁLISIS COMPLETADO!")
    print("="*70)
    print("Archivos generados:")
    print("- resultados_datos_reales.png")
    print("- reporte_datos_reales.txt")
    print("- datos_sensor_procesados.csv")
    print("- datos_referencia_procesados.csv")

if __name__ == "__main__":
    main()
