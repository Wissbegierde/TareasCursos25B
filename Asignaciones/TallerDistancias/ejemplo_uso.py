"""
Ejemplo de uso del Taller de Distancias
Este script muestra cómo usar el calibrador con datos reales o simulados.
"""

import numpy as np
import matplotlib.pyplot as plt
import os
from taller_distancias import CalibradorSensores, generar_datos_ejemplo

def ejemplo_datos_simulados():
    """
    Ejemplo usando datos simulados para demostrar las capacidades del calibrador.
    """
    print("="*60)
    print("EJEMPLO 1: DATOS SIMULADOS")
    print("="*60)
    
    # Generar datos de ejemplo
    datos_ref, datos_sensor = generar_datos_ejemplo(n_puntos=300, ruido=0.2)
    
    # Crear calibrador
    calibrador = CalibradorSensores(datos_ref, datos_sensor)
    
    # Análisis completo
    print("Realizando análisis completo...")
    
    # 1. Encontrar mejor ventana
    ventanas = np.linspace(0.5, 6.0, 25)
    mejor_ventana, menor_distancia, distancias = calibrador.encontrar_mejor_ventana(ventanas)
    
    # 2. Ajuste lineal
    modelo = calibrador.ajuste_lineal(mejor_ventana)
    
    # 3. Alcance de validez
    tolerancia = 0.2  # 20% de error máximo
    alcance = calibrador.determinar_alcance_validez(tolerancia)
    
    # 4. Visualizaciones
    calibrador.visualizar_resultados(mejor_ventana)
    
    # 5. Reporte
    reporte = calibrador.generar_reporte(mejor_ventana)
    print(reporte)
    
    return calibrador

def ejemplo_datos_reales():
    """
    Ejemplo de cómo cargar y usar datos reales.
    """
    print("="*60)
    print("EJEMPLO 2: DATOS REALES")
    print("="*60)
    print("Para usar datos reales, sigue estos pasos:")
    print("1. Descarga los datos del repositorio:")
    print("   https://github.com/nunezluis/MisCursos/tree/main/MisMateriales/Asignaciones/TallerDistancias/DatosDistancias")
    print("2. Carga los archivos CSV:")
    print("   - 'Datos Estaciones AMB' (datos de referencia)")
    print("   - 'mediciones...' (datos del sensor IoT)")
    print("3. Usa el siguiente código:")
    print()
    print("""
# Cargar datos reales
import pandas as pd

# Cargar datos de referencia
df_ref = pd.read_csv('Datos Estaciones AMB.csv')
tiempos_ref = df_ref['tiempo'].values  # Ajustar nombre de columna
concentraciones_ref = df_ref['PM2.5'].values  # Ajustar nombre de columna
datos_ref = np.column_stack([tiempos_ref, concentraciones_ref])

# Cargar datos del sensor
df_sensor = pd.read_csv('mediciones_sensor.csv')
tiempos_sensor = df_sensor['tiempo'].values
concentraciones_sensor = df_sensor['PM2.5'].values
datos_sensor = np.column_stack([tiempos_sensor, concentraciones_sensor])

# Crear calibrador
calibrador = CalibradorSensores(datos_ref, datos_sensor)

# Realizar análisis
ventanas = np.linspace(0.5, 5.0, 20)
mejor_ventana, _, _ = calibrador.encontrar_mejor_ventana(ventanas)
modelo = calibrador.ajuste_lineal(mejor_ventana)
alcance = calibrador.determinar_alcance_validez(0.15)
calibrador.visualizar_resultados(mejor_ventana)
    """)

def analisis_sensibilidad():
    """
    Análisis de sensibilidad del método a diferentes parámetros.
    """
    print("="*60)
    print("ANÁLISIS DE SENSIBILIDAD")
    print("="*60)
    
    # Generar datos base
    datos_ref, datos_sensor = generar_datos_ejemplo(n_puntos=200, ruido=0.1)
    
    # Probar diferentes niveles de ruido
    niveles_ruido = [0.05, 0.1, 0.2, 0.3, 0.5]
    resultados = []
    
    for ruido in niveles_ruido:
        # Generar datos con diferente ruido
        _, datos_sensor_ruido = generar_datos_ejemplo(n_puntos=200, ruido=ruido)
        calibrador = CalibradorSensores(datos_ref, datos_sensor_ruido)
        
        # Encontrar mejor ventana
        ventanas = np.linspace(0.5, 4.0, 15)
        mejor_ventana, menor_distancia, _ = calibrador.encontrar_mejor_ventana(ventanas)
        
        # Ajuste lineal
        modelo = calibrador.ajuste_lineal(mejor_ventana)
        
        resultados.append({
            'ruido': ruido,
            'mejor_ventana': mejor_ventana,
            'distancia': menor_distancia,
            'r_cuadrado': modelo['r_cuadrado'],
            'rmse': modelo['rmse']
        })
    
    # Mostrar resultados
    print("Resultados del análisis de sensibilidad:")
    print("-" * 80)
    print(f"{'Ruido':<8} {'Ventana':<10} {'Distancia':<12} {'R²':<8} {'RMSE':<8}")
    print("-" * 80)
    
    for r in resultados:
        print(f"{r['ruido']:<8.2f} {r['mejor_ventana']:<10.2f} {r['distancia']:<12.4f} "
              f"{r['r_cuadrado']:<8.3f} {r['rmse']:<8.3f}")
    
    # Visualizar sensibilidad
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle('Análisis de Sensibilidad', fontsize=16, fontweight='bold')
    
    ruidos = [r['ruido'] for r in resultados]
    ventanas = [r['mejor_ventana'] for r in resultados]
    distancias = [r['distancia'] for r in resultados]
    r_cuadrados = [r['r_cuadrado'] for r in resultados]
    rmses = [r['rmse'] for r in resultados]
    
    # Ventana óptima vs ruido
    axes[0, 0].plot(ruidos, ventanas, 'bo-', linewidth=2, markersize=8)
    axes[0, 0].set_xlabel('Nivel de Ruido')
    axes[0, 0].set_ylabel('Ventana Óptima (h)')
    axes[0, 0].set_title('Ventana Óptima vs Ruido')
    axes[0, 0].grid(True, alpha=0.3)
    
    # Distancia vs ruido
    axes[0, 1].plot(ruidos, distancias, 'ro-', linewidth=2, markersize=8)
    axes[0, 1].set_xlabel('Nivel de Ruido')
    axes[0, 1].set_ylabel('Distancia Euclidiana')
    axes[0, 1].set_title('Distancia vs Ruido')
    axes[0, 1].grid(True, alpha=0.3)
    
    # R² vs ruido
    axes[1, 0].plot(ruidos, r_cuadrados, 'go-', linewidth=2, markersize=8)
    axes[1, 0].set_xlabel('Nivel de Ruido')
    axes[1, 0].set_ylabel('R²')
    axes[1, 0].set_title('Calidad del Ajuste vs Ruido')
    axes[1, 0].grid(True, alpha=0.3)
    
    # RMSE vs ruido
    axes[1, 1].plot(ruidos, rmses, 'mo-', linewidth=2, markersize=8)
    axes[1, 1].set_xlabel('Nivel de Ruido')
    axes[1, 1].set_ylabel('RMSE')
    axes[1, 1].set_title('Error Cuadrático vs Ruido')
    axes[1, 1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    os.makedirs('Asignaciones/TallerDistancias', exist_ok=True)
    plt.savefig('Asignaciones/TallerDistancias/analisis_sensibilidad.png', 
                dpi=300, bbox_inches='tight')
    plt.show()
    
    print(f"\nGráfico de sensibilidad guardado como 'analisis_sensibilidad.png'")

def main():
    """
    Función principal que ejecuta todos los ejemplos.
    """
    print("TALLER DE DISTANCIAS - EJEMPLOS DE USO")
    print("=" * 50)
    
    # Ejecutar ejemplos
    calibrador = ejemplo_datos_simulados()
    ejemplo_datos_reales()
    analisis_sensibilidad()
    
    print("\n" + "="*60)
    print("¡ANÁLISIS COMPLETADO!")
    print("="*60)
    print("Archivos generados:")
    print("- resultados_calibracion.png")
    print("- analisis_sensibilidad.png")
    print("- reporte_calibracion.txt")
    print("\nPara usar con tus propios datos:")
    print("1. Modifica el código en ejemplo_datos_reales()")
    print("2. Ajusta los nombres de las columnas según tus datos")
    print("3. Ejecuta el análisis completo")

if __name__ == "__main__":
    main()
