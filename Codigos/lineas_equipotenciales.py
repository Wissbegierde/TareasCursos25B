"""
Programa para graficar líneas equipotenciales u = constante de una función armónica
y de su conjugada v en los mismos ejes.

"""

import numpy as np
import matplotlib.pyplot as plt
from sympy import symbols, cos, sin, exp, lambdify
from matplotlib.colors import ListedColormap

# Definir símbolos
x, y = symbols('x y')

def graficar_lineas_equipotenciales(u_expr, v_expr, x_range=(-3, 3), y_range=(-3, 3), 
                                     n_contours=20, titulo="Líneas Equipotenciales"):
    """
    Grafica las líneas equipotenciales de u y v en los mismos ejes.
    
    Parámetros:
    -----------
    u_expr : expresión SymPy
        Función armónica u(x, y)
    v_expr : expresión SymPy
        Conjugada armónica v(x, y)
    x_range : tuple
        Rango de x para la gráfica (x_min, x_max)
    y_range : tuple
        Rango de y para la gráfica (y_min, y_max)
    n_contours : int
        Número de curvas de nivel
    titulo : str
        Título de la gráfica
    """
    # Convertir expresiones SymPy a funciones numéricas
    u_func = lambdify((x, y), u_expr, 'numpy')
    v_func = lambdify((x, y), v_expr, 'numpy')
    
    # Crear malla de puntos
    x_vals = np.linspace(x_range[0], x_range[1], 500)
    y_vals = np.linspace(y_range[0], y_range[1], 500)
    X, Y = np.meshgrid(x_vals, y_vals)
    
    # Evaluar funciones en la malla
    U = u_func(X, Y)
    V = v_func(X, Y)
    
    # Crear figura
    fig, ax = plt.subplots(figsize=(10, 10))
    
    # Graficar líneas equipotenciales de u (en azul)
    contours_u = ax.contour(X, Y, U, levels=n_contours, colors='blue', 
                            linestyles='solid', linewidths=1.5, alpha=0.7)
    ax.clabel(contours_u, inline=True, fontsize=8, fmt='%1.1f')
    
    # Graficar líneas equipotenciales de v (en rojo)
    contours_v = ax.contour(X, Y, V, levels=n_contours, colors='red', 
                            linestyles='dashed', linewidths=1.5, alpha=0.7)
    ax.clabel(contours_v, inline=True, fontsize=8, fmt='%1.1f')
    
    # Configurar gráfica
    ax.set_xlabel('x', fontsize=12)
    ax.set_ylabel('y', fontsize=12)
    ax.set_title(titulo, fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.set_aspect('equal')
    ax.legend(['u = constante (azul)', 'v = constante (rojo)'], 
              loc='upper right', fontsize=10)
    
    plt.tight_layout()
    return fig, ax

# ============================================================================
# Ejercicio 7: Aplicar el programa a las siguientes funciones
# ============================================================================

if __name__ == "__main__":
    print("=" * 70)
    print("Programa para graficar líneas equipotenciales")
    print("=" * 70)
    
    # a) u = x^2 - y^2, v = 2xy
    print("\na) u = x² - y², v = 2xy")
    u_a = x**2 - y**2
    v_a = 2*x*y
    fig_a, ax_a = graficar_lineas_equipotenciales(
        u_a, v_a, 
        x_range=(-3, 3), 
        y_range=(-3, 3),
        n_contours=15,
        titulo="a) u = x² - y², v = 2xy"
    )
    plt.savefig('Codigos/lineas_equipotenciales_a.png', dpi=300, bbox_inches='tight')
    print("   Gráfica guardada como: lineas_equipotenciales_a.png")
    
    # b) u = x^3 - 3xy^2, v = 3x^2y - y^3
    print("\nb) u = x³ - 3xy², v = 3x²y - y³")
    u_b = x**3 - 3*x*y**2
    v_b = 3*x**2*y - y**3
    fig_b, ax_b = graficar_lineas_equipotenciales(
        u_b, v_b,
        x_range=(-2, 2),
        y_range=(-2, 2),
        n_contours=20,
        titulo="b) u = x³ - 3xy², v = 3x²y - y³"
    )
    plt.savefig('Codigos/lineas_equipotenciales_b.png', dpi=300, bbox_inches='tight')
    print("   Gráfica guardada como: lineas_equipotenciales_b.png")
    
    # c) u = e^x cos(y), v = e^x sin(y)
    print("\nc) u = eˣ cos(y), v = eˣ sin(y)")
    u_c = exp(x) * cos(y)
    v_c = exp(x) * sin(y)
    fig_c, ax_c = graficar_lineas_equipotenciales(
        u_c, v_c,
        x_range=(-2, 2),
        y_range=(-3, 3),
        n_contours=20,
        titulo="c) u = eˣ cos(y), v = eˣ sin(y)"
    )
    plt.savefig('Codigos/lineas_equipotenciales_c.png', dpi=300, bbox_inches='tight')
    print("   Gráfica guardada como: lineas_equipotenciales_c.png")
    
    print("\n" + "=" * 70)
    print("Todas las gráficas han sido generadas exitosamente.")
    print("=" * 70)
    
    # Mostrar todas las gráficas
    plt.show()

