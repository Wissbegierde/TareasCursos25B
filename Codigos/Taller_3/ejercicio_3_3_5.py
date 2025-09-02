"""
Taller 3: Productos Tensoriales de Espacios Vectoriales de Polinomios
====================================================================

Este código resuelve el problema de construcción de espacios tensoriales
a partir de espacios vectoriales de polinomios de grado ≤ 2.
"""

import numpy as np
import sympy as sp
from sympy import symbols, expand, simplify, collect
import matplotlib.pyplot as plt

# Configuración de símbolos
x, y = symbols('x y')

class PolinomioLegendre:
    """
    Clase para manejar polinomios de Legendre de grado ≤ 2
    """
    
    def __init__(self):
        # Polinomios de Legendre de grado 0, 1, 2
        self.P0 = 1
        self.P1 = x
        self.P2 = (3*x**2 - 1) / 2
        
        # Base de polinomios de Legendre
        self.base_legendre = [self.P0, self.P1, self.P2]
        
        # Base de monomios
        self.base_monomios = [1, x, x**2] #espacio de polinomios P2(x), esto se puede usar para G2(y)

    
    def expresar_en_legendre(self, polinomio):
        """
        Expresa un polinomio en términos de la base de Legendre
        
        Args:
            polinomio: Polinomio a expresar
            
        Returns:
            dict: Coeficientes en la base de Legendre
        """
        # Expandir el polinomio
        pol_expandido = expand(polinomio)
        
        # Coeficientes en la base de monomios [1, x, x²] - convertir a float
        coef_monomios = [float(pol_expandido.coeff(x, i)) for i in range(3)]
        
        # Cálculo directo de los coeficientes en base de Legendre
        # P0 = 1, P1 = x, P2 = (3x²-1)/2
        # Entonces: x² = (2P2 + 1)/3 = (2P2 + P0)/3
        
        # Para un polinomio a₀ + a₁x + a₂x²:
        # a₀ + a₁x + a₂x² = a₀ + a₁x + a₂(2P2 + P0)/3
        #                  = a₀ + a₁x + (2a₂/3)P2 + (a₂/3)P0
        #                  = (a₀ + a₂/3)P0 + a₁P1 + (2a₂/3)P2
        
        a0, a1, a2 = coef_monomios
        
        coef_P0 = a0 + a2/3
        coef_P1 = a1
        coef_P2 = 2*a2/3
        
        return {
            'P0': coef_P0,
            'P1': coef_P1,
            'P2': coef_P2
        }
    
    def expresar_en_monomios(self, polinomio):
        """
        Expresa un polinomio en términos de la base de monomios
        
        Args:
            polinomio: Polinomio a expresar
            
        Returns:
            list: Coeficientes en la base de monomios [1, x, x²]
        """
        pol_expandido = expand(polinomio)
        return [float(pol_expandido.coeff(x, i)) for i in range(3)]

def producto_tensorial(p1, p2): #por ejemplo le pasamos [1,x,x^2] y [1,y,y^2]
    """
    Calcula el producto tensorial de dos polinomios
    
    Args:
        p1: Primer polinomio (en x)
        p2: Segundo polinomio (en y)
        
    Returns:
        sympy expression: Producto tensorial p1(x) ⊗ p2(y)
    """
    return expand(p1 * p2)

def componentes_tensor_monomios(tensor, grado_max=2):
    """
    Encuentra las componentes del tensor en la base de monomios
    
    Args:
        tensor: Producto tensorial p(x,y)
        grado_max: Grado máximo de los polinomios
        
    Returns:
        numpy array: Matriz de componentes c^ij
    """
    # Expandir el tensor
    tensor_expandido = expand(tensor)
    
    # Inicializar matriz de componentes
    componentes = np.zeros((grado_max + 1, grado_max + 1))
    
    # Extraer coeficientes para cada término x^i * y^j
    for i in range(grado_max + 1):
        for j in range(grado_max + 1):
            # Coeficiente del término x^i * y^j
            coef = tensor_expandido.coeff(x, i).coeff(y, j)
            componentes[i, j] = float(coef)
    
    return componentes

def componentes_tensor_legendre(tensor, grado_max=2):
    """
    Encuentra las componentes del tensor en la base de Legendre
    
    Args:
        tensor: Producto tensorial p(x,y)
        grado_max: Grado máximo de los polinomios
        
    Returns:
        numpy array: Matriz de componentes c̃^ij
    """
    # Expandir el tensor
    tensor_expandido = expand(tensor)
    
    # Inicializar matriz de componentes
    componentes = np.zeros((grado_max + 1, grado_max + 1))
    
    # Para cada término x^i * y^j, convertirlo a base de Legendre
    for i in range(grado_max + 1):
        for j in range(grado_max + 1):
            coef_xy = tensor_expandido.coeff(x, i).coeff(y, j)
            if coef_xy != 0:
                # Convertir x^i a base de Legendre
                coef_x_legendre = convertir_monomio_a_legendre(i)
                # Convertir y^j a base de Legendre  
                coef_y_legendre = convertir_monomio_a_legendre(j)
                
                # El coeficiente del tensor es el producto de los coeficientes
                for k, coef_k in enumerate(coef_x_legendre):
                    for l, coef_l in enumerate(coef_y_legendre):
                        componentes[k, l] += float(coef_xy) * coef_k * coef_l
    
    return componentes

def convertir_monomio_a_legendre(grado):
    """
    Convierte x^grado a base de Legendre
    
    Args:
        grado: Grado del monomio (0, 1, o 2)
        
    Returns:
        list: Coeficientes en base de Legendre [P0, P1, P2]
    """
    if grado == 0:  # x^0 = 1
        return [1, 0, 0]
    elif grado == 1:  # x^1 = x
        return [0, 1, 0]
    elif grado == 2:  # x^2 = (2P2 + 1)/3
        return [1/3, 0, 2/3]
    else:
        raise ValueError("Solo se soportan grados 0, 1, 2")

def main():
    """
    Función principal que resuelve el problema paso a paso
    """
    print("=" * 80)
    print("TALLER 3: PRODUCTOS TENSORIALES DE ESPACIOS VECTORIALES")
    print("=" * 80)
    
    # Inicializar la clase de polinomios de Legendre
    legendre = PolinomioLegendre() #crea una variable de la clase PolinomioLegendre
    
    # PARTE (a): Expresar p^P(x) = x² + x + 3 en base de Legendre
    print("\n" + "="*50)
    print("PARTE (a): Expresión en base de Legendre")
    print("="*50)
    
    p_P = x**2 + x + 3
    print(f"Polinomio original: p^P(x) = {p_P}")
    
    coef_legendre = legendre.expresar_en_legendre(p_P)
    print(f"\nCoeficientes en base de Legendre:")
    print(f"P0 (grado 0): {coef_legendre['P0']:.6f}")
    print(f"P1 (grado 1): {coef_legendre['P1']:.6f}")
    print(f"P2 (grado 2): {coef_legendre['P2']:.6f}")
    
    # Verificación manual
    print(f"\nVerificación manual:")
    print(f"P0 = {legendre.P0}")
    print(f"P1 = {legendre.P1}")
    print(f"P2 = {legendre.P2}")
    
    p_reconstruido = (coef_legendre['P0'] * legendre.P0 + 
                     coef_legendre['P1'] * legendre.P1 + 
                     coef_legendre['P2'] * legendre.P2)
    print(f"Reconstruido: {expand(p_reconstruido)}")
    print(f"Original: {p_P}")
    print(f"¿Coinciden? {simplify(expand(p_reconstruido) - p_P) == 0}")
    print(f"Diferencia: {simplify(expand(p_reconstruido) - p_P)}")
    
    # PARTE (b): Construir el tensor p^P⊗G(x,y)
    print("\n" + "="*50)
    print("PARTE (b): Construcción del tensor")
    print("="*50)
    
    p_G = y + 1
    print(f"p^P(x) = {p_P}")
    print(f"p^G(y) = {p_G}")
    
    tensor = producto_tensorial(p_P, p_G)
    print(f"\nTensor p^P⊗G(x,y) = p^P(x) ⊗ p^G(y) = {tensor}")
    
    # PARTE (c): Componentes en base de monomios
    print("\n" + "="*50)
    print("PARTE (c): Componentes en base de monomios")
    print("="*50)
    
    componentes_monomios = componentes_tensor_monomios(tensor)
    print("Matriz de componentes c^ij en base de monomios:")
    print("Base: {1, x, x²} ⊗ {1, y, y²}")
    print("\n     |  1  |  y  | y²")
    print("-----|-----|-----|-----")
    for i, fila in enumerate(componentes_monomios):
        base_x = ["1", "x", "x²"][i]
        print(f"  {base_x}  |", end="")
        for j, coef in enumerate(fila):
            print(f" {coef:4.1f} |", end="")
        print()
    
    # PARTE (d): Componentes en base de Legendre
    print("\n" + "="*50)
    print("PARTE (d): Componentes en base de Legendre")
    print("="*50)
    
    componentes_legendre = componentes_tensor_legendre(tensor)
    print("Matriz de componentes c̃^ij en base de Legendre:")
    print("Base: {P0, P1, P2} ⊗ {P0, P1, P2}")
    print("\n     |  P0  |  P1  |  P2")
    print("-----|------|------|------")
    for i, fila in enumerate(componentes_legendre):
        base_x = ["P0", "P1", "P2"][i]
        print(f"  {base_x}  |", end="")
        for j, coef in enumerate(fila):
            print(f" {coef:5.3f} |", end="")
        print()
    
    # Visualización de las matrices
    print("\n" + "="*50)
    print("VISUALIZACIÓN DE LAS MATRICES")
    print("="*50)
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # Matriz de monomios
    im1 = ax1.imshow(componentes_monomios, cmap='viridis', aspect='equal')
    ax1.set_title('Componentes en Base de Monomios\nc^ij')
    ax1.set_xlabel('Base y: {1, y, y²}')
    ax1.set_ylabel('Base x: {1, x, x²}')
    
    # Añadir valores en la matriz
    for i in range(3):
        for j in range(3):
            ax1.text(j, i, f'{componentes_monomios[i,j]:.1f}', 
                    ha='center', va='center', color='white', fontweight='bold')
    
    # Matriz de Legendre
    im2 = ax2.imshow(componentes_legendre, cmap='viridis', aspect='equal')
    ax2.set_title('Componentes en Base de Legendre\nc̃^ij')
    ax2.set_xlabel('Base y: {P0, P1, P2}')
    ax2.set_ylabel('Base x: {P0, P1, P2}')
    
    # Añadir valores en la matriz
    for i in range(3):
        for j in range(3):
            ax2.text(j, i, f'{componentes_legendre[i,j]:.3f}', 
                    ha='center', va='center', color='white', fontweight='bold')
    
    plt.tight_layout()
    plt.savefig('tensor_components.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    print("\n" + "="*50)
    print("RESUMEN DE RESULTADOS")
    print("="*50)
    print(f"1. p^P(x) = x² + x + 3 en base de Legendre:")
    print(f"   = {coef_legendre['P0']:.3f}P0 + {coef_legendre['P1']:.3f}P1 + {coef_legendre['P2']:.3f}P2")
    print(f"2. Tensor: p^P⊗G(x,y) = {tensor}")
    print(f"3. Componentes en monomios: matriz 3×3 con elementos c^ij")
    print(f"4. Componentes en Legendre: matriz 3×3 con elementos c̃^ij")
    print(f"5. Gráfico guardado como 'tensor_components.png'")

if __name__ == "__main__":
    main()
