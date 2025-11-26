"""
Script para graficar las curvas del ejercicio 7:
- Rectas originales en el plano z: y = x - 1 y y = 0
- Imágenes transformadas en el plano w bajo w = 1/z
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle

# Configuración de la figura
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# ============================================================================
# Plano z: Curvas originales
# ============================================================================
ax1 = axes[0]

# Recta y = x - 1
x1 = np.linspace(-2, 3, 100)
y1 = x1 - 1
ax1.plot(x1, y1, 'b-', linewidth=2, label=r'$y = x - 1$')

# Recta y = 0 (eje x)
x2 = np.linspace(-2, 3, 100)
y2 = np.zeros_like(x2)
ax1.plot(x2, y2, 'r-', linewidth=2, label=r'$y = 0$')

# Marcar el punto z = 1
ax1.plot(1, 0, 'ko', markersize=8, label=r'$z_0 = 1$')
ax1.annotate('$z_0 = 1$', xy=(1, 0), xytext=(1.2, 0.2), fontsize=10)

# Configurar el plano z
ax1.set_xlabel('$x$ (parte real)', fontsize=12)
ax1.set_ylabel('$y$ (parte imaginaria)', fontsize=12)
ax1.set_title('Plano $z$: Curvas originales', fontsize=14, fontweight='bold')
ax1.grid(True, alpha=0.3)
ax1.legend(loc='upper right', fontsize=10)
ax1.set_aspect('equal')
ax1.set_xlim(-2, 3)
ax1.set_ylim(-2, 2)
ax1.axhline(y=0, color='k', linewidth=0.5, alpha=0.3)
ax1.axvline(x=0, color='k', linewidth=0.5, alpha=0.3)

# ============================================================================
# Plano w: Curvas transformadas
# ============================================================================
ax2 = axes[1]

# Círculo: u^2 + v^2 - u - v = 0
# Completando cuadrados: (u - 1/2)^2 + (v - 1/2)^2 = (1/√2)^2
center = (0.5, 0.5)
radius = 1/np.sqrt(2)
circle = Circle(center, radius, fill=False, edgecolor='b', linewidth=2, label=r'$u^2 + v^2 - u - v = 0$')
ax2.add_patch(circle)

# Recta v = 0 (eje u)
u = np.linspace(-1, 2, 100)
v = np.zeros_like(u)
ax2.plot(u, v, 'r-', linewidth=2, label=r'$v = 0$')

# Marcar el punto w = 1 (imagen de z = 1)
ax2.plot(1, 0, 'ko', markersize=8, label=r'$w_0 = 1$')
ax2.annotate('$w_0 = 1$', xy=(1, 0), xytext=(1.2, 0.2), fontsize=10)

# Configurar el plano w
ax2.set_xlabel('$u$ (parte real)', fontsize=12)
ax2.set_ylabel('$v$ (parte imaginaria)', fontsize=12)
ax2.set_title('Plano $w$: Imágenes bajo $w = 1/z$', fontsize=14, fontweight='bold')
ax2.grid(True, alpha=0.3)
ax2.legend(loc='upper right', fontsize=10)
ax2.set_aspect('equal')
ax2.set_xlim(-1, 2)
ax2.set_ylim(-1, 1.5)
ax2.axhline(y=0, color='k', linewidth=0.5, alpha=0.3)
ax2.axvline(x=0, color='k', linewidth=0.5, alpha=0.3)

plt.tight_layout()
plt.savefig('Codigos/transformacion_inversa_ejercicio7.png', dpi=300, bbox_inches='tight')
print("Gráfica guardada como: Codigos/transformacion_inversa_ejercicio7.png")
plt.show()

