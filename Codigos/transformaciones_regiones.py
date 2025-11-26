"""
Script para graficar las regiones del ejercicio 11 y sus imágenes bajo transformaciones:
a) |z| ≤ 1/2, -π/8 < Arg(z) < π/8 bajo w = z²
b) 1 < |z| < 3, 0 < Arg(z) < π/2 bajo w = z³
c) 2 ≤ Im(z) ≤ 5 bajo w = z³
d) |z - 1/2| ≤ 1/2 bajo w = 1/z
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Wedge, Rectangle, Polygon
import matplotlib.patches as mpatches

# ============================================================================
# Ejercicio 11a: |z| ≤ 1/2, -π/8 < Arg(z) < π/8 bajo w = z²
# ============================================================================
fig1, (ax1a, ax1b) = plt.subplots(1, 2, figsize=(14, 6))

# Plano z: Cuña
theta = np.linspace(-np.pi/8, np.pi/8, 100)
r_max = 0.5
# Frontera del círculo
circle_z = Circle((0, 0), r_max, fill=False, edgecolor='b', linewidth=2)
ax1a.add_patch(circle_z)
# Sector angular
wedge = Wedge((0, 0), r_max, -22.5, 22.5, fill=True, alpha=0.3, color='blue', edgecolor='b', linewidth=2)
ax1a.add_patch(wedge)
# Rayos del sector
ax1a.plot([0, r_max*np.cos(-np.pi/8)], [0, r_max*np.sin(-np.pi/8)], 'b-', linewidth=2)
ax1a.plot([0, r_max*np.cos(np.pi/8)], [0, r_max*np.sin(np.pi/8)], 'b-', linewidth=2)

ax1a.set_xlabel('$x$ (parte real)', fontsize=12)
ax1a.set_ylabel('$y$ (parte imaginaria)', fontsize=12)
ax1a.set_title('Plano $z$: $|z| \\leq 1/2$, $-\\pi/8 < \\text{Arg}(z) < \\pi/8$', fontsize=12, fontweight='bold')
ax1a.grid(True, alpha=0.3)
ax1a.set_aspect('equal')
ax1a.set_xlim(-0.6, 0.6)
ax1a.set_ylim(-0.6, 0.6)
ax1a.axhline(y=0, color='k', linewidth=0.5, alpha=0.3)
ax1a.axvline(x=0, color='k', linewidth=0.5, alpha=0.3)

# Plano w: Imagen bajo w = z² (radio máximo 1/4, ángulo -π/4 a π/4)
r_max_w = 0.25
theta_w = np.linspace(-np.pi/4, np.pi/4, 100)
circle_w = Circle((0, 0), r_max_w, fill=False, edgecolor='r', linewidth=2)
ax1b.add_patch(circle_w)
wedge_w = Wedge((0, 0), r_max_w, -45, 45, fill=True, alpha=0.3, color='red', edgecolor='r', linewidth=2)
ax1b.add_patch(wedge_w)
ax1b.plot([0, r_max_w*np.cos(-np.pi/4)], [0, r_max_w*np.sin(-np.pi/4)], 'r-', linewidth=2)
ax1b.plot([0, r_max_w*np.cos(np.pi/4)], [0, r_max_w*np.sin(np.pi/4)], 'r-', linewidth=2)

ax1b.set_xlabel('$u$ (parte real)', fontsize=12)
ax1b.set_ylabel('$v$ (parte imaginaria)', fontsize=12)
ax1b.set_title('Plano $w$: Imagen bajo $w = z^2$', fontsize=12, fontweight='bold')
ax1b.grid(True, alpha=0.3)
ax1b.set_aspect('equal')
ax1b.set_xlim(-0.3, 0.3)
ax1b.set_ylim(-0.3, 0.3)
ax1b.axhline(y=0, color='k', linewidth=0.5, alpha=0.3)
ax1b.axvline(x=0, color='k', linewidth=0.5, alpha=0.3)

plt.tight_layout()
plt.savefig('Codigos/ejercicio11a.png', dpi=300, bbox_inches='tight')
print("Gráfica guardada como: Codigos/ejercicio11a.png")
plt.close()

# ============================================================================
# Ejercicio 11b: 1 < |z| < 3, 0 < Arg(z) < π/2 bajo w = z³
# ============================================================================
fig2, (ax2a, ax2b) = plt.subplots(1, 2, figsize=(14, 6))

# Plano z: Sector de anillo
# Círculos exteriores e interiores
circle_z1 = Circle((0, 0), 1, fill=False, edgecolor='b', linewidth=2, linestyle='--')
circle_z2 = Circle((0, 0), 3, fill=False, edgecolor='b', linewidth=2)
ax2a.add_patch(circle_z1)
ax2a.add_patch(circle_z2)
# Sector angular (primer cuadrante)
theta_b = np.linspace(0, np.pi/2, 100)
# Rayos
ax2a.plot([0, 3*np.cos(0)], [0, 3*np.sin(0)], 'b-', linewidth=2)
ax2a.plot([0, 3*np.cos(np.pi/2)], [0, 3*np.sin(np.pi/2)], 'b-', linewidth=2)
# Relleno del sector
r_vals = np.linspace(1, 3, 50)
for r in [1.5, 2, 2.5]:
    x_arc = r * np.cos(theta_b)
    y_arc = r * np.sin(theta_b)
    ax2a.fill_between(x_arc, y_arc, alpha=0.2, color='blue')

ax2a.set_xlabel('$x$ (parte real)', fontsize=12)
ax2a.set_ylabel('$y$ (parte imaginaria)', fontsize=12)
ax2a.set_title('Plano $z$: $1 < |z| < 3$, $0 < \\text{Arg}(z) < \\pi/2$', fontsize=12, fontweight='bold')
ax2a.grid(True, alpha=0.3)
ax2a.set_aspect('equal')
ax2a.set_xlim(-0.5, 3.5)
ax2a.set_ylim(-0.5, 3.5)
ax2a.axhline(y=0, color='k', linewidth=0.5, alpha=0.3)
ax2a.axvline(x=0, color='k', linewidth=0.5, alpha=0.3)

# Plano w: Imagen bajo w = z³ (1 < |w| < 27, 0 < Arg(w) < 3π/2)
circle_w1 = Circle((0, 0), 1, fill=False, edgecolor='r', linewidth=2, linestyle='--')
circle_w2 = Circle((0, 0), 27, fill=False, edgecolor='r', linewidth=2)
ax2b.add_patch(circle_w1)
ax2b.add_patch(circle_w2)
# Sector angular (0 a 3π/2)
theta_w_b = np.linspace(0, 3*np.pi/2, 200)
ax2b.plot([0, 27*np.cos(0)], [0, 27*np.sin(0)], 'r-', linewidth=2)
ax2b.plot([0, 27*np.cos(3*np.pi/2)], [0, 27*np.sin(3*np.pi/2)], 'r-', linewidth=2)
ax2b.plot([0, 27*np.cos(np.pi/2)], [0, 27*np.sin(np.pi/2)], 'r-', linewidth=2)

ax2b.set_xlabel('$u$ (parte real)', fontsize=12)
ax2b.set_ylabel('$v$ (parte imaginaria)', fontsize=12)
ax2b.set_title('Plano $w$: Imagen bajo $w = z^3$', fontsize=12, fontweight='bold')
ax2b.grid(True, alpha=0.3)
ax2b.set_aspect('equal')
ax2b.set_xlim(-30, 30)
ax2b.set_ylim(-30, 30)
ax2b.axhline(y=0, color='k', linewidth=0.5, alpha=0.3)
ax2b.axvline(x=0, color='k', linewidth=0.5, alpha=0.3)

plt.tight_layout()
plt.savefig('Codigos/ejercicio11b.png', dpi=300, bbox_inches='tight')
print("Gráfica guardada como: Codigos/ejercicio11b.png")
plt.close()

# ============================================================================
# Ejercicio 11c: 2 ≤ Im(z) ≤ 5 bajo w = z³
# ============================================================================
fig3, (ax3a, ax3b) = plt.subplots(1, 2, figsize=(14, 6))

# Plano z: Banda horizontal
ax3a.axhline(y=2, color='b', linewidth=2, label='$y = 2$')
ax3a.axhline(y=5, color='b', linewidth=2, label='$y = 5$')
ax3a.fill_between([-3, 3], 2, 5, alpha=0.3, color='blue')

ax3a.set_xlabel('$x$ (parte real)', fontsize=12)
ax3a.set_ylabel('$y$ (parte imaginaria)', fontsize=12)
ax3a.set_title('Plano $z$: $2 \\leq \\text{Im}(z) \\leq 5$', fontsize=12, fontweight='bold')
ax3a.grid(True, alpha=0.3)
ax3a.set_aspect('equal')
ax3a.set_xlim(-3, 3)
ax3a.set_ylim(0, 6)
ax3a.axhline(y=0, color='k', linewidth=0.5, alpha=0.3)
ax3a.axvline(x=0, color='k', linewidth=0.5, alpha=0.3)
ax3a.legend(fontsize=10)

# Plano w: Imagen bajo w = z³ (curvas cúbicas)
x_vals = np.linspace(-3, 3, 500)
# Para y = 2
u_2 = x_vals**3 - 3*x_vals*4
v_2 = 3*x_vals**2*2 - 8
# Para y = 5
u_5 = x_vals**3 - 3*x_vals*25
v_5 = 3*x_vals**2*5 - 125

ax3b.plot(u_2, v_2, 'r-', linewidth=2, label='Imagen de $y = 2$')
ax3b.plot(u_5, v_5, 'r-', linewidth=2, label='Imagen de $y = 5$')
# Relleno aproximado
ax3b.fill_between(u_2, v_2, v_5, where=(u_2 <= u_5), alpha=0.3, color='red', interpolate=True)

ax3b.set_xlabel('$u$ (parte real)', fontsize=12)
ax3b.set_ylabel('$v$ (parte imaginaria)', fontsize=12)
ax3b.set_title('Plano $w$: Imagen bajo $w = z^3$', fontsize=12, fontweight='bold')
ax3b.grid(True, alpha=0.3)
ax3b.set_aspect('equal')
ax3b.legend(fontsize=10)
ax3b.axhline(y=0, color='k', linewidth=0.5, alpha=0.3)
ax3b.axvline(x=0, color='k', linewidth=0.5, alpha=0.3)

plt.tight_layout()
plt.savefig('Codigos/ejercicio11c.png', dpi=300, bbox_inches='tight')
print("Gráfica guardada como: Codigos/ejercicio11c.png")
plt.close()

# ============================================================================
# Ejercicio 11d: |z - 1/2| ≤ 1/2 bajo w = 1/z
# ============================================================================
fig4, (ax4a, ax4b) = plt.subplots(1, 2, figsize=(14, 6))

# Plano z: Disco centrado en 1/2
circle_z_d = Circle((0.5, 0), 0.5, fill=True, alpha=0.3, color='blue', edgecolor='b', linewidth=2)
ax4a.add_patch(circle_z_d)

ax4a.set_xlabel('$x$ (parte real)', fontsize=12)
ax4a.set_ylabel('$y$ (parte imaginaria)', fontsize=12)
ax4a.set_title('Plano $z$: $|z - 1/2| \\leq 1/2$', fontsize=12, fontweight='bold')
ax4a.grid(True, alpha=0.3)
ax4a.set_aspect('equal')
ax4a.set_xlim(-0.5, 1.5)
ax4a.set_ylim(-0.8, 0.8)
ax4a.axhline(y=0, color='k', linewidth=0.5, alpha=0.3)
ax4a.axvline(x=0, color='k', linewidth=0.5, alpha=0.3)
ax4a.plot(0.5, 0, 'ko', markersize=5)
ax4a.annotate('Centro: $1/2$', xy=(0.5, 0), xytext=(0.6, 0.3), fontsize=10)

# Plano w: Imagen bajo w = 1/z (semiplano u ≥ 1)
ax4b.axvline(x=1, color='r', linewidth=2, label='$u = 1$')
ax4b.fill_between([1, 2], -1, 1, alpha=0.3, color='red')

ax4b.set_xlabel('$u$ (parte real)', fontsize=12)
ax4b.set_ylabel('$v$ (parte imaginaria)', fontsize=12)
ax4b.set_title('Plano $w$: Imagen bajo $w = 1/z$ (semiplano $u \\geq 1$)', fontsize=12, fontweight='bold')
ax4b.grid(True, alpha=0.3)
ax4b.set_aspect('equal')
ax4b.set_xlim(-0.5, 2)
ax4b.set_ylim(-1, 1)
ax4b.axhline(y=0, color='k', linewidth=0.5, alpha=0.3)
ax4b.axvline(x=0, color='k', linewidth=0.5, alpha=0.3)
ax4b.legend(fontsize=10)

plt.tight_layout()
plt.savefig('Codigos/ejercicio11d.png', dpi=300, bbox_inches='tight')
print("Gráfica guardada como: Codigos/ejercicio11d.png")
plt.close()

print("\n" + "="*70)
print("Todas las gráficas del ejercicio 11 han sido generadas exitosamente.")
print("="*70)

