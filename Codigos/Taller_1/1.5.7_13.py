# Calcular el trabajo hecho hecho en contra 
# de este campo de fuerza al moverse alrededor de un circulo
# de radio 1 y en el playo x-y
# a) desde 0 a pi, en sentido contrario a la agujas del reloj
# b) desde 0 a -pi en sentido de las agujas del reloj

import sympy
from sympy import *
from sympy.vector import *

# Sistema de coordenadas
R = CoordSys3D('R')

# Campo de fuerza
F = -(R.y / (R.x**2 + R.y**2)) * R.i + (R.x / (R.x**2 + R.y**2)) * R.j

# Parámetro
t = symbols('t', real=True)

# Curva en el plano x-y (radio 1)
r = cos(t)*R.i + sin(t)*R.j + 0*R.k   # Vector posición
dr_dt = diff(r, t)                    # Derivada r'(t)

# Sustitución de x, y, z en el campo
F_sub = F.subs({R.x: cos(t), R.y: sin(t), R.z: 0})

# Producto punto y luego integración en contra del campo de fuerza
W1 = -integrate(F_sub.dot(dr_dt), (t, 0, pi))
W2 = -integrate(F_sub.dot(dr_dt), (t, 0, -pi))

print("a) Trabajo 0→π CCW:", W1)
print("b) Trabajo 0→−π CW:", W2)

