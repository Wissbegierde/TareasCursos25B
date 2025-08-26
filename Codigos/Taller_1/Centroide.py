import sympy 
from sympy import * 

# Ejemplo de tres vértices para un triángulo, cambiar a gusto
x = Matrix([0, 0])
y = Matrix([4, 0])
z = Matrix([2, 3])


def is_triangle(a, b, c):
    det = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
    return det != 0

def centroide(a, b, c):
    if not is_triangle(a, b, c):
        raise ValueError("Los puntos no forman un triángulo válido.")
    cx = (a[0] + b[0] + c[0]) / 3
    cy = (a[1] + b[1] + c[1]) / 3
    return Point2D(cx, cy)

print("Centroide del triángulo formado por los puntos a, b, c:")
print(centroide(x, y, z))
