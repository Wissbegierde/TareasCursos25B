import sympy
from sympy import *

# Vectores a considerar (como ImmutableMatrix)
a = Matrix([1, 2, 3])
b = Matrix([4, 5, 6])
c = Matrix([3, 2, 1])
d = Matrix([6, 5, 4])

# Vectores base
e1 = Matrix([1, 0, 0])
e2 = Matrix([0, 1, 0])
e3 = Matrix([0, 0, 1])

e4 = MatrixBase

# ...resto del código igual...
print("a)----------------------")
sum_1 = a + b + c + d
print("Suma de a+b+c+d:", sum_1)
sum_2 = a + b - c - d
print("Suma de a+b-c-d:", sum_2)
sum_3 = a - b + c - d
print("Suma de a-b+c-d:", sum_3)
sum_4 = -a + b - c + d
print("Suma de -a+b-c+d:", sum_4)

def angle_between(v1, v2):
    dot_product = v1.dot(v2)
    magnitude_v1 = v1.norm()
    magnitude_v2 = v2.norm()
    cos_theta = dot_product / (magnitude_v1 * magnitude_v2)
    angle_rad = acos(cos_theta)
    angle_deg = angle_rad * 180 / pi
    return angle_rad.evalf(), angle_deg.evalf()  # Devuelve el valor numérico en radianes y grados

print("b)-------------------------RADIANES---------GRADOS")
print("Angulo entre a y e1:", angle_between(a, e1))
print("Angulo entre a y e2:", angle_between(a, e2))
print("Angulo entre a y e3:", angle_between(a, e3))
print("Angulo entre b y e1:", angle_between(b, e1))
print("Angulo entre b y e2:", angle_between(b, e2))
print("Angulo entre b y e3:", angle_between(b, e3))
print("Angulo entre c y e1:", angle_between(c, e1))
print("Angulo entre c y e2:", angle_between(c, e2))
print("Angulo entre c y e3:", angle_between(c, e3))
print("Angulo entre d y e1:", angle_between(d, e1))
print("Angulo entre d y e2:", angle_between(d, e2))
print("Angulo entre d y e3:", angle_between(d, e3))

#c) magnitud de los vectores
print("c)----------------------")
print("Magnitud de a:", a.norm().evalf())
print("Magnitud de b:", b.norm().evalf())
print("Magnitud de c:", c.norm().evalf())
print("Magnitud de d:", d.norm().evalf())

# d) El angulo entre a y b y entre c y d
print("d)-------------------------RADIANES---------GRADOS")
print("Angulo entre a y b:", angle_between(a, b))
print("Angulo entre c y d:", angle_between(c, d))

#e) proyeccion de a sobre b 
print("e)----------------------")
def projection(v, u):
    return (v.dot(u) / u.norm()**2) * u
print("Proyección de a sobre b:", projection(a, b))

# f) si los vectores a,b,c,d coplanares
print("f)----------------------")
M = Matrix.hstack(a, b, c, d)
print("Matriz de vectores:")
print(M)
# Verificar rango
if M.rank() <= 2:
    print("Son coplanares")
else:
    print("No son coplanares")

# g) Encuentre (a+b) * (c+d)
print("g)----------------------")
print("Producto punto de a y b:", (a+b).dot(c+d))

# h) Encuentre el producto cruzado de a y b, b y c, c y d, y los angulos que forman con d
print("h)----------------------")
print("Producto cruzado de a y b: ", a.cross(b))
print("Angulo con respecto a d: ", angle_between(a.cross(b), d))
print("Producto cruzado de b y c: ", b.cross(c))
print("Angulo con respecto a d: ", angle_between(b.cross(c), d))
print("Producto cruzado de c y d: ", c.cross(d))
print("Angulo con respecto a d: ", angle_between(c.cross(d), d))
# i) c*(a x b)
print("i)----------------------")
print("Producto punto de c con el producto cruzado de a y b:", c.dot(a.cross(b)))