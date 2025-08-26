from typing import List, Union
from sympy import symbols, Matrix, simplify, integrate, sqrt

# parámetros globales
t = symbols('t')
a, b = -1, 1

def _make_basis(n: int) -> List[Matrix]:
    """Construye la base canónica [1, t, ..., t^(n-1)] como vectores columna."""
    return [Matrix([t**i]) for i in range(n)]

def product(f: Matrix, g: Matrix):
    """
    Producto interno con peso sqrt(1-t^2):
    <f,g> = ∫_{a}^{b} f(t)*g(t)*√(1-t^2) dt
    """
    return simplify(integrate(f[0] * g[0] * sqrt(1 - t**2), (t, a, b)))

def gram_schmidt(basis_or_n: Union[List[Matrix], int]) -> List[Matrix]:
    """
    Gram–Schmidt (sin normalizar) usando el producto definido en este módulo.
    Acepta un entero n (construye la base) o una lista de vectores.
    """
    if isinstance(basis_or_n, int):
        basis = _make_basis(basis_or_n)
    else:
        basis = basis_or_n

    ort: List[Matrix] = []
    for v in basis:
        u = v
        for o in ort:
            coeff = product(v, o) / product(o, o)
            u = simplify(u - coeff * o)
        ort.append(simplify(u))
    return ort

def orthonormalize(basis_or_n: Union[List[Matrix], int]) -> List[Matrix]:
    """
    Normaliza la salida de gram_schmidt (norma 1 con el mismo producto).
    Acepta entero n o lista de vectores.
    """
    if isinstance(basis_or_n, int):
        basis = _make_basis(basis_or_n)
    else:
        basis = basis_or_n

    ort = gram_schmidt(basis)
    normed: List[Matrix] = []
    for u in ort:
        norm_sq = simplify(product(u, u))
        if norm_sq == 0:
            normed.append(u)
        else:
            normed.append(simplify(u / sqrt(norm_sq)))
    return normed

def is_ort_base(vectores_or_n: Union[List[Matrix], int], numeric: bool = False, tol: float = 1e-12) -> bool:
    """
    Comprueba ortogonalidad usando el producto de este módulo.
    Acepta entero n (construye base) o lista de vectores.
    """
    if isinstance(vectores_or_n, int):
        vectores = _make_basis(vectores_or_n)
    else:
        vectores = vectores_or_n

    m = len(vectores)
    for i in range(m):
        for j in range(i+1, m):
            ip = product(vectores[i], vectores[j])
            if numeric:
                if abs(float(ip.evalf())) > tol:
                    return False
            else:
                if simplify(ip) != 0:
                    return False
    return True

if __name__ == "__main__":
    # demo rápido — cambia n aquí si quieres
    n = 5
    base = _make_basis(n)
    print("Base original (grado < n):")
    for v in base:
        print(v.T)
    print("¿Es ortogonal? ->", is_ort_base(base, numeric=True))

    base_ort = gram_schmidt(n)           # construye internamente la base y ortogonaliza
    print("\nBase ortogonalizada (peso sqrt(1-t^2)):")
    for v in base_ort:
        print(v.T)

    base_ortonorm = orthonormalize(n)
    print("\nBase ortonormalizada:")
    for v in base_ortonorm:
        print(v.T)

    print("\nComprobación (ortonorm):", is_ort_base(base_ortonorm, numeric=True))
