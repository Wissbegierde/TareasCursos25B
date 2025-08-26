from typing import List, Union
from sympy import *
from ejercicio_2_3_6_6 import *
import matplotlib.pyplot as plt
import numpy as np
# parámetros globales
t = symbols('t')
a, b = -1, 1

def func_product(f: Matrix, g: Matrix):
    """Producto interno sin peso: ∫_{a}^{b} f(t) g(t) dt"""
    return simplify(integrate(f[0] * g[0], (t, a, b)))

def _make_basis(n: int) -> List[Matrix]:
    """Construye la base canónica [1, t, ..., t^(n-1)]."""
    return [Matrix([t**i]) for i in range(n)]

def gram_schmidt_custom(basis_or_n: Union[List[Matrix], int], product_func) -> List[Matrix]:
    """
    Gram–Schmidt que usa product_func.
    Si el primer argumento es int, construye la base hasta grado n-1 internamente.
    """
    if isinstance(basis_or_n, int):
        basis = _make_basis(basis_or_n)
    else:
        basis = basis_or_n

    ort: List[Matrix] = []
    for v in basis:
        u = v
        for o in ort:
            coeff = product_func(v, o) / product_func(o, o)
            u = simplify(u - coeff * o)
        ort.append(simplify(u))
    return ort

def orthonormalize_custom(basis_or_n: Union[List[Matrix], int], product_func) -> List[Matrix]:
    """
    Normaliza según product_func. Acepta entero n o lista de vectores.
    """
    if isinstance(basis_or_n, int):
        basis = _make_basis(basis_or_n)
    else:
        basis = basis_or_n

    ort = gram_schmidt_custom(basis, product_func)
    normed: List[Matrix] = []
    for u in ort:
        norm_sq = simplify(product_func(u, u))
        if norm_sq == 0:
            normed.append(u)
        else:
            normed.append(simplify(u / sqrt(norm_sq)))
    return normed

def is_ort_base_custom(vectores_or_n: Union[List[Matrix], int], product_func, numeric: bool = False, tol: float = 1e-12) -> bool:
    """
    Comprueba ortogonalidad usando product_func.
    Acepta entero n (construye la base) o lista de vectores.
    """
    if isinstance(vectores_or_n, int):
        vectores = _make_basis(vectores_or_n)
    else:
        vectores = vectores_or_n

    m = len(vectores)
    for i in range(m):
        for j in range(i+1, m):
            ip = product_func(vectores[i], vectores[j])
            if numeric:
                if abs(float(ip.evalf())) > tol:
                    return False
            else:
                if simplify(ip) != 0:
                    return False
    return True

def expansion_in_legendre(h, n: int):
    """
    Expande h en los primeros n polinomios de Legendre (usando orthonormalize_custom y func_product).
    Devuelve (approx_matrix, coeffs_list).
    """
    P = orthonormalize_custom(n, func_product)  # lista de sympy.Matrix

    # convertir h a Matrix columna si se pasa expresión
    if not isinstance(h, Matrix):
        h = Matrix([h])

    if len(P) == 0:
        return Matrix.zeros(1, 1), []

    # crear approx como Matrix de ceros con la misma forma que los p
    approx = Matrix.zeros(*P[0].shape)
    coeffs = []
    for p in P:
        coeff = simplify(func_product(h, p))   # si P es ortonormal, solo <h,p>
        coeffs.append(coeff)
        approx = approx + coeff * p

    return simplify(approx), coeffs

def _chebyshev_u_polys(n):
    """Devuelve lista de polinomios simbólicos U_0..U_{n-1} (segunda especie)."""
    U = []
    if n == 0:
        return U
    U0 = 1
    U.append(U0)
    if n == 1:
        return U
    U1 = 2*t
    U.append(U1)
    for k in range(2, n):
        Uk = simplify(2*t*U[-1] - U[-2])
        U.append(Uk)
    return U

def expansion_in_chebyshev(h, n: int, Nquad: int = 2001):
    """
    Expande h en los primeros n polinomios de Chebyshev (2ª especie, peso sqrt(1-t^2)).
    - h: expresión sympy o Matrix([expr])
    - n: número de polinomios
    - Nquad: número de puntos para la cuadratura (debe ser impar preferiblemente)
    Devuelve: (approx_matrix_sympy, coeffs_list_sympyfloats)
    """
    # convertir h a expresión si viene en Matrix
    if isinstance(h, Matrix):
        h_expr = h[0]
    else:
        h_expr = h

    # construir polinomios simbólicos U_0..U_{n-1}
    U_sym = _chebyshev_u_polys(n)
    if len(U_sym) == 0:
        return Matrix([0]), []

    # preparar cuadratura en [-1,1]
    x = np.linspace(-1.0, 1.0, Nquad)
    w = np.sqrt(np.clip(1.0 - x**2, 0.0, None))   # peso sqrt(1-x^2), clip evita nan por redondeo


    f_h = lambdify(t, h_expr, 'numpy')
    U_num = [lambdify(t, Ui, 'numpy') for Ui in U_sym]

    # evaluar en la malla
    h_vals = f_h(x)
    coeffs = []
    approx_num_vals = np.zeros_like(x, dtype=float)

    for Ui_sym, Ui_num in zip(U_sym, U_num):
        ui_vals = Ui_num(x)
        # numérico: ∫ h * Ui * weight
        num = np.trapezoid(h_vals * ui_vals * w, x)
        # den = ∫ Ui^2 * weight
        den = np.trapezoid(ui_vals * ui_vals * w, x)
        # evitar división por cero
        if abs(den) < 1e-16:
            ci = 0.0
        else:
            ci = num / den
        coeffs.append(Float(ci))            # guardamos como sympy.Float para impresión
        approx_num_vals += ci * ui_vals

    # construir aproximación simbólica (coeficientes numéricos como Float)
    approx_sym = sum(c * U for c, U in zip(coeffs, U_sym))
    approx_sym = simplify(approx_sym)
    approx_mat = Matrix([approx_sym])

    return approx_mat, coeffs



if __name__ == "__main__":
    
    # a) base {1,t,t^2}
    base3 = _make_basis(3)
    print("a) Base {1, t, t^2}:")
    for v in base3:
        print(v.T)
    print("¿Es ortogonal? ->", is_ort_base_custom(base3, func_product, numeric=True))

    # b) primeros n polinomios ortogonales (Legendre)
    n = 10
    ort_n = gram_schmidt_custom(n, func_product)
    print(f"\nb) Primeros {n} polinomios ortogonales (Gram–Schmidt):")
    for i, v in enumerate(ort_n):
        print(f"P_{i}(t) =", simplify(v[0]))
    print(f"es ortogonal? ->", is_ort_base_custom(ort_n, func_product, numeric=True))

    # c) modificar para peso sqrt(1-t^2) (Chebyshev) y encuentra primeros n
    n = 3
    ort_n_cheb = gram_schmidt(n)
    print(f"\nc) Primeros {n} polinomios ortogonales con peso sqrt(1-t^2):")
    for i, v in enumerate(ort_n_cheb):
        print(f"T_{i}(t) =", simplify(v[0]))
    print(f"es ortogonal? ->", is_ort_base(ort_n_cheb, numeric=True))

    # d) suponga la funcion h(t) = sen(3t)(1-t^2)
    h = sin(3*t)*(1-t**2)
    h_matrix = Matrix([h])
    t_vals = np.linspace(-1, 1, 400)

    # I. Expanda la función h(x) en términos de la base de monomios y de polinomios de
    #Legendre, grafique, compare y encuentre el grado de los polinomios en los cuales
    #difieren las expansiones.
    n = 5
    legendre_approx, coeffs_leg = expansion_in_legendre(h_matrix, n)
    print("coef Legendre:", coeffs_leg[:5])
    print(f"\nExpansión de h(t) en términos de los primeros {n} polinomios de Legendre:")
    print(legendre_approx[:5])

    f_h   = lambdify(t, h, 'numpy')
    f_leg = lambdify(t, legendre_approx[0], 'numpy')

    # II. Expanda la función h(x) en términos de la base de monomios y de polinomios de
    #Chebyshev, grafique, compare y encuentre el grado de los polinomios en los cuales
    #difieren las expansiones.
    chebyshev_approx, coeffs_cheb = expansion_in_chebyshev(h_matrix, n)
    print("coef Chebyshev:", coeffs_cheb[:5])
    print(f"\nExpansión de h(t) en términos de los primeros {n} polinomios de Chebyshev:")
    print(chebyshev_approx[:5])

    f_cheb = lambdify(t, chebyshev_approx[0], 'numpy')

    # III. Expanda la función h(x) en términos de la base de monomios y de polinomios de
    #legendre y de chebyshev, grafique, compare y encuentre el grado de los polinomios en los cuales
    #difieren las expansiones.

    # IV. Estime en cada caso el error que se comete como función del grado del polinomio (o
    #monomio) de la expansión.
    grados = list(range(1, 11))  # hasta grado 10
    errores_leg, errores_cheb = [], []

    for k in grados:
        leg_k, _ = expansion_in_legendre(h_matrix, k)
        cheb_k, _ = expansion_in_chebyshev(h_matrix, k)
        f_leg_k  = lambdify(t, leg_k[0], 'numpy')
        f_cheb_k = lambdify(t, cheb_k[0], 'numpy')
        x = np.linspace(-1, 1, 2001)
        h_vals = f_h(x)
        err_leg  = np.sqrt(np.trapz((h_vals - f_leg_k(x))**2, x))
        err_cheb = np.sqrt(np.trapz((h_vals - f_cheb_k(x))**2, x))
        errores_leg.append(err_leg)
        errores_cheb.append(err_cheb)

    # --- Subplots (I, II, III, IV) ---
    fig, axs = plt.subplots(4, 1, figsize=(9,16), sharex=False)

    # I: h vs Legendre
    axs[0].plot(t_vals, f_h(t_vals),   label='h(t)', color='black', linewidth=2)
    axs[0].plot(t_vals, f_leg(t_vals), label=f'Legendre (n={n})', color='C0', linestyle='--')
    axs[0].set_title('I. h(t) vs Legendre')
    axs[0].axhline(0, color='gray', lw=0.6, ls='--')
    axs[0].legend(); axs[0].grid(True)

    # II: h vs Chebyshev
    axs[1].plot(t_vals, f_h(t_vals),    label='h(t)', color='black', linewidth=2)
    axs[1].plot(t_vals, f_cheb(t_vals), label=f'Chebyshev (n={n})', color='C1', linestyle='--')
    axs[1].set_title('II. h(t) vs Chebyshev')
    axs[1].axhline(0, color='gray', lw=0.6, ls='--')
    axs[1].legend(); axs[1].grid(True)

    # III: Legendre vs Chebyshev
    axs[2].plot(t_vals, f_leg(t_vals),  label=f'Legendre (n={n})', color='C0', linestyle='--')
    axs[2].plot(t_vals, f_cheb(t_vals), label=f'Chebyshev (n={n})', color='C1', linestyle='-.')
    axs[2].set_title('III. Legendre vs Chebyshev')
    axs[2].axhline(0, color='gray', lw=0.6, ls='--')
    axs[2].legend(); axs[2].grid(True)

    # IV: Errores
    axs[3].plot(grados, errores_leg, 'o-', label='Error Legendre')
    axs[3].plot(grados, errores_cheb, 's--', label='Error Chebyshev')
    axs[3].set_title('IV. Error de aproximación en función del grado')
    axs[3].set_xlabel('Grado n')
    axs[3].set_ylabel('Error L2')
    axs[3].grid(True); axs[3].legend()

    plt.tight_layout()
    plt.show()
