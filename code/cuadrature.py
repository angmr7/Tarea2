import numpy as np

def gaussxw(N):
    """Calcular puntos de muestra y pesos para la cuadratura gaussiana.

    Esta función calcula los puntos de muestra y los pesos para la
    cuadratura gaussiana de orden N. Utiliza el método de Newton para
    encontrar las raíces del polinomio de Legendre.

    Examples:
        >>> gaussxw(2)
        (array([0.57735027, 0.57735027]), array([1., 1.]))

    Args:
        N (int): El orden de la cuadratura gaussiana.

    Returns:
        tuple: Una tupla que contiene dos arreglos de numpy:
            - x (numpy.ndarray): Los puntos de muestra.
            - w (numpy.ndarray): Los pesos correspondientes a los puntos de muestra.
    """
    a = np.linspace(3, 4 * (N - 1), N) / ((4 * N) + 2)
    x = np.cos(np.pi * a + 1 / (8 * N * N * np.tan(a)))

    epsilon = 1e-15
    delta = 1.0
    while delta > epsilon:
        p0 = np.ones(N, dtype=float)
        p1 = np.copy(x)
        for k in range(1, N):
            p0, p1 = p1, ((2 * k + 1) * x * p1 - k * p0) / (k + 1)
        dp = (N + 1) * (p0 - x * p1) / (1 - x * x)
        dx = p1 / dp
        x -= dx
        delta = np.max(np.abs(dx))

    w = 2 * (N + 1) * (N + 1) / (N * N * (1 - x * x) * dp * dp)
    return x, w

def gaussxwab(a, b, x, w):
    """Escalar puntos de muestra y pesos a un intervalo específico.

    Esta función escala los puntos de muestra y los pesos obtenidos de
    la cuadratura gaussiana al intervalo [a, b].

    Examples:
        >>> gaussxwab(0,2, x0_points_N_2, weights_N_2)
        (array([1.57735027, 1.57735027]), array([1., 1.]))

    Args:
        a (float): El inicio del intervalo.
        b (float): El final del intervalo.
        x (numpy.ndarray): Los puntos de muestra.
        w (numpy.ndarray): Los pesos.

    Returns:
        tuple: Una tupla que contiene dos arreglos de numpy:
            - scaled_x (numpy.ndarray): Los puntos de muestra escalados.
            - scaled_w (numpy.ndarray): Los pesos escalados.
    """
    return 0.5 * (b - a) * x + 0.5 * (b + a), 0.5 * (b - a) * w

def int_function(x):
    """Calcular el valor de la función integrando.

    Esta función calcula el valor de la función a integrar,
    que es x^6 - x^2 * sin(2 * x).

    Examples:
        >>> int_function(3)
        731.5147394837903

    Args:
        x (float o numpy.ndarray): El valor o valores de entrada.

    Returns:
        float o numpy.ndarray: El valor calculado de la función.
    """
    return x**6 - x**2 * np.sin(2 * x)

N_values = [2, 3, 4, 5, 11, 12, 13]

for N in N_values:
    x0_points, weights = gaussxw(N)
    scaled_x0_points, scaled_weights = gaussxwab(1, 3, x0_points, weights)
    result = np.sum(scaled_weights * int_function(scaled_x0_points))
    print(f"Resultado de la integración con N={N}, es {result}")

