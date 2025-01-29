import numpy as np

def gaussxw(N):
    """Calculate sample points and weights for Gaussian quadrature.

    This function computes the sample points and weights for the
    Gaussian quadrature of order N. It uses the Newton's method to
    find the roots of the Legendre polynomial.

    Examples:
        >>> gaussxw(2)
        (array([0.57735027, 0.57735027]), array([1., 1.]))

    Args:
        N (int): The order of the Gaussian quadrature.

    Returns:
        tuple: A tuple containing two numpy arrays:
            - x (numpy.ndarray): The sample points.
            - w (numpy.ndarray): The weights corresponding to the sample points.
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
    """Scale sample points and weights to a specific interval.

    This function scales the sample points and weights obtained from
    the Gaussian quadrature to the interval [a, b].

    Examples:
        >>> gaussxwab(0,2, x0_points_N_2, weights_N_2)
        (array([1.57735027, 1.57735027]), array([1., 1.]))

    Args:
        a (float): The start of the interval.
        b (float): The end of the interval.
        x (numpy.ndarray): The sample points.
        w (numpy.ndarray): The weights.

    Returns:
        tuple: A tuple containing two numpy arrays:
            - scaled_x (numpy.ndarray): The scaled sample points.
            - scaled_w (numpy.ndarray): The scaled weights.
    """
    return 0.5 * (b - a) * x + 0.5 * (b + a), 0.5 * (b - a) * w

def int_function(x):
    """Compute the value of the integrand function.

    This function computes the value of the function to be integrated,
    which is x^6 - x^2 * sin(2 * x).

    Examples:
        >>> int_function(3)
        731.5147394837903

    Args:
        x (float or numpy.ndarray): The input value(s).

    Returns:
        float or numpy.ndarray: The computed value of the function.
    """
    return x**6 - x**2 * np.sin(2 * x)

N_values = [2, 3, 4, 5, 11, 12, 13]

for N in N_values:
    x0_points, weights = gaussxw(N)
    scaled_x0_points, scaled_weights = gaussxwab(1, 3, x0_points, weights)
    result = np.sum(scaled_weights * int_function(scaled_x0_points))
    print(f"Resultado de la integración con N={N}, es {result}")

