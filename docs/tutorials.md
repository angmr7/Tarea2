# Tutoriales

## Ejemplo de uso

En el siguiente ejemplo queremos resolver la integral 

$$\int_0^3 (x^6 - x^2 sin(2x)) dx$$

---

Se empieza eligiendo el valor de N para el cual se quiere resolver, en este caso escojamos N=4

Se calculan los puntos de muestreo y los pesos

```python
x0_points, weights = gaussxw(2)
```

Se escalan los puntos de muestreo y los pesos al intervalo de integracion [1,3] para este caso


```python
scaled_x0_points, scaled_weights = gaussxwab(1, 3, x0_points, weights)
```

Se evalúa la función en los puntos escalados y se calcula la integral aproximada, haciendo uso de la ecuación.

$$\int_a^b {\rm{d}}x f(x) \approx \sum_{k=1}^{N+1} w_k f(x_k)$$


```python
result = np.sum(scaled_weights * int_function(scaled_x0_points))
```

