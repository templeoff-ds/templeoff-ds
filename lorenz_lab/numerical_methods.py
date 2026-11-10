"""Самостоятельные явные одношаговые методы для векторной задачи Коши."""
from dataclasses import dataclass
import numpy as np

@dataclass
class Solution:
    t: np.ndarray
    y: np.ndarray  # строки — моменты времени, столбцы — компоненты
    n_steps: int
    nfev: int

def euler_step(fun, t, u, h):
    return u + h*fun(t, u)

def rk4_step(fun, t, u, h):
    # Половинные стадии и веса 1,2,2,1 / 6 — классическая таблица RK4.
    k1 = fun(t, u)
    k2 = fun(t+h/2, u+h*k1/2)
    k3 = fun(t+h/2, u+h*k2/2)
    k4 = fun(t+h, u+h*k3)
    return u + h*(k1+2*k2+2*k3+k4)/6

def time_grid(t_start, t_end, h):
    if not np.all(np.isfinite([t_start, t_end, h])) or h <= 0 or t_end <= t_start:
        raise ValueError('Нужны конечные t_end>t_start и h>0')
    n = int(round((t_end-t_start)/h))
    # Эти допуски проверяют геометрию сетки, а не ошибку интегрирования.
    if n < 1 or not np.isclose(n*h, t_end-t_start, rtol=1e-11, atol=1e-13):
        raise ValueError('Длина интервала должна делиться на шаг; сетка не усекается')
    return np.linspace(t_start, t_end, n+1)

def integrate(fun, initial_state, t_start, t_end, h, method='rk4'):
    t = time_grid(t_start, t_end, h)
    u = np.asarray(initial_state, dtype=float)
    if u.ndim != 1 or not np.all(np.isfinite(u)):
        raise ValueError('Начальное состояние — конечный одномерный вектор')
    if method not in ('euler', 'rk4'):
        raise ValueError('Неизвестный метод: '+method)
    step = euler_step if method == 'euler' else rk4_step
    y = np.empty((len(t), len(u))); y[0] = u
    with np.errstate(over='raise', invalid='raise', divide='raise'):
        for i in range(len(t)-1):
            y[i+1] = step(fun, t[i], y[i], h)
            if not np.all(np.isfinite(y[i+1])):
                raise FloatingPointError(f'Неконечное состояние на шаге {i+1}')
    n = len(t)-1
    return Solution(t, y, n, n*(1 if method == 'euler' else 4))
