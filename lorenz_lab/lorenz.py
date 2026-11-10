"""Правая часть, Якобиан и равновесия без зависимости от интегратора."""
from dataclasses import dataclass
import numpy as np

@dataclass(frozen=True)
class Parameters:
    sigma: float = 10.0
    rho: float = 28.0
    beta: float = 8.0 / 3.0

    def __post_init__(self):
        if not np.all(np.isfinite([self.sigma, self.rho, self.beta])):
            raise ValueError('Параметры должны быть конечными')
        if self.sigma <= 0 or self.beta <= 0:
            raise ValueError('В этой работе нужны sigma>0 и beta>0')

def rhs(t, state, params=Parameters()):
    """t не используется: система автономная. Возвращается новый массив."""
    x, y, z = state
    return np.array([params.sigma * (y-x), x*(params.rho-z)-y,
                     x*y-params.beta*z], dtype=float)

def jacobian(state, params=Parameters()):
    x, y, z = state
    return np.array([[-params.sigma, params.sigma, 0.0],
                     [params.rho-z, -1.0, -x],
                     [y, x, -params.beta]])

def equilibria(params=Parameters()):
    result = {'O': np.zeros(3)}
    if params.rho > 1:
        a = np.sqrt(params.beta*(params.rho-1))
        result.update({'C+': np.array([a, a, params.rho-1]),
                       'C-': np.array([-a, -a, params.rho-1])})
    return result
