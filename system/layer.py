from math import inf
from typing import Optional, Union

import numericalunits as nu

class Layer:
    def __init__(self, 
                 name: str,
                 thickness: Union[float, str],
                 ex: complex,
                 ez: Optional[complex] = None,
                 mu: float = 1.0):
        """
        Класс для представления слоя в оптической системе.
        
        Args:
            name: Название слоя
            thickness: Толщина слоя в нм или 'inf' для бесконечной толщины
            ex: Диэлектрическая проницаемость по x (может быть комплексной)
            ez: Диэлектрическая проницаемость по z (по умолчанию равна ex)
            mu: Магнитная проницаемость (по умолчанию 1)
        """
        self.name = name
        self._thickness = thickness
        self.ex = ex
        self.ez = ez if ez is not None else ex
        self.mu = mu
        
    @property
    def thickness(self) -> float:
        """Возвращает толщину слоя в метрах."""
        if isinstance(self._thickness, str) and self._thickness.lower() == 'inf':
            return inf
        return self._thickness * nu.nm
    
    @property
    def is_infinite(self) -> bool:
        """Проверяет, является ли слой бесконечно толстым."""
        return self.thickness == inf
    
    def __repr__(self) -> str:
        return (f"Layer(name='{self.name}', thickness={self._thickness}, "
                f"ex={self.ex}, ez={self.ez}, mu={self.mu})")