from copy import deepcopy
from dataclasses import dataclass
from typing import Any, Dict, List, Optional

from system.layer import Layer

INF = float('inf')
@dataclass
class OpticalSystemConfig:
    """
    Конфигурация оптической системы.
    
    Attributes:
        w: Угловая частота (рад/с)
        layers: Список слоев системы
    """
    w: float
    layers: List[Layer]
    
    def __post_init__(self):
        """Проверка корректности конфигурации после инициализации."""
        assert len(self.layers) >= 2, "Система должна содержать минимум 2 слоя"
        assert self.layers[0].is_infinite and self.layers[-1].is_infinite, \
            "Первый и последний слои должны быть бесконечными"
    
    @property
    def d_list(self) -> List[float]:
        """Список толщин слоев в метрах."""
        return [layer.thickness for layer in self.layers]
    
    @property
    def ex_list(self) -> List[complex]:
        """Список диэлектрических проницаемостей по x."""
        return [layer.ex for layer in self.layers]
    
    @property
    def ez_list(self) -> List[complex]:
        """Список диэлектрических проницаемостей по z."""
        return [layer.ez for layer in self.layers]
    
    @property
    def mu_list(self) -> List[float]:
        """Список магнитных проницаемостей."""
        return [layer.mu for layer in self.layers]
    
    @property
    def layer_bottom_list(self) -> List[float]:
        """Список нижних границ слоев в метрах."""
        bottoms = [-INF, 0.0]
        for layer in self.layers[1:-1]:
            bottoms.append(bottoms[-1] + layer.thickness)
        return bottoms
    
    def __getitem__(self, key: str) -> Any:
        """Позволяет получать параметры через config['w'], config['d_list'] и т.д."""
        if hasattr(self, key):
            return getattr(self, key)
        elif key in ['d_list', 'ex_list', 'ez_list', 'mu_list', 'layer_bottom_list']:
            return getattr(self, key)
        raise KeyError(f"У OpticalSystemConfig нет поля {key}")

class OpticalState:
    """Состояние оптической системы (поля и волновые числа)."""
    def __init__(self, 
                 kx: Optional[complex] = None,
                 kz_list: Optional[List[complex]] = None,
                 H_up_list: List[complex] = None,
                 H_down_list: List[complex] = None,
                 Ex_up_list: List[complex] = None,
                 Ex_down_list: List[complex] = None,
                 Ez_up_list: List[complex] = None,
                 Ez_down_list: List[complex] = None,
                 Sx_list: List[complex] = None,
                 Sx_total: Optional[complex] = None):
        self.kx = kx
        self.kz_list = kz_list or []
        self.H_up_list = H_up_list or []
        self.H_down_list = H_down_list or []
        self.Ex_up_list = Ex_up_list or []
        self.Ex_down_list = Ex_down_list or []
        self.Ez_up_list = Ez_up_list or []
        self.Ez_down_list = Ez_down_list or []
        self.Sx_list = Sx_list or []
        self.Sx_total = Sx_total
        
    def __getitem__(self, key: str) -> Any:
        if hasattr(self, key):
            return getattr(self, key)
        raise KeyError(f"У OpticalState нет поля {key}")
        
    def __setitem__(self, key: str, value: Any):
        if hasattr(self, key):
            setattr(self, key, value)
        else:
            raise KeyError(f"У OpticalState нет поля {key}")

class OpticalSystem:
    """Многослойная оптическая система для моделирования поверхностных мод."""
    def __init__(self, config: OpticalSystemConfig, state: Optional[OpticalState] = None):
        self.config = config
        self.state = state or OpticalState()
        
    def __getitem__(self, key: str) -> Any:
        """Позволяет получать параметры через system['w'], system['d_list'] и т.д."""
        if hasattr(self.config, key):
            return getattr(self.config, key)
        if hasattr(self.state, key):
            return getattr(self.state, key)
        raise KeyError(f"У OpticalSystem нет поля {key}")
        
    def __setitem__(self, key: str, value: Any):
        """Позволяет устанавливать параметры через system['w'] = ..."""
        if hasattr(self.config, key):
            setattr(self.config, key, value)
        elif hasattr(self.state, key):
            setattr(self.state, key, value)
        else:
            raise KeyError(f"У OpticalSystem нет поля {key}")
        
    def copy(self) -> 'OpticalSystem':
        """Возвращает глубокую копию текущего объекта OpticalSystem."""
        return OpticalSystem(
            config=deepcopy(self.config),
            state=deepcopy(self.state)
        )
    
    def update(self, data: Dict[str, Any]) -> 'OpticalSystem':
        """Обновляет атрибуты системы из словаря."""
        for key, value in data.items():
            self[key] = value
        return self