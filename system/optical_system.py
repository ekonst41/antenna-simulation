from typing import List, Any

from copy import deepcopy

class OpticalSystemConfig:
    def __init__(self, w: float, d_list: List[float], ex_list: List[complex],
                 ez_list: List[complex], mu_list: List[complex]=None):
        """
        Определяет изначальную конфигурацию системы.

        Args:
            w (float): угловая частота (рад/с)
            d_list (list): толщины слоёв (м), первый и последний должны быть inf
            ex_list (list): ε_x каждого слоя
            ez_list (list): ε_z каждого слоя
            mu_list (list or None): μ_y каждого слоя (по умолчанию [1] * N)
        """
        self.w = w
        self.d_list = d_list
        self.ex_list = ex_list
        self.ez_list = ez_list
        self.mu_list = mu_list if mu_list is not None else [1.0] * len(ex_list)
        self.N = len(self.d_list)
        
        assert len(self.ex_list) == self.N and len(self.ez_list) == self.N and len(self.mu_list) == self.N
        assert self.N >= 2
        assert self.d_list[0] == float('inf') and self.d_list[-1] == float('inf')
        
        self.layer_bottom_list = [-float('inf'), 0]
        for i in range(1, self.N - 1):
            self.layer_bottom_list.append(self.layer_bottom_list[-1] + self.d_list[i])
            
        
    def __getitem__(self, key: str) -> Any:
        if hasattr(self, key):
            return getattr(self, key)
        else: 
            raise KeyError(f"У OpticalSystemConfig нет поля {key}")
        
    def __setitem__(self, key: str, value: Any):
        if hasattr(self, key):
            setattr(self, key, value)
        else:
            raise KeyError(f"У OpticalSystemConfig нет поля {key}")
        
class OpticalState:
    def __init__(self, kx=None, kz_list=None, H_up_list=[], H_down_list=[],
                 Ex_up_list=[], Ex_down_list=[], Ez_up_list=[], Ez_down_list=[],
                 Sx_list=[], Sx_total=None):
        self.kx = kx
        self.kz_list = kz_list
        self.H_up_list = H_up_list
        self.H_down_list = H_down_list
        self.Ex_up_list = Ex_up_list
        self.Ex_down_list = Ex_down_list
        self.Ez_up_list = Ez_up_list
        self.Ez_down_list = Ez_down_list
        self.Sx_list = Sx_list
        self.Sx_total = Sx_total
        
    def __getitem__(self, key: str) -> Any:
        if hasattr(self, key):
            return getattr(self, key)
        else: 
            raise KeyError(f"У OpticalState нет поля {key}")
        
    def __setitem__(self, key: str, value: Any):
        if hasattr(self, key):
            setattr(self, key, value)
        else:
            raise KeyError(f"У OpticalState нет поля {key}")

class OpticalSystem:
    def __init__(self, config: OpticalSystemConfig, state: OpticalState=None):
        """
        Определяет многослойную оптическую систему для моделирования поверхностных мод.

        Args:
            config (OpticalSystemConfig): изначальная конфигурация системы
            state (OpticalState): состояение полей в ситемы
        """
        assert config is not None
        self.config = config
        self.state = state if state else OpticalState()
        
    def __getitem__(self, key: str) -> Any:
        """
        Позволяет получать параметры через system['w'], system['d_list'] и т.д.
        """
        if hasattr(self.config, key):
            return getattr(self.config, key)
        elif hasattr(self.state, key):
            return getattr(self.state, key)
        else: 
            raise KeyError(f"У OpticalSystem нет поля {key}")
        
    def __setitem__(self, key: str, value: Any):
        """
        Позволяет устанавливать параметры через system['w'] = ...
        """
        if hasattr(self.config, key):
            setattr(self.config, key, value)
        elif hasattr(self.state, key):
            setattr(self.state, key, value)
        else:
            raise KeyError(f"У OpticalSystem нет поля {key}")
        
    def copy(self):
        """
        Возвращает глубокую копию текущего объекта OpticalSystem.
        """
        copied = OpticalSystem(
            config=deepcopy(self.config),
            state=deepcopy(self.state)
        )
        return copied
    
    def update(self, data: dict):
        """
        Обновляет атрибуты системы из словаря.
        
        Args:
            data (dict): словарь с новыми значениями полей и параметров
            
        Raises:
            KeyError: если передан ключ, которого нет в классе
        """
        for key, value in data.items():
            if hasattr(self.config, key):
                setattr(self.config, key, value)
            elif hasattr(self.state, key):
                setattr(self.state, key, value)
            else:
                raise KeyError(f"У OpticalSystem нет поля {key}")
        return self