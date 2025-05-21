from typing import List, Any

from copy import deepcopy

class OpticalSystem:
    def __init__(self, w: float, d_list: List[float], ex_list: List[complex],
                 ez_list: List[complex], mu_list: List[complex]=None):
        """
        Определяет многослойную оптическую систему для моделирования поверхностных мод.

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
            
        # === Динамически добавляемые поля ===
        self.kx = None
        self.kz_list = None
        self.H_up_list = []
        self.H_down_list = []
        self.Ex_up_list = []
        self.Ex_down_list = []
        self.Ez_up_list = []
        self.Ez_down_list = []
        self.Sx_list = []
        self.Sx_total = None
        
    def __getitem__(self, key: str) -> Any:
        """
        Позволяет получать параметры через system['w'], system['d_list'] и т.д.
        """
        if hasattr(self, key):
            return getattr(self, key)
        else: 
            raise KeyError(f"У OpticalSystem нет поля {key}")
        
    def __setitem__(self, key: str, value: Any):
        """
        Позволяет устанавливать параметры через system['w'] = ...
        """
        if hasattr(self, key):
            setattr(self, key, value)
        else:
            raise KeyError(f"У OpticalSystem нет поля {key}")
        
    def copy(self):
        """
        Возвращает глубокую копию текущего объекта OpticalSystem.
        """
        copied = OpticalSystem(
            w=self.w,
            d_list=deepcopy(self.d_list),
            ex_list=deepcopy(self.ex_list),
            ez_list=deepcopy(self.ez_list),
            mu_list=deepcopy(self.mu_list)
        )
        copied.kz_list = deepcopy(self.kz_list)
        copied.layer_bottom_list = deepcopy(self.layer_bottom_list)

        copied.H_up_list = deepcopy(self.H_up_list)
        copied.H_down_list = deepcopy(self.H_down_list)
        copied.Ex_up_list = deepcopy(self.Ex_up_list)
        copied.Ex_down_list = deepcopy(self.Ex_down_list)
        copied.Ez_up_list = deepcopy(self.Ez_up_list)
        copied.Ez_down_list = deepcopy(self.Ez_down_list)
        copied.Sx_list = deepcopy(self.Sx_list)
        copied.Sx_total = deepcopy(self.Sx_total)

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
            if hasattr(self, key):
                setattr(self, key, value)
            else:
                raise KeyError(f"У OpticalSystem нет поля {key}")
        return self