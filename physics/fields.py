from __future__ import division, print_function

from physics.modes import find_layer

import cmath
from typing import Dict, Optional

inf = float('inf')

__all__ = ['calculate_Hy', 'calculate_Ex', 'calculate_Ez', 'calculate_Sx']

def calculate_Hy(z: float, params: Dict, x: float = 0, layer: Optional[int] = None):
    """
    Рассчитывает комплексное значение H-поля (компонента y) в точке (x, z).

    Args:
        z (float): z-координата точки.
        params (dict): параметры моды, содержащий:
            - d_list (list): список толщин слоёв,
            - H_up_list (list): амплитуды волны, идущей вверх,
            - H_down_list (list): амплитуды волны, идущей вниз,
            - kz_list (list): волновые числа в каждом слое,
            - kx (complex): горизонтальное волновое число,
            - layer_bottom_list (list): позиции нижних границ слоёв.
        x (float): x-координата точки (по умолчанию 0).
        layer (int or None): номер слоя для расчёта поля.
            Если None, определяется автоматически.

    Returns:
        complex: Значение Hy в указанной точке.
    """
    N = len(params['d_list'])
    layer_idx = layer if layer is not None else find_layer(z, params)
    
    H_up = params['H_up_list'][layer_idx]
    H_down = params['H_down_list'][layer_idx]
    kz = params['kz_list'][layer_idx]
    kx = params['kx']

    layer_bottom = params['layer_bottom_list'][layer_idx]
    layer_top = inf if layer_idx == N - 1 else params['layer_bottom_list'][layer_idx + 1]

    up_term = _get_up_term(H_up, kz, z, layer_bottom, x, kx)
    down_term = _get_down_term(H_down, kz, layer_top, z, x, kx)

    return up_term + down_term


def _get_up_term(amplitude, kz, z, bottom, x, kx):
    """Вспомогательная функция для верхней компоненты."""
    if amplitude == 0:
        return 0
    return amplitude * cmath.exp(1j * kz * (z - bottom) + 1j * kx * x)


def _get_down_term(amplitude, kz, top, z, x, kx):
    """Вспомогательная функция для нижней компоненты."""
    if amplitude == 0:
        return 0
    return amplitude * cmath.exp(1j * kz * (top - z) + 1j * kx * x)

def calculate_Ex(z: float, params: Dict, x: float = 0, layer: Optional[int] =None):
    """
    Рассчитывает x-компоненту электрического поля в точке (x, z).

    Args:
        z (float): z-координата точки.
        params (dict): параметры моды, содержащий:
            - d_list (list): список толщин слоёв,
            - Ex_up_list (list): амплитуды волны, идущей вверх,
            - Ex_down_list (list): амплитуды волны, идущей вниз,
            - kz_list (list): волновые числа в каждом слое,
            - layer_bottom_list (list): позиции нижних границ слоёв.
        x (float): x-координата точки (по умолчанию 0).
        layer (int or None): номер слоя для расчёта поля.
            Если None, определяется автоматически.

    Returns:
        complex: Значение Ex в указанной точке.
    """
    N = len(params['d_list'])
    layer_idx = layer if layer is not None else find_layer(z, params)

    Ex_up = params['Ex_up_list'][layer_idx]
    Ex_down = params['Ex_down_list'][layer_idx]
    kz = params['kz_list'][layer_idx]
    kx = params['kx']

    layer_bottom = params['layer_bottom_list'][layer_idx]
    layer_top = inf if layer_idx == N - 1 else params['layer_bottom_list'][layer_idx + 1]

    up_term = _get_up_term(Ex_up, kz, z, layer_bottom, x, kx)
    down_term = _get_down_term(Ex_down, kz, layer_top, z, x, kx)

    return up_term + down_term

def calculate_Ez(z: float, params: Dict, x: float = 0, layer: Optional[int] = None):
    """
    Рассчитывает z-компоненту электрического поля в точке (x, z).

    Args:
        z (float): z-координата точки.
        params (dict): параметры моды, содержащий:
            - Ez_up_list (list): амплитуды Ez, идущей вверх,
            - Ez_down_list (list): амплитуды Ez, идущей вниз.
        x (float): x-координата точки (по умолчанию 0).
        layer (int or None): номер слоя для расчёта поля.
            Если None, определяется автоматически.

    Returns:
        complex: Значение Ez в указанной точке.
    """
    N = len(params['d_list'])
    layer_idx = layer if layer is not None else find_layer(z, params)

    Ez_up = params['Ez_up_list'][layer_idx]
    Ez_down = params['Ez_down_list'][layer_idx]
    kz = params['kz_list'][layer_idx]
    kx = params['kx']

    layer_bottom = params['layer_bottom_list'][layer_idx]
    layer_top = inf if layer_idx == N - 1 else params['layer_bottom_list'][layer_idx + 1]

    up_term = _get_up_term(Ez_up, kz, z, layer_bottom, x, kx)
    down_term = _get_down_term(Ez_down, kz, layer_top, z, x, kx)

    return up_term + down_term

def calculate_Sx(z: float, params: Dict, x: float = 0, layer: Optional[int] = None):
    """
    Рассчитывает x-компоненту вектора Пойнтинга в точке (x, z).

    Args:
        z (float): z-координата точки.
        params (dict): параметры моды, содержащий:
            - Ez_up_list (list), Ez_down_list (list),
            - H_up_list (list), H_down_list (list),
            - kz_list (list), kx (complex),
            - layer_bottom_list (list).
        x (float): x-координата точки (по умолчанию 0).
        layer (int or None): номер слоя для расчёта.
            Если None, определяется автоматически.

    Returns:
        complex: Значение Sx в указанной точке.
    """
    Ez_here = calculate_Ez(z, params, x=x, layer=layer)
    Hy_here = calculate_Hy(z, params, x=x, layer=layer)
    return -0.5 * Ez_here * Hy_here.conjugate()