from physics.fields import calculate_Ex, calculate_Ez, calculate_Hy, calculate_Sx
import numericalunits as nu
import scipy
import math, cmath
import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy
inf = float("inf")

def floats_are_equal(a, b, tol=1e-5):
    """
    Проверяет, равны ли два числа с заданной относительной точностью.

    Args:
        a (complex or float): первое число
        b (complex or float): второе число
        tol (float): относительная погрешность

    Returns:
        bool: True, если числа совпадают в пределах допуска
    """
    return abs(a - b) <= tol * (abs(a) + abs(b))


def assert_floats_are_equal(a, b, tol=1e-5):
    """
    Проверяет, равны ли два числа, иначе вызывает AssertionError.

    Args:
        a (complex or float): первое число
        b (complex or float): второе число
        tol (float): относительная погрешность
    """
    assert floats_are_equal(a, b, tol), f"{a} != {b}"
    
from physics.fields import calculate_Ex, calculate_Ez, calculate_Hy


def _check_boundary_conditions(params, tol=1e-5):
    """
    Проверяет выполнение граничных условий для Ex, Ez, Hy между слоями.

    Args:
        params (dict): параметры моды
        tol (float): относительная погрешность сравнения

    Returns:
        str or True: сообщение об ошибке или True
    """
    N = len(params['d_list'])

    for i in range(N - 1):
        z = params['layer_bottom_list'][i + 1]  # граница между слоями i и i+1

        # Проверяем Ex
        ex_under = calculate_Ex(z, params, layer=i)
        ex_over = calculate_Ex(z, params, layer=i + 1)
        if not floats_are_equal(ex_under, ex_over, tol):
            return f'Ex не совпадает на границе {i}-{i+1}: {ex_under}, {ex_over}'

        # Проверяем Ez
        ez_under = calculate_Ez(z, params, layer=i)
        ez_over = calculate_Ez(z, params, layer=i + 1)
        if not floats_are_equal(ez_under, ez_over, tol):
            return f'Ez не совпадает на границе {i}-{i+1}: {ez_under}, {ez_over}'

        # Проверяем Hy
        hy_under = calculate_Hy(z, params, layer=i)
        hy_over = calculate_Hy(z, params, layer=i + 1)
        if not floats_are_equal(hy_under, hy_over, tol):
            return f'Hy не совпадает на границе {i}-{i+1}: {hy_under}, {hy_over}'

    return True

def _check_kz_relations(params, tol=1e-8):
    """
    Проверяет соотношение kz^2 == w²*μ*ε/c0² - kx²*ε/ε_z.

    Args:
        params (dict): параметры моды
        tol (float): допуск проверки

    Returns:
        None or raises: ошибка при неверном соотношении
    """
    w = params['w']
    kx = params['kx']
    kz_list = params['kz_list']
    ex_list = params['ex_list']
    ez_list = params['ez_list']
    mu_list = params['mu_list']

    for i in range(len(kz_list)):
        kz = kz_list[i]
        ex = ex_list[i]
        ez = ez_list[i]
        mu = mu_list[i]

        lhs = kz ** 2
        rhs = (w ** 2 * mu * ex / nu.c0 ** 2) - (kx ** 2 * ex / ez)

        assert_floats_are_equal(lhs, rhs, tol=tol)
        assert kz.imag >= 0, f"Im(kz) < 0 @ слой {i}"

        if (i == 0 or i == len(kz_list) - 1) and kz.imag == 0:
            raise ValueError(f"kz не затухает @ слой {i}")

    return True

from physics.fields import calculate_Sx


def _check_poynting_vector(params, tol=1e-5):
    """
    Проверяет согласованность вектора Пойнтинга через интегрирование.

    Args:
        params (dict): параметры моды
        tol (float): допуск проверки

    Returns:
        None or raises: ошибка при несоответствии
    """
    Sx_list = params['Sx_list']
    N = len(Sx_list)
    scale_factor = max(abs(calculate_Sx(0, params, layer=0)), abs(calculate_Sx(0, params, layer=1)))
    assert scale_factor != 0

    for i in range(N):
        lower_z = params['layer_bottom_list'][i] if i != 0 else -20 / abs(params['kz_list'][i].imag)
        upper_z = params['layer_bottom_list'][i + 1] if i != N - 1 else 20 / abs(params['kz_list'][i].imag)

        def integrand_re(z):
            return (calculate_Sx(z, params, layer=i) / scale_factor).real

        def integrand_im(z):
            return (calculate_Sx(z, params, layer=i) / scale_factor).imag

        integral_re, _ = scipy.integrate.quad(integrand_re, lower_z, upper_z)
        integral_im, _ = scipy.integrate.quad(integrand_im, lower_z, upper_z)

        Sx_integrated = (integral_re + 1j * integral_im) * scale_factor
        Sx_expected = Sx_list[i]

        assert floats_are_equal(Sx_integrated, Sx_expected, tol), \
            f"Sx не совпадает @ слой {i}: {Sx_integrated} != {Sx_expected}"

    total_Sx = sum(Sx_list)
    assert floats_are_equal(total_Sx, params['Sx_total'], tol=1e-8), \
        "Sx_total ≠ сумме по слоям"

    return True

def check_mode(params, thorough=False, tol=1e-5):
    """
    Проверяет корректность найденной моды.

    Args:
        params (dict): параметры моды
        thorough (bool): провести полную проверку (медленнее)
        tol (float): относительная погрешность

    Returns:
        True или строка с первой обнаруженной ошибкой
    """
    bc_result = _check_boundary_conditions(params, tol)
    if bc_result is not True:
        return bc_result

    _check_kz_relations(params)

    if thorough:
        _check_poynting_vector(params, tol)

    return True
