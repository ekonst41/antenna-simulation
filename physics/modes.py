from __future__ import division, print_function
import numpy as np
import math, cmath
import scipy.optimize
import scipy.integrate
import matplotlib.pyplot as plt
from math import pi
from copy import deepcopy
import numericalunits as nu

from tmm.tmm import bc_matrix

inf = float('inf')


def find_kzs(params):
    """
    Рассчитывает вертикальное волновое число kz для каждого слоя.

    Args:
        params (dict): параметры системы, включая:
            - w (float): угловая частота,
            - kx (complex): горизонтальное волновое число,
            - ex_list (list): комплексная проницаемость по x,
            - ez_list (list): комплексная проницаемость по z,
            - mu_list (list): магнитная проницаемость.

    Returns:
        dict: обновлённый params с добавленным 'kz_list'
    """
    w = params['w']
    kx = params['kx']
    ex_list = params['ex_list']
    ez_list = params['ez_list']
    mu_list = params['mu_list']
    N = len(ez_list)

    assert N == len(ex_list) == len(ez_list) == len(mu_list) >= 2
    assert w > 0
    for list_name in ['ex_list', 'ez_list', 'mu_list']:
        for i in range(N):
            assert params[list_name][i].imag >= 0

    kz_list = [
        cmath.sqrt(w**2 * ex_list[i] * mu_list[i] / nu.c0**2 - kx**2 * ex_list[i] / ez_list[i])
        for i in range(N)
    ]
    kz_list = [(-kz if kz.imag < 0 else kz) for kz in kz_list]

    new_params = deepcopy(params)
    new_params['kz_list'] = kz_list
    return new_params

def find_kzs(params):
    """
    "params" is a dictionary containing w (angular frequency), kx (angular
    wavenumber), ex_list (unitless permittivity of each layer in x-direction),
    ez_list (ditto in z direction), mu_list (unitless permeability in
    y-direction).
    
    This function returns a new dictionary containing all those data PLUS
    kz_list, a list of kz in each layer.
    """
    w = params['w'] # angular frequency (w looks like omega)
    kx = params['kx']
    ex_list = params['ex_list']
    ez_list = params['ez_list']
    mu_list = params['mu_list']
    N = len(ez_list)
    assert N == len(ex_list) == len(ez_list) == len(mu_list) >= 2
    assert w > 0
    for list_name in ['ex_list', 'ez_list', 'mu_list']:
        for i in range(N):
            assert params[list_name][i].imag >= 0

    kz_list = [cmath.sqrt(w**2 * ex_list[i] * mu_list[i] / nu.c0**2
                         - kx**2 * ex_list[i] / ez_list[i]) for i in range(N)]
    # Imaginary parts should be nonnegative
    kz_list = [(-kz if kz.imag < 0 else kz) for kz in kz_list]
    
    new_params = deepcopy(params)
    new_params['kz_list'] = kz_list
    return new_params


def _get_null_vector(mat):
    """Возвращает собственный вектор матрицы, соответствующий нулевому собственному значению."""
    eigenvals, eigenvecs = np.linalg.eig(mat)
    idx = np.argmin(np.abs(eigenvals))
    return eigenvecs[:, idx]


def _calculate_H_components(null_vector, N):
    """Возвращает списки H_up и H_down."""
    H_up_list = [0]
    H_up_list.extend(null_vector[i] for i in range(1, 2*N-2, 2))
    H_down_list = [null_vector[i] for i in range(0, 2*N-2, 2)]
    H_down_list.append(0)
    return H_up_list, H_down_list


def _calculate_Ex_components(H_up_list, H_down_list, kz_list, w, ex_list):
    """Рассчитывает Ex_up и Ex_down."""
    N = len(H_up_list)
    Ex_up_list = [H_up * kz / (w * ex * nu.eps0) for H_up, kz, ex in zip(H_up_list, kz_list, ex_list)]
    Ex_down_list = [-H_down * kz / (w * ex * nu.eps0) for H_down, kz, ex in zip(H_down_list, kz_list, ex_list)]
    return Ex_up_list, Ex_down_list


def _calculate_Ez_components(H_up_list, H_down_list, kx, w, ez_list):
    """Рассчитывает Ez_up и Ez_down."""
    N = len(H_up_list)
    Ez_up_list = [-H_up * kx / (w * ez * nu.eps0) for H_up, ez in zip(H_up_list, ez_list)]
    Ez_down_list = [-H_down * kx / (w * ez * nu.eps0) for H_down, ez in zip(H_down_list, ez_list)]
    return Ez_up_list, Ez_down_list


def _normalize_fields(Ez_up_list, scale_factor=1e9):
    """Нормирует поля так, чтобы максимальное значение Ez_up было равно 1."""
    largest_Ez_index = np.argmax(np.abs(Ez_up_list))
    norm_factor = 1 / Ez_up_list[largest_Ez_index]
    Ez_up_normed = [e * norm_factor for e in Ez_up_list]
    return norm_factor, Ez_up_normed


def _calculate_Sx_contributions(Ez_up_list, Ez_down_list, H_up_list, H_down_list, kz_list, d_list):
    """Рассчитывает вклад Sx из каждого слоя."""
    N = len(Ez_up_list)
    Sx_list = []
    for i in range(N):
        Ez_up = Ez_up_list[i]
        Ez_down = Ez_down_list[i]
        H_up_star = H_up_list[i].conjugate()
        H_down_star = H_down_list[i].conjugate()
        kz = kz_list[i]
        d = d_list[i]

        Sx = 0
        if Ez_up * H_up_star != 0:
            Sx += ((-Ez_up * H_up_star) / (4 * kz.imag) *
                   (1 - cmath.exp(-2 * kz.imag * d)))
        if Ez_down * H_down_star != 0:
            Sx += ((-Ez_down * H_down_star) / (4 * kz.imag) *
                   (1 - cmath.exp(-2 * kz.imag * d)))
        if Ez_down * H_up_star != 0:
            Sx += ((-Ez_down * H_up_star) / (4j * kz.real) *
                   (1 - cmath.exp(-2j * kz.real * d)) *
                   cmath.exp(1j * kz * d))
        if Ez_up * H_down_star != 0:
            Sx += ((-Ez_up * H_down_star) / (4j * kz.real) *
                   (1 - cmath.exp(-2j * kz.real * d)) *
                   cmath.exp(1j * kz * d))

        Sx_list.append(Sx)
    return Sx_list

def find_all_params_from_kx(params):
    """
    Рассчитывает все параметры моды при заданном kx.

    Args:
        params (dict): параметры системы, включая:
            - kx (complex): горизонтальное волновое число,
            - w (float): угловая частота,
            - d_list (list): толщины слоёв.

    Returns:
        dict: обновлённый params с полями H_up/down, Ex/up/down, Ez/up/down, Sx и др.
    """
    new_params = find_kzs(deepcopy(params))
    w = new_params['w']
    d_list = new_params['d_list']
    kx = new_params['kx']
    kz_list = new_params['kz_list']
    ex_list = new_params['ex_list']
    ez_list = new_params['ez_list']
    mu_list = new_params['mu_list']
    N = len(mu_list)

    mat = bc_matrix(new_params)
    null_vector = _get_null_vector(mat)
    
    H_up_list, H_down_list = _calculate_H_components(null_vector, N)
    Ex_up_list, Ex_down_list = _calculate_Ex_components(H_up_list, H_down_list, kz_list, w, ex_list)
    Ez_up_list, Ez_down_list = _calculate_Ez_components(H_up_list, H_down_list, kx, w, ez_list)

    # Нормировка
    _, Ez_up_normed = _normalize_fields(Ez_up_list)
    scale_factor = 1 / Ez_up_normed[np.argmax(np.abs(Ez_up_normed))]
    
    def apply_scale(lst):
        return [x * scale_factor for x in lst]

    H_up_list = apply_scale(H_up_list)
    H_down_list = apply_scale(H_down_list)
    Ex_up_list = apply_scale(Ex_up_list)
    Ex_down_list = apply_scale(Ex_down_list)
    Ez_up_list = apply_scale(Ez_up_list)
    Ez_down_list = apply_scale(Ez_down_list)

    # Расчёт вектора Пойнтинга
    Sx_list = _calculate_Sx_contributions(
        Ez_up_list, Ez_down_list, H_up_list, H_down_list, kz_list, d_list
    )
    Sx_total = sum(Sx_list)

    # Сохраняем результаты в params
    new_params.update({
        'H_up_list': H_up_list,
        'H_down_list': H_down_list,
        'Ex_up_list': Ex_up_list,
        'Ex_down_list': Ex_down_list,
        'Ez_up_list': Ez_up_list,
        'Ez_down_list': Ez_down_list,
        'Sx_list': Sx_list,
        'Sx_total': Sx_total
    })

    # Генерируем список нижних границ слоёв
    layer_bottom_list = [-inf, 0]
    for i in range(1, N - 1):
        layer_bottom_list.append(layer_bottom_list[-1] + d_list[i])

    new_params['layer_bottom_list'] = layer_bottom_list
    return new_params

# def find_all_params_from_kx(params):
#     """
#     params is a dictionary containing kx and other simulation parameters like
#     w, d_list, etc. It is assumed that this kx really is a mode!
    
#     This function calculates kz_list, H_up_list, H_down_list, Ex_up_list,
#     Ex_down_list, Ez_up_list, Ez_down_list.
#     It returns a new parameter dictionary containing all the old information
#     plus those newly-calculated parameters.
    
#     This is linear optics, so you can scale the E and H up or down by any
#     constant factor. (And Poynting vector by the square of that factor.)
#     I chose the normalization that makes the maximum of Ez_up_list equal to
#     1 V/nm. (This is arbitrary.)
    
#     layer_bottom_list[i] is the z-coordinate of the bottom of layer i. Assume
#     that layer 0 is z<0,
#     layer 1 is 0 < z < d_list[1],
#     layer 2 is d_list[1] < z < d_list[1] + d_list[2], etc.
#     """
#     new_params = find_kzs(deepcopy(params))
#     w = new_params['w']
#     d_list = new_params['d_list']
#     kx = new_params['kx']
#     kz_list = new_params['kz_list']
#     ex_list = new_params['ex_list']
#     ez_list = new_params['ez_list']
#     mu_list = new_params['mu_list']
#     N = len(mu_list)
    
#     mat = bc_matrix(new_params)
#     eigenvals, eigenvecs = np.linalg.eig(mat)
#     which_eigenval_is_zero = np.argmin(np.abs(eigenvals))
#     null_vector = eigenvecs[:,which_eigenval_is_zero]
#     if False:
#         print('null vector:')
#         print(null_vector)
#         print('matrix entry absolute values:')
#         print(np.abs(mat))
#         print('abs(mat . null_vector) should be 0:')
#         print(np.abs(np.dot(mat, null_vector)))
#         print('calculated eigenvalue:')
#         print(eigenvals[which_eigenval_is_zero])
#     H_up_list = [0]
#     H_up_list.extend(null_vector[i] for i in range(1, 2*N-2, 2))
#     H_down_list = [null_vector[i] for i in range(0, 2*N-2, 2)]
#     H_down_list.append(0)
#     assert N == len(H_up_list) == len(H_down_list)
    
#     Ex_up_list = [H_up_list[i] * kz_list[i] / (w * ex_list[i] * nu.eps0)
#                                                             for i in range(N)]
#     Ex_down_list = [-H_down_list[i] * kz_list[i] / (w * ex_list[i] * nu.eps0)
#                                                             for i in range(N)]
#     Ez_up_list = [-H_up_list[i] * kx / (w * ez_list[i] * nu.eps0)
#                                                             for i in range(N)]
#     Ez_down_list = [-H_down_list[i] * kx / (w * ez_list[i] * nu.eps0)
#                                                             for i in range(N)]
    
#     # normalize E and H.
#     largest_Ez_up_index = np.argmax(np.abs(np.array(Ez_up_list)))
#     scale_factor = (1 * nu.V/nu.nm) / Ez_up_list[largest_Ez_up_index]
#     for X_list in [H_up_list, H_down_list, Ex_up_list, Ex_down_list,
#                    Ez_up_list, Ez_down_list]:
#         for i in range(N):
#             X_list[i] *= scale_factor
#     new_params['H_up_list'] = H_up_list
#     new_params['H_down_list'] = H_down_list
#     new_params['Ex_up_list'] = Ex_up_list
#     new_params['Ex_down_list'] = Ex_down_list
#     new_params['Ez_up_list'] = Ez_up_list
#     new_params['Ez_down_list'] = Ez_down_list
    
#     # x-component of complex Poynting vector, integrated over a layer
#     Sx_list = []
#     for i in range(N):
#         Ez_up = Ez_up_list[i]
#         Ez_down = Ez_down_list[i]
#         H_up_star = H_up_list[i].conjugate()
#         H_down_star = H_down_list[i].conjugate()
#         kz = kz_list[i]
#         d = d_list[i]
#         Sx = 0
#         # add each term only if it's nonzero, to avoid 0 * nan in top and
#         # bottom layers
#         if Ez_up * H_up_star != 0:
#             Sx += ((-Ez_up * H_up_star) / (4 * kz.imag)
#                     * (1 - cmath.exp(-2 * kz.imag * d)))
#         if Ez_down * H_down_star != 0:
#             Sx += ((-Ez_down * H_down_star) / (4 * kz.imag)
#                     * (1 - cmath.exp(-2 * kz.imag * d)))
#         if Ez_down * H_up_star != 0:
#             Sx += ((-Ez_down * H_up_star) / (4j * kz.real)
#                    * (1 - cmath.exp(-2j * kz.real * d))
#                    * cmath.exp(1j * kz * d))
#         if Ez_up * H_down_star != 0:
#             Sx += ((-Ez_up * H_down_star) / (4j * kz.real)
#                    * (1 - cmath.exp(-2j * kz.real * d))
#                    * cmath.exp(1j * kz * d))
#         Sx_list.append(Sx)
#     new_params['Sx_list'] = Sx_list
#     # x-component of complex Poynting vector, integrated over all layers
#     Sx_total = sum(Sx_list)
#     new_params['Sx_total'] = Sx_total
    
#     layer_bottom_list = [-inf, 0]
#     for i in range(1,N-1):
#         layer_bottom_list.append(layer_bottom_list[-1] + d_list[i])
    
#     new_params['layer_bottom_list'] = layer_bottom_list
#     return new_params

def find_layer(z, params):
    """
    Return the layer index (0 through N-1) in which you find the z-coordinate
    z. At a layer boundary, returns either one of those two layers arbitrarily.
    """
    N = len(params['d_list'])
    for i in range(N):
        if z <= params['layer_bottom_list'][i]:
            return i-1
    return N-1