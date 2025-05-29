from matplotlib import pyplot as plt

from physics.fields import calculate_Sx, calculate_Sz
from physics.modes import find_kzs, find_all_params_from_kx
from visualization.visualize import visualize
from tests.runner import run
from copy import deepcopy
import numpy as np
import numericalunits as nu


def make_complex_grid(start: complex, finish: complex, N_real, N_imag):
    re_vals = (start.real + finish.real) / 2 if N_real == 1 else np.linspace(start.real, finish.real, num=N_real)
    im_vals = (start.imag + finish.imag) / 2 if N_imag == 1 else np.linspace(start.imag, finish.imag, num=N_imag)
    grid_2d = np.meshgrid(re_vals, im_vals)
    complex_grid = grid_2d[0] + 1j * grid_2d[1]
    return complex_grid.ravel()

def find_r_and_t(kz_list, ez_list):
    r = []
    t = []
    for i in range(len(ez_list) - 1):
        ei, ki = ez_list[i], kz_list[i]
        ej, kj = ez_list[i + 1], kz_list[i + 1]
        ni, nj = np.sqrt(ei), np.sqrt(ej)
        r.append((ej * ki - ei * kj) / (ei * ki + ej * kj))
        t.append(2 * ni * nj * ki / (ej * ki + ei * kj))
    return r, t

def visualize_r_and_t(kx, r, t):
    kx = np.array(kx)
    r = np.array(r)
    t = np.array(t)

    # Значения для осей
    re_kx = kx.real * nu.nm
    im_kx = kx.imag * nu.nm

    # Создаем сетку из 2x2 графиков
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # 1) r vs Re(kx)
    ax = axes[0, 0]
    ax.plot(re_kx, r.real, label='Re(r)')
    ax.plot(re_kx, r.imag, label='Im(r)')
    ax.set_xlabel('Re(kx)')
    ax.set_ylabel('r')
    ax.set_title('r vs Re(kx)')
    ax.legend()
    ax.grid(True)

    # 2) r vs Im(kx)
    ax = axes[0, 1]
    ax.plot(im_kx, r.real, label='Re(r)')
    ax.plot(im_kx, r.imag, label='Im(r)')
    ax.set_xlabel('Im(kx)')
    ax.set_ylabel('r')
    ax.set_title('r vs Im(kx)')
    ax.legend()
    ax.grid(True)

    # 3) t vs Re(kx)
    ax = axes[1, 0]
    ax.plot(re_kx, t.real, label='Re(t)')
    ax.plot(re_kx, t.imag, label='Im(t)')
    ax.set_xlabel('Re(kx)')
    ax.set_ylabel('t')
    ax.set_title('t vs Re(kx)')
    ax.legend()
    ax.grid(True)

    # 4) t vs Im(kx)
    ax = axes[1, 1]
    ax.plot(im_kx, t.real, label='Re(t)')
    ax.plot(im_kx, t.imag, label='Im(t)')
    ax.set_xlabel('Im(kx)')
    ax.set_ylabel('t')
    ax.set_title('t vs Im(kx)')
    ax.legend()
    ax.grid(True)

    plt.tight_layout()
    plt.show()

def run_reflection_test(path):
    kx_list, params  = run(path, visual=False, return_params=True)
    kx_orig = kx_list[0]
    kx_list = make_complex_grid(kx_orig * 1/2, kx_orig * 3/2, 5, 1)
    r_list = []
    tx_list = []
    tz_list = []
    for kx in kx_list:
        new_params = deepcopy(params)
        new_params['kx'] = kx
        out = find_all_params_from_kx(new_params)
        boundary_pos = sum(params['d_list'][1:-1])
        delta = boundary_pos * 10**(-4)
        Sx_left = calculate_Sx(-delta, out)
        Sz_left = calculate_Sz(-delta, out)
        Sx_right = calculate_Sx(boundary_pos + delta, out)
        Sz_right = calculate_Sz(boundary_pos + delta, out)
        #print(Sx_right.real / Sx_left.real, Sz_right.real / Sz_left.real)
        visualize(out, new_params, type=['Sz', 'S', 'H'])
        tx_list.append(Sx_right / Sx_left)
        tz_list.append(Sz_right / Sz_left)
    visualize_r_and_t(kx_list, tx_list, tz_list)
    '''kz_list = new_params['kz_list']
        ez_list = new_params['ez_list']
        r, t = find_r_and_t(kz_list, ez_list)
        r_list.append(r[1])
        t_list.append(t[1])
    visualize(kx_list, r_list, t_list)'''





