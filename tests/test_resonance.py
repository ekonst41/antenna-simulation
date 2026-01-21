from copy import deepcopy
import numpy as np
from matplotlib import pyplot as plt
import numericalunits as nu

from physics.fields import calculate_Sx
from physics.modes import find_all_params_from_kx
from tests.runner import run
from visualization.visualize import visualize


def make_complex_grid(start: complex, finish: complex, N_real, N_imag):
    re_vals = (start.real + finish.real) / 2 if N_real == 1 else np.linspace(start.real, finish.real, num=N_real)
    im_vals = (start.imag + finish.imag) / 2 if N_imag == 1 else np.linspace(start.imag, finish.imag, num=N_imag)
    grid_2d = np.meshgrid(re_vals, im_vals)
    complex_grid = grid_2d[0] + 1j * grid_2d[1]
    return complex_grid.ravel()

def visualize_resonance(h_list, k_list, k_orig):
    k_list = [kx.real * nu.um for kx in k_list]
    h_list = [h.real for h in h_list]
    plt.figure(figsize=(10, 6))
    plt.plot(k_list, h_list,
         linestyle='-',
         color='blue',
         label='h(k)')

    plt.xlabel('Волновой вектор (k)')
    plt.ylabel('Параметр (h)')
    plt.title('Зависимость h от k')
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.xlim(left=k_list[0], right=k_list[-1])
    plt.axvline(x=k_orig.real * nu.um,
                color='r',
                linestyle='--',
                linewidth=1.5,
                alpha=0.4)
    plt.legend()
    plt.show()


def run_resonance_test(path):
    kx_list, params = run(path=path, num_to_visualize=0, visualisation_type=['H', 'S'], return_params=True)
    kx_orig = kx_list[2]
    print(kx_orig * nu.um)
    kx_list = make_complex_grid(kx_orig * 3/4, kx_orig * 5/4, 11, 1)
    #print(kx_list)
    H_list = []
    for kx in kx_list:
        new_params = deepcopy(params)
        new_params['kx'] = kx
        out = find_all_params_from_kx(new_params)
        start_of_diel = new_params['d_list'][1]
        H = calculate_Sx(start_of_diel, out)
        #print(H)
        H_list.append(H)
        visualize(out, new_params, type=['S'])
    visualize_resonance(H_list, kx_list, kx_orig)



