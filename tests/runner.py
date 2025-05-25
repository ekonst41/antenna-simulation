from __future__ import division, print_function
from pathlib import Path
from copy import deepcopy
from typing import List, Optional
import numpy as np

from utils.config_loader import load_config
import numericalunits as nu
from physics.modes import find_all_params_from_kx
from physics.zeros import ModeFinder
from system.optical_system import OpticalSystemConfig, OpticalState, OpticalSystem

from visualization.visualize import visualize

def run(path: str,
        search_domain_kx: Optional[List[float]] = None,
        grid_points: int = 20,
        iterations: int = 10,
        reduction_factor: int = 9,
        visual = True,
        num_to_visualize: int = 10**18,
        visualisation_type: List = ['H'],
        save_kx: bool = False,
        save_kx_path = None
        ):

    """
    Эта функция:
    - Создаёт структуру,
    - Ищет комплексные моды с помощью `find_kx()`,
    - Строит графики.

    Returns:
        None
    """

    CONFIG_PATH = Path(__file__).parent.parent / "config" / path

    params_dict = load_config(CONFIG_PATH)
    config = OpticalSystemConfig(**params_dict)

    params = OpticalSystem(config=config)

    print(f"Created optical system with {len(params['layers'])} layers")

    finder = ModeFinder()
    kx_list = finder.find_kx_modes(params, show_progress=False,
                                   search_domain=search_domain_kx,
                                   grid_points=grid_points, iterations=iterations, reduction_factor=reduction_factor,
                                   plot_full_region=False)
    print('kx_list -- ' + str(len(kx_list)) + ' items')
    print('---')
    for kx in kx_list[:5]:
        print(f'{round(kx.real / nu.um**-1, 3)} + {round(kx.imag / nu.um**-1, 3)} i')
    print('...')
    print('---')

    if save_kx:
        wavelength = round(2 * np.pi * nu.c0 / (params['w'] * nu.nm))
        save_kx_array_to_file(kx_list, wavelength, save_kx_path)

    for i in range(min(0 if not visual else num_to_visualize, len(kx_list))):
        new_params = deepcopy(params)
        new_params['kx'] = kx_list[i]
        out = find_all_params_from_kx(new_params)
        visualize(out, new_params, type=visualisation_type)


def save_kx_array_to_file(array, wl, save_path):
    fname = f"{save_path}/wavelength_{wl}nm.npz"
    np.savez(fname, wavelength=wl, kx=array)
    print(f"Saved to file {fname}")

