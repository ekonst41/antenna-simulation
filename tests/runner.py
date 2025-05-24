from __future__ import division, print_function
from pathlib import Path
from copy import deepcopy
from typing import List, Optional

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
        num_to_visualize: int = 10**18
        ):

    """
    Эта функция:
    - Создаёт структуру,
    - Ищет комплексные моды с помощью `find_kx()`,
    - Сравнивает их с опубликованными значениями,
    - Строит графики Hy(z) для каждой моды.

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
                                   plot_full_region=True)
    print('kx_list -- ' + str(len(kx_list)) + ' items')
    print('---')
    for kx in kx_list:
        print(f'{round(kx.real / nu.um**-1, 3)} + {round(kx.imag / nu.um**-1, 3)} i')
    print('---')


    for i in range(min(num_to_visualize, len(kx_list))):
        new_params = deepcopy(params)
        new_params['kx'] = kx_list[i]
        out = find_all_params_from_kx(new_params)
        visualize(out, new_params)
