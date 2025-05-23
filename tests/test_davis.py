from __future__ import division, print_function
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from math import pi
from copy import deepcopy

from utils.config_loader import load_config
import numericalunits as nu
from physics.modes import find_all_params_from_kx
# from physics.zeros import find_kx
from physics.zeros import ModeFinder
from physics.fields import calculate_Hy
from utils.checks import check_mode

from system.optical_system import OpticalSystemConfig, OpticalState, OpticalSystem

CONFIG_PATH = Path(__file__).parent.parent / "config" / "davis.yaml"
            
def test_davis():
    """
    Воспроизводит результаты T.J. Davis (2009) по поверхностным плазмонам
    в многослойной системе.

    http://dx.doi.org/10.1016/j.optcom.2008.09.043

    Эта функция:
    - Создаёт структуру из статьи (7 слоёв: glass, Au, MgF2, Au, MgF2, Au, glass),
    - Ищет комплексные моды с помощью `find_kx()`,
    - Сравнивает их с опубликованными значениями,
    - Строит графики Hy(z) для каждой моды.

    Returns:
        None
    """
    params_dict = load_config(CONFIG_PATH)
    config = OpticalSystemConfig(**params_dict)

    params = OpticalSystem(config=config)
        
    finder = ModeFinder()
    kx_list = finder.find_kx_modes(params, show_progress=False,
                      search_domain=[-0.05/nu.nm, 0.05/nu.nm, 0, 0.4/nu.nm],
                      grid_points=20, iterations=10, reduction_factor=9,
                      plot_full_region=True)
    print('kx_list -- ' + str(len(kx_list)) + ' entries...')
    print(['(%.5g+%.5gj) rad/um' % (kx.real / nu.um**-1, kx.imag / nu.um**-1)
                                                          for kx in kx_list])
    
    # The modes discoved by Davis (Table 1 of paper)
    davis_modes = [x * nu.nm**-1 for x in [1.2969e-2 + 2.7301e-5j,
                                           1.2971e-2 + 2.7644e-5j,
                                           3.0454e-2 + 3.7872e-4j,
                                           3.2794e-2 + 4.6749e-4j,
                                          -2.1254e-4 + 5.4538e-2j,
                                           1.2634e-3 + 5.4604e-2j]]
    for i in range(len(davis_modes)):
        davis_kx = davis_modes[i]
        print('Looking for "Mode ' + str(i+1) + '" in Davis paper -- kx =',
              '(%.5g+%.5gj) rad/um' % (davis_kx.real / nu.um**-1, davis_kx.imag / nu.um**-1))
        which_kx = np.argmin(abs(np.array(kx_list) - davis_kx))
        my_kx = kx_list[which_kx]
        print('Relative error: ',
              abs(my_kx - davis_kx) / (abs(my_kx) + abs(davis_kx)))

    print('---')
    print('Are the last two modes missing? They were for me. Re-trying with a')
    print('smaller search domain (zoomed towards kx=0). (By the way, ')
    print('using a larger number for grid_points would also work here.)')
    print('---')

    kx_list2 = finder.find_kx_modes(params, show_progress=False,
                      search_domain=[-0.05/nu.nm, 0.05/nu.nm, 0, 0.1/nu.nm],
                      grid_points=20, iterations=10, reduction_factor=9,
                      plot_full_region=True)
    print('kx_list2 -- ' + str(len(kx_list2)) + ' entries...')
    print(['(%.5g+%.5gj) rad/um' % (kx.real / nu.um**-1, kx.imag / nu.um**-1)
                                                          for kx in kx_list2])
    
    for i in range(len(davis_modes)):
        davis_kx = davis_modes[i]
        print('Looking for "Mode ' + str(i+1) + '" in Davis paper -- kx =',
              '(%.5g+%.5gj) rad/um' % (davis_kx.real / nu.um**-1, davis_kx.imag / nu.um**-1))
        which_kx = np.argmin(abs(np.array(kx_list2) - davis_kx))
        my_kx = kx_list2[which_kx]
        print('Relative error: ',
              abs(my_kx - davis_kx) / (abs(my_kx) + abs(davis_kx)))
        
        new_params = deepcopy(params)
        new_params['kx'] = my_kx
        out = find_all_params_from_kx(new_params)
        print(out)
        plt.figure()
        plt.title('"Mode ' + str(i+1) + '" in Davis paper -- Plot of Re(Hy) and Im(Hy)')
        zs = np.linspace(-300 * nu.nm, 500 * nu.nm, num=400)
        Hs = np.array([calculate_Hy(z, out) for z in zs])
        plt.plot(zs / nu.nm, Hs.real / max(abs(Hs)),
                 zs / nu.nm, Hs.imag / max(abs(Hs)))
        plt.xlabel('z (nm)')
        plt.ylabel('Hy (arbitrary units)')
        plt.show()