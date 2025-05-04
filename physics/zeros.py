from tmm.tmm import bc_matrix
from physics.modes import find_kzs

import numpy as np
import cmath
import math
from numpy import pi
import scipy.optimize as opt
import matplotlib.pyplot as plt
import numericalunits as nu

inf = float('inf')

def find_all_zeros(min_re, max_re, min_im, max_im, fn,
                   grid_points, iterations, reduction_factor,
                   plot_full_region, show_progress):
    """
    Ищет все нули переданной комплексной функции fn(z) в заданной области комплексной плоскости.

    Args:
        min_re (float): Минимальное значение действительной части z.
        max_re (float): Максимальное значение действительной части z.
        min_im (float): Минимальное значение мнимой части z.
        max_im (float): Максимальное значение мнимой части z.
        fn (function): Функция от одного комплексного аргумента.
        grid_points (int, optional): Число точек по каждой оси для начального поиска. По умолчанию 20.
        iterations (int, optional): Число итераций сужения области поиска. По умолчанию 9.
        reduction_factor (int, optional): Во сколько раз уменьшать область поиска на каждой итерации. По умолчанию 9.
        plot_full_region (bool, optional): Если True, строит графики |fn(z)| и интеграла вокруг прямоугольника. По умолчанию True.
        show_progress (bool, optional): Если True, выводит промежуточные сообщения о прогрессе. По умолчанию False.

    Returns:
        list of complex: Список найденных кандидатов на нули функции fn(z).
    """
    # Check arguments
    assert reduction_factor > 1 and max_re > min_re and max_im > min_im
    assert (max_re.imag == 0 and min_re.imag == 0
            and max_im.imag == 0 and min_im.imag == 0)
    # Edge-point rejection (see below) relies on the following assumption:
    assert grid_points > 2 * reduction_factor
    
    
    if plot_full_region:
    
        def inverse_fn(z):
            """ 1 / fn(z) """
            f = fn(z)
            return inf if f == 0 else 1/f
            
        def contour_int(z, d_re, d_im):
            """
            Approximate the contour integral of inverse_fn around a point z,
            using a rectangle of half-width d_re (in real direction) and
            half-height d_im. Just a nice plot that makes zeros stand out.
            """
            assert d_re.imag == 0 and d_im.imag == 0 and d_re > 0 and d_im > 0
            below = inverse_fn(z - 1j * d_im)
            above = inverse_fn(z + 1j * d_im)
            left = inverse_fn(z - d_re)
            right = inverse_fn(z + d_re)
            return (below * (2 * d_re) + right * (2j * d_im)
                    + above * (-2 * d_re) + left * (-2j * d_im))
        
        res, re_step = np.linspace(min_re, max_re, num=100, retstep=True)
        ims, im_step = np.linspace(min_im, max_im, num=100, retstep=True)
        
        fig = plt.figure()
        direct_plot = fig.add_subplot(111)
        data = [[math.log10(abs(fn(re + 1j * im))) for re in res] for im in ims]
        direct_plot.imshow(data, extent=(min_re * nu.um, max_re * nu.um,
                                         min_im * nu.um, max_im * nu.um),
                           origin='lower')
        direct_plot.set_xlabel('Re(kx) [rad/um]')
        direct_plot.set_ylabel('Im(kx) [rad/um]')
        direct_plot.set_title('log(|fn(z)|) -- Looking for minima (blue)')

        fig = plt.figure()
        contour_plot = fig.add_subplot(111)
        data = [[-math.log10(abs(contour_int(re + 1j * im, re_step, im_step)))
                                                  for re in res] for im in ims]
        contour_plot.imshow(data, extent=(min_re * nu.um, max_re * nu.um,
                                         min_im * nu.um, max_im * nu.um),
                           origin='lower')
        contour_plot.set_xlabel('Re(kx) [rad/um]')
        contour_plot.set_ylabel('Im(kx) [rad/um]')
        contour_plot.set_title(
             '-log(|contour integral of 1/fn(z) around a little rectangle|)\n'
           + ' -- This plot highlights zeros in fn(z), but also lines of\n'
           + 'discontinuity (where top or bottom kz is pure-imaginary)')
    
    # "regions" is a list where each entry has the form
    # [min_re, max_re, min_im, max_im]. Each entry describes a region in which we
    # are seeking local minima.
    regions = [[min_re, max_re, min_im, max_im]]
    
    region_width_re = max_re - min_re
    region_width_im = max_im - min_im
    
    for iteration_number in range(iterations):
        # all_local_mins will be a list of (x, y) for every local minimum in
        # every region. This is used to generate the next iteration.
        all_local_mins = []
        for region_index in range(len(regions)):
            min_re_now, max_re_now, min_im_now, max_im_now = regions[region_index]
            results_grid = []
            re_list, re_step = np.linspace(min_re_now, max_re_now, num=grid_points, retstep=True)
            im_list, im_step = np.linspace(min_im_now, max_im_now, num=grid_points, retstep=True)
            fn_to_minimize = lambda z : abs(fn(z))
            
            results_grid = [[fn_to_minimize(re + 1j * im) for im in im_list]
                                                             for re in re_list]
            results_grid = np.array(results_grid)
            # local_mins will be a list of (i,j) where (re_list[i], im_list[j])
            # is a local minimum on the results_grid
            local_mins = []
            for i in range(grid_points):
                for j in range(grid_points):
                    is_min = all(results_grid[i2, j2] >= results_grid[i,j]
                                    for i2 in [i-1, i, i+1]
                                      for j2 in [j-1, j, j+1]
                                        if (0 <= i2 < grid_points
                                             and 0 <= j2 < grid_points))
                    if is_min:
                        local_mins.append((i,j))
            # local_mins_OK is the subset of local_mins that passes the
            # the edge-rejection test.
            # The edge-rejection test says that after the 0'th iteration, any
            # point at an edge is probably not a true minimum.
            
            local_mins_OK = []
            for (i,j) in local_mins:
                z_now = re_list[i] + 1j * im_list[j]
                if iteration_number >= 2 and (i == 0 or j == 0 or
                                    i == grid_points-1 or j == grid_points-1):
                    # Rejecting an edge point...
                    if show_progress:
                        print('----')
                        print('Deleting edge point: region #'
                              + str(region_index+1) + '  (i,j)=', (i,j),
                              '  kx in rad/um=',
                              z_now / nu.um**-1,
                              '  fn(z)=', fn(z_now))
                else:
                    local_mins_OK.append((i,j))
            
            # Add local_mins_OK entries into all_local_mins
            for (i,j) in local_mins_OK:
                all_local_mins.append(re_list[i] + 1j * im_list[j])
            
            if show_progress:
                print('----')
                print('iter #' + str(iteration_number)
                    + ' , region #' + str(region_index+1) + ' of ' + str(len(regions))
                    + ' , ' + str(len(local_mins_OK)) + ' minima')
                if len(local_mins_OK) > 0:
                    print('For each, here is ((i, j), kx in rad/um, fn(kx)):')
                    print([((i, j), (re_list[i] + 1j * im_list[j]) / nu.um**-1,
                                                  fn(re_list[i] + 1j * im_list[j]))
                                                      for (i,j) in local_mins_OK])

        # Now we've gone through every region.
        # Delete redundant minima that showed up in overlapping regions.
        all_local_mins_norepeat = []
        def is_repeat(z1, z2):
            return ((abs((z1 - z2).real) <= 0.5 * re_step) and
                    (abs((z1 - z2).imag) <= 0.5 * im_step))
        for z_now in all_local_mins:
            if not any(is_repeat(z_now, z) for z in all_local_mins_norepeat):
                all_local_mins_norepeat.append(z_now)
        if show_progress:
            num_deleted = len(all_local_mins) - len(all_local_mins_norepeat)
            if num_deleted > 0:
                print('----')
                print('After iter #' + str(iteration_number)
                    + ', deleted ' + str(num_deleted) + ' redundant point(s)')

        all_local_mins = all_local_mins_norepeat
        
        if show_progress:
            print('----')
            print('** After iter #' + str(iteration_number) + ', we have '
                  + str(len(all_local_mins)) + ' candidate minima')
        
        region_width_re /= reduction_factor
        region_width_im /= reduction_factor
        
        regions = [[z.real - region_width_re / 2, z.real + region_width_re / 2,
                    z.imag - region_width_im / 2, z.imag + region_width_im / 2]
                        for z in all_local_mins]
    
    # Done with main algorithm. Show the discovered minima on the plots as
    # white X's. Note: Zeros outside the plot region will not be seen here,
    # but the function still returns them.
    if plot_full_region:
        # Keep the image filling the plot area
        direct_plot.autoscale(False)
        contour_plot.autoscale(False)
        for z in all_local_mins:
            direct_plot.plot(z.real * nu.um, z.imag * nu.um, 'wx')
            contour_plot.plot(z.real * nu.um, z.imag * nu.um, 'wx')
    return all_local_mins

def find_kx(input_params, search_domain=None, show_progress=False,
            grid_points=20, iterations=9, reduction_factor=9,
            plot_full_region=True):
    """
    Ищет возможные значения комплексного волнового числа kx для поверхностных плазмонов.

    Args:
        input_params (dict): Параметры модели: частота w, толщины d_list, диэлектрические проницаемости ex_list, ez_list, и т.п.
        search_domain (list, optional): Область поиска в виде [min_re, max_re, min_im, max_im]. По умолчанию None.
        show_progress (bool, optional): Если True, выводит информацию о ходе поиска. По умолчанию False.
        grid_points (int, optional): Число точек на сетке. По умолчанию 20.
        iterations (int, optional): Число итераций. По умолчанию 9.
        reduction_factor (int, optional): Коэффициент уменьшения области поиска. По умолчанию 9.
        plot_full_region (bool, optional): Если True, строит графики. По умолчанию True.

    Returns:
        list of complex: Список значений kx, соответствующих возможным поверхностным плазмонам.
    """
    w = input_params['w']
    d_list = input_params['d_list']
    ex_list = input_params['ex_list']
    ez_list = input_params['ez_list']
    mu_list = input_params['mu_list']
    N = len(mu_list)
    assert N == len(d_list) == len(ex_list) == len(ez_list)
    # error(z) approaches 0 as kx = z approaches a true plasmon mode.
    # It's proportional to the determinant of the boundary-condition matrix, 
    # which equals zero at modes.
    def error(kx):
        if kx == 0:
            return inf
        temp_params = input_params.copy()
        temp_params['kx'] = kx
        should_be_zero = np.linalg.det(bc_matrix(find_kzs(temp_params)))
        return should_be_zero / kx**(N+1)
        # "return should_be_zero" is also OK but has an overall slope that
        # makes it harder to find zeros; also, there's a false-positive at k=0.
    
    # choose the region in which to search for minima. My heuristic is:
    # The upper limit of kx should be large enough that
    # 2 * pi * i * kzm * d ~ 20 for the thinnest layer we have, or 3 times
    # the light-line, whichever is bigger.
    if search_domain is None:
        kx_re_max = max(max(abs((20 / (2 * pi * d_list[i]))
                        * cmath.sqrt(ez_list[i] / ex_list[i])) for i in range(1,N)),
                    3 * w / nu.c0)
        kx_re_min = -kx_re_max
        kx_im_min = 0
        kx_im_max = abs(kx_re_max)
    else:
        kx_re_min = search_domain[0]
        kx_re_max = search_domain[1]
        kx_im_min = search_domain[2]
        kx_im_max = search_domain[3]
    
    # Main part of function: Call find_all_zeros()
    kx_list = find_all_zeros(kx_re_min, kx_re_max, kx_im_min, kx_im_max, error,
                           show_progress=show_progress, grid_points=grid_points,
                           iterations=iterations,
                           reduction_factor=reduction_factor,
                           plot_full_region=plot_full_region)
    
    # sort and remove "repeats" with opposite signs
    kx_list = sorted(kx_list, key=(lambda kx : abs(kx)))
    i=0
    while i < len(kx_list) - 1:
        if abs(kx_list[i] + kx_list[i+1]) <= 1e-6 * (abs(kx_list[i]) + abs(kx_list[i+1])):
            kx_list.pop(i)
        else:
            i += 1
    
    # Fix amplifying waves
    kx_list = [(-kx if (kx.imag < 0 or (kx.imag==0 and kx.real < 0)) else kx)
                                                            for kx in kx_list]
    
    return kx_list