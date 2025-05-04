from __future__ import division, print_function
import numpy as np
import cmath
import numericalunits as nu

inf = float('inf')

def bc_matrix(params):
    """
    Calculate the "boundary condition matrix". This is a matrix M such that
    
    M * [[H0down],[H1up],[H1down],...] = [[0],[0],...]
    
    IF the boundary conditions are all satisfied. (See online docs for
    definitions and what's going on.)
    
    params should contain ex_list, ez_list, kx, kz_list, d_list (thickness of
    each layer, first and last should be inf.)
    """
    w = params['w']
    kx = params['kx']
    d_list = params['d_list']
    ex_list = params['ex_list']
    ez_list = params['ez_list']
    kz_list = params['kz_list']
    N = len(d_list)
    assert N == len(d_list) == len(ex_list) == len(ez_list) == len(kz_list)
    assert N >= 2
    assert d_list[0] == d_list[-1] == inf
    
    # delta = e^{i * kz * d}, i.e. phase change across each layer
    # delta[0] and delta[-1] are undefined and are not used.
    delta_list = [cmath.exp(1j * kz_list[i] * d_list[i]) for i in range(N)]
    
    Ex_up_over_H_up_list = [kz_list[i] / (w * ex_list[i] * nu.eps0)
                                                           for i in range(N)]
    Ex_down_over_H_down_list = [-a for a in Ex_up_over_H_up_list]
    Ez_up_over_H_up_list = [-kx / (w * ez_list[i] * nu.eps0) for i in range(N)]
    Ez_down_over_H_down_list = Ez_up_over_H_up_list[:]
    
    mat = np.zeros((2*N-2, 2*N-2), dtype=complex)
    
    for row_now in range(N-1):
        # This row concerns continuity of Ex across the boundary between
        # layer_under and layer_over (under and over the boundary respectively)
        layer_under = row_now
        layer_over = layer_under + 1
        # up_under_index is the column index in mat that gets multiplied by
        # H_{up} in layer_under.
        up_under_index = 2 * layer_under - 1
        down_under_index = 2 * layer_under
        up_over_index = 2 * layer_over - 1
        down_over_index = 2 * layer_over
        
        if layer_under != 0:
            assert 0 <= up_under_index < 2*N-2
            mat[row_now, up_under_index] = (
                  Ex_up_over_H_up_list[layer_under] * delta_list[layer_under])
        mat[row_now, down_under_index] = Ex_down_over_H_down_list[layer_under]
        mat[row_now, up_over_index] = -Ex_up_over_H_up_list[layer_over]
        if layer_over != N-1:
            assert 0 <= down_over_index < 2*N-2
            mat[row_now, down_over_index] = (
                -Ex_down_over_H_down_list[layer_over] * delta_list[layer_over])

    for row_now in range(N-1, 2*N-2):
        # This row concerns continuity of eps_z * Ez across the boundary between
        # layer_under and layer_over (under and over the boundary respectively)
        layer_under = row_now - (N-1)
        layer_over = layer_under + 1
        # up_under_index is the column index in mat that gets multiplied by
        # H_{up} in layer_under.
        up_under_index = 2 * layer_under - 1
        down_under_index = 2 * layer_under
        up_over_index = 2 * layer_over - 1
        down_over_index = 2 * layer_over
        
        if layer_under != 0:
            assert 0 <= up_under_index < 2*N-2
            mat[row_now, up_under_index] = (ez_list[layer_under] *
                   Ez_up_over_H_up_list[layer_under] * delta_list[layer_under])
        mat[row_now, down_under_index] = (ez_list[layer_under] *
                                         Ez_down_over_H_down_list[layer_under])
        mat[row_now, up_over_index] = (-ez_list[layer_over] * 
                                              Ez_up_over_H_up_list[layer_over])
        if layer_over != N-1:
            assert 0 <= down_over_index < 2*N-2
            mat[row_now, down_over_index] = (-ez_list[layer_over] *
                 Ez_down_over_H_down_list[layer_over] * delta_list[layer_over])
    
    return mat