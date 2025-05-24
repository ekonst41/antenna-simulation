import numpy as np
import matplotlib.pyplot as plt
import numericalunits as nu

from physics.fields import calculate_Hy

left_bound = -300
right_bound = 500
num = 400

def visualize(out, params):
    plt.figure()
    kx = params['kx']
    plt.title(f'Re(Hy) and Im(Hy) for mode kx = {round(kx.real / nu.um**-1,4)} + {round(kx.imag / nu.um**-1,4)}i')
    zs = np.linspace(left_bound * nu.nm, right_bound * nu.nm, num=num)
    Hs = np.array([calculate_Hy(z, out) for z in zs])
    plt.plot(zs / nu.nm, Hs.real / max(abs(Hs)), label='Re(H)', color='blue')
    plt.plot(zs / nu.nm, Hs.imag / max(abs(Hs)), label='Im(H)', color='red')
    plt.axvline(x=0, color='r', linestyle='--')
    plt.axvline(x=sum(params['d_list'][1:-1]) / nu.nm, color='r', linestyle='--')
    plt.xlabel('z (nm)')
    plt.ylabel('Hy (arbitrary units)')
    plt.legend()
    plt.show()