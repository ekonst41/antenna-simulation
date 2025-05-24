import numpy as np
import matplotlib.pyplot as plt
import numericalunits as nu
from physics.fields import calculate_Hy

def visualize(out, params, left_bound=-300, right_bound=500, num=400):
    plt.style.use('seaborn-v0_8-pastel')
    
    plt.figure(figsize=(10, 6), dpi=100)
    
    colors = {
        're': '#88CCEE',
        'im': '#CC6677',
        'bg': '#F0F0F0',
        'grid': '#DDDDDD',
        'boundary': '#332288'
    }
    
    ax = plt.gca()
    ax.set_facecolor(colors['bg'])
    
    kx = params['kx']
    zs = np.linspace(left_bound * nu.nm, right_bound * nu.nm, num=num)
    Hs = np.array([calculate_Hy(z, out) for z in zs])

    Hs_real = Hs.real / max(abs(Hs))
    Hs_imag = Hs.imag / max(abs(Hs))
    
    plt.plot(zs / nu.nm, Hs_real, 
             label=r"$Re(H_y)$", 
             color=colors['re'], 
             linewidth=2.5,
             alpha=0.8)
    
    plt.plot(zs / nu.nm, Hs_imag, 
             label=r"$Im(H_y)$", 
             color=colors['im'], 
             linewidth=2.5,
             alpha=0.8)
    
    boundary_pos = sum(params['d_list'][1:-1]) / nu.nm
    plt.axvline(x=0, 
                color=colors['boundary'], 
                linestyle='--', 
                linewidth=1.5,
                alpha=0.7,
                label="Границы структуры")
    
    plt.axvline(x=boundary_pos, 
                color=colors['boundary'], 
                linestyle='--', 
                linewidth=1.5,
                alpha=0.7)
    
    title = r"Распределение поля $H_y$ для моды $k_x = {:.4f} + {:.4f}i$".format(
        kx.real / nu.um**-1, 
        kx.imag / nu.um**-1
    )
    plt.title(title, pad=20, fontsize=14)
    
    plt.xlabel(r"$z$, нм", fontsize=12)
    plt.ylabel(r"Нормированная амплитуда $H_y$", fontsize=12)

    legend = plt.legend(frameon=True, fontsize=11)
    legend.get_frame().set_facecolor('white')
    legend.get_frame().set_alpha(0.8)
    
    plt.grid(True, color=colors['grid'], linestyle='-', alpha=0.6)
    
    plt.xlim(left_bound, right_bound)
    
    plt.tight_layout()
    
    plt.show()