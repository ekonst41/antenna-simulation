import yaml
import numpy as np
import numericalunits as nu
from runner import run
from visualization.visualize_dispersion import plot_kx_dispersion

wavelengths = np.linspace(400, 800, num=21)
epsilon = []
w_plasma = 1.39 * 10**16
gamma = 91.55 * 10**12
path = '1d_bragg.yaml'

def change_epsilon(wl):
    cycle_w = 2 * np.pi * nu.c0 / (wl * 10**(-9))
    eps = complex(round(1 - (w_plasma  / cycle_w)**2, 4),
                  round(gamma * w_plasma**2 / (cycle_w**3), 4))
    return eps

config = yaml.load(open(f"config/{path}"), Loader=yaml.Loader)
nu.reset_units('SI')

def

def change_wavelength(wavelengths):
    for w in wavelengths:
        config["system"]["w"] = f"{w} nm"
        for layer in config["system"]["layers"]:
            if layer.get("name") == "Au":
                eps = change_epsilon(w)
                epsilon.append(eps)
                layer["ex"] = eps
                layer["ez"] = eps
        with open("config/config_tmp.yaml", "w") as f:
            yaml.dump(config, f, sort_keys=False)
        print(f"\n=== Running for λ = {w} nm ===")
        run("config_tmp.yaml", visual=False, save_kx=True, save_kx_path='results/dispersion/data')

    plot_kx_dispersion('results/dispersion/data', wavelengths, epsilon)
