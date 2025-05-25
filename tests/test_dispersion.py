import yaml
from runner import run
import numericalunits as nu
import numpy as np
from visualization.visualize_dispersion import plot_kx_dispersion

wavelengths = np.linspace(400, 800, num=21)
w_plasma = 1.39 * 10**16
gamma = 91.55 * 10**12

config = yaml.load(open("config/air_metall_air.yaml"), Loader=yaml.Loader)
search_domain_kx = [-0.05/nu.nm, 0.05/nu.nm, 0, 0.4/nu.nm]
nu.reset_units('SI')

for w in wavelengths:
    config["system"]["w"] = f"{w} nm"
    for layer in config["system"]["layers"]:
        if layer.get("name") == "Au":
            cycle_w = 2 * np.pi * nu.c0 / (w * 10**(-9))
            eps = complex(round(1 - (w_plasma  / cycle_w)**2, 4),
                          round(gamma * w_plasma**2 / (cycle_w**3), 4))
            layer["ex"] = eps
            layer["ez"] = eps
    with open("config/config_tmp.yaml", "w") as f:
        yaml.dump(config, f, sort_keys=False)
    print(f"\n=== Running for λ = {w} nm ===")
    run("config_tmp.yaml", visual=False, save_kx=True, save_kx_path='results/dispersion/data')

plot_kx_dispersion('results/dispersion/data', wavelengths)
