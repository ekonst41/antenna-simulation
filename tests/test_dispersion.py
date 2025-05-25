import yaml
from runner import run
import numericalunits as nu
import numpy as np
from visualization.visualize_dispersion import plot_kx_dispersion

wavelengths = np.linspace(400, 800, num=21)
config = yaml.load(open("config/air_metall_air.yaml"), Loader=yaml.Loader)
search_domain_kx = [-0.05/nu.nm, 0.05/nu.nm, 0, 0.4/nu.nm]

for w in wavelengths:
    config["system"]["w"] = f"{w} nm"
    with open("config/config_tmp.yaml", "w") as f:
        yaml.dump(config, f, sort_keys=False)
    print(f"\n=== Running for λ = {w} nm ===")
    run("config_tmp.yaml", search_domain_kx=search_domain_kx, visual=False, save_kx=True, save_kx_path='results/dispersion/data')

plot_kx_dispersion('results/dispersion/data')
