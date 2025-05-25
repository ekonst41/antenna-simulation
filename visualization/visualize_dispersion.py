import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

def plot_kx_dispersion(folder: str):
    folder_path = Path(folder)
    npz_files = sorted(folder_path.glob("wavelength_*nm.npz"),
                       key=lambda p: int(p.stem.split('_')[1].replace('nm', '')))

    wavelengths = []
    re_kx = []
    im_kx = []

    for file in npz_files:
        data = np.load(file)
        wl = data["wavelength"].item()
        kx = data["kx"][0]
        wavelengths.append(wl)
        re_kx.append(kx.real)
        im_kx.append(kx.imag)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))

    ax1.plot(wavelengths, re_kx, marker='o')
    ax1.set_xlabel("Длина волны, nm")
    ax1.set_ylabel("Re($k_x$)")
    ax1.set_title("Действительная часть")

    ax2.plot(wavelengths, im_kx, marker='o')
    ax2.set_xlabel("Длина волны, nm")
    ax2.set_ylabel("Im($k_x$)")
    ax2.set_title("Мнимая часть")

    plt.tight_layout()
    plt.show()
