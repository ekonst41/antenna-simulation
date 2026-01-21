from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

def plot_kx_dispersion(folder: str, wavelength: float, epsilon):
    folder_path = Path(folder)
    npz_files = [
        folder_path / f"wavelength_{int(w)}nm.npz"
        for w in sorted(wavelength)
        if (folder_path / f"wavelength_{int(w)}nm.npz").exists()
    ]

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

    #print(re_kx)
    #print(im_kx)
    kx_th = [2 * np.pi / wavelengths[i] * np.sqrt(epsilon[i] / (1 + epsilon[i])) for i in range(len(wavelengths))]
    kx_th_re = [kx.real for kx in kx_th]
    kx_th_im = [kx.imag for kx in kx_th]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))

    ax1.plot(wavelengths, kx_th_re, color='lightblue', zorder=1, label='Теор зависимость')
    ax1.scatter(wavelengths, re_kx, marker='o', color='blue', zorder=2, label='Результаты')
    ax1.set_xlabel("Длина волны, nm")
    ax1.set_ylabel("Re($k_x$)")
    ax1.set_title("Действительная часть")

    ax2.plot(wavelengths, kx_th_im, color='lightblue', zorder=1, label='Теор зависимость')
    ax2.scatter(wavelengths, im_kx, marker='o', color='blue', zorder=2, label='Результаты')
    ax2.set_xlabel("Длина волны, nm")
    ax2.set_ylabel("Im($k_x$)")
    ax2.set_title("Мнимая часть")

    plt.tight_layout()
    plt.legend()
    plt.show()
