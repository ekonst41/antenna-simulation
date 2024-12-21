import numpy as np
import matplotlib.pyplot as plt
import imageio
import os
from tqdm import tqdm

# Параметры моделирования
grid_size = (200, 200)  # Размер сетки
dx = 0.01  # Шаг сетки (м)
dy = 0.01  # Шаг сетки (м)
dt = dx / (2 * 3e8)
total_time_steps = 500

# TE волна
Ex = np.zeros(grid_size)
Ey = np.zeros(grid_size)
Hz = np.zeros(grid_size)

epsilon_0 = 8.854e-12
mu_0 = 4e-7 * np.pi

for t in tqdm(range(total_time_steps)):
    for i in range(1, grid_size[0]):
        for j in range(1, grid_size[1]):
            Hz[i, j] += (dt / (mu_0 * dx)) * (Ey[i, j] - Ey[i-1, j]) - \
                        (dt / (mu_0 * dy)) * (Ex[i, j] - Ex[i, j-1])

    # Обновление электрического поля (Ex и Ey)
    for i in range(grid_size[0] - 1):
        for j in range(grid_size[1] - 1):
            Ex[i, j] += (dt / (epsilon_0 * dy)) * (Hz[i, j+1] - Hz[i, j])
            Ey[i, j] -= (dt / (epsilon_0 * dx)) * (Hz[i+1, j] - Hz[i, j])

    # Источник волны (плоская волна)
    Ex[100, 100] += np.sin(2 * np.pi * t / 100)

    if t % 10 == 0 and t != 0:
        plt.figure(figsize=(6, 6))
        plt.imshow(Ex, cmap='hot', extent=[0, grid_size[0]*dx, 0, grid_size[1]*dy])
        plt.colorbar(label='Ex (V/m)')
        plt.title(f"Временной шаг {t}")
        plt.xlabel("X (м)")
        plt.ylabel("Y (м)")
        plt.axis("equal")
        plt.tight_layout()
        plt.show()

