import numpy as np
import matplotlib.pyplot as plt

Nx, Ny = 100, 100
dx, dy = 0.001, 0.001
dt = 1e-12

steps = 1000

epsilon0 = 8.854e-12
mu0 = 4 * np.pi * 1e-7
epsilon = 1
mu = 1

Ez = np.zeros((Nx, Ny))
Hx = np.zeros((Nx, Ny))
Hy = np.zeros((Nx, Ny))


def update_fields():
  for i in range(Nx - 1):
    for j in range(Ny - 1):
      Hx[i, j] -= (dt / (mu0 * mu * dy)) * (Ez[i, j + 1] - Ez[i, j])
      Hy[i, j] += (dt  / (mu0 * mu * dx)) * (Ez[i + 1, j] - Ez[i, j])
      
      
  for i in range(1, Nx - 1):
    for j in range(1, Ny - 1):
      Ez[i, j] += (dt / (epsilon0 * epsilon * dx)) * (Hy[i, j] - Hy[i - 1, j]) - \
        (dt / (epsilon0 * epsilon * dy)) * (Hx[i, j] - Hx[i, j - 1])
        
        
        
        
def source(t: float):
  return np.exp(-((t - 30 * dt) / (5 * dt)) ** 2)


# for t in range(steps):
#     update_fields()
#     Ez[Nx // 2, Ny // 2] += source(t * dt) # Источник в центре

#     if t % 10 == 0:
#         plt.imshow(Ez, cmap='RdBu', extent=[0, Nx * dx, 0, Ny * dy])
#         plt.colorbar(label='Ez (V/m)')
#         plt.title(f'Время: {t * dt:.2e} с')
#         plt.xlabel('x (м)')
#         plt.ylabel('y (м)')
#         plt.pause(0.01)
#         plt.clf()

# plt.show()

# Магнитное поле
for t in range(steps):
    update_fields()
    Ez[Nx // 2, Ny // 2] += source(t * dt) # Источник в центре

    # Визуализация векторов магнитного поля
    if t % 10 == 0:
        plt.figure(figsize=(8, 8))
        print("Hx: ", Hx)
        print("Hy: ", Hy)
        skip = (slice(None, None, 5), slice(None, None, 5))  # Шаг для отображения векторов
        plt.quiver(np.arange(0, Nx * dx, dx)[skip[0]], np.arange(0, Ny * dy, dy)[skip[1]],
                   Hx[skip], Hy[skip], scale=10, scale_units='xy', angles='xy', color='k')
        plt.title(f'Magnetic Field (Hx, Hy) at t = {t * dt:.2e} s')
        plt.xlabel('x (м)')
        plt.ylabel('y (м)')
        plt.axis('equal')
        
        plt.pause(0.01)
        plt.clf()

plt.show()
