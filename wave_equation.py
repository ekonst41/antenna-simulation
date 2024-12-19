import numpy as np
import matplotlib.pyplot as plt

Nx, Nz = 100, 100
a, c = 1, 1
dx, dz = 0.001, 0.001
dt = 1e-12

steps = 1000

epsilon0 = 8.854e-12
mu0 = 4 * np.pi * 1e-7
epsilon = 1
mu = 1
c_light = 1 / np.sqrt(epsilon0 * mu0)

n, p = 1, 1
omega = 2 * np.pi * 2.12e8
k = omega / c_light
E0 = 1

Ey = np.zeros((Nx, Nz))
Hx = np.zeros((Nx, Nz))
Hz = np.zeros((Nx, Nz))


def update_fields():
  for i in range(Nx - 1):
    for j in range(Nz - 1):
      Hx[i, j] -= (dt / (mu0 * mu * dz)) * (Ey[i, j + 1] - Ey[i, j])
      Hz[i, j] += (dt  / (mu0 * mu * dx)) * (Ey[i + 1, j] - Ey[i, j])
      
      
  for i in range(1, Nx - 1):
    for j in range(1, Nz - 1):
      Ey[i, j] += (dt / (epsilon0 * epsilon * dx)) * (Hz[i, j] - Hz[i - 1, j]) - \
        (dt / (epsilon0 * epsilon * dz)) * (Hx[i, j] - Hx[i, j - 1])
        
        
        
        
def source(x:float, z:float, t: float):
  frequency = 2.12e8 # (n, m, p) = (1, 1, 0) для размера 1м на 1м
  #return np.exp(-((t - 30 * dt) / (5 * dt)) ** 2)
  return E0 * np.sin(n * np.pi * x / a) * np.sin(p * np.pi * z / c) * np.cos(omega * t)


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
    Ey[Nx // 2, Nz // 2] += source(a / 2, c / 2, t * dt) # Источник в центре

    # Визуализация векторов магнитного поля
    if t % 10 == 0:
        plt.figure(figsize=(8, 8))
        print("Hx: ", Hx)
        print("Hy: ", Hz)
        skip = (slice(None, None, 5), slice(None, None, 5))  # Шаг для отображения векторов
        plt.quiver(np.arange(0, Nx * dx, dx)[skip[0]], np.arange(0, Nz * dz, dz)[skip[1]],
                   Hx[skip] * 100, Hz[skip] * 100, scale=10, scale_units='xy', angles='xy', color='k')
        plt.title(f'Magnetic Field (Hx, Hy) at t = {t * dt:.2e} s')
        plt.xlabel('x (м)')
        plt.ylabel('y (м)')
        plt.axis('equal')
        
        plt.pause(0.01)
        plt.clf()

plt.show()
