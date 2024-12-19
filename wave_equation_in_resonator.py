import numpy as np

frames = []

a = 1.0 # длина резонатора
h = 1.0 # высота резонатора
d = 0.2 # толщина волновода
l = 1.0 # длина волновода
E_0 = 0.01 # амплитуда поля
c = 3e8 # скорость света

steps_pl = 100

Nx = int(a * steps_pl + 2 * l * steps_pl)
Ny = int(h * steps_pl)

dx = 1.0 / steps_pl
dy = 1.0 / steps_pl
dt = 1e-11

E = np.zeros((Ny, Nx))
E_prev = np.zeros((Ny, Nx))
E_next = np.zeros((Ny, Nx))


def out(x:float, y:float):
    if (x <= l or x >= l+a) and (y <= (h - d)/2 or y >= (h + d)/2):
        return True
    elif l < x < l+a and (y <= 0 or y >= h):
        return True
    return False


def starting_position(x:float, y:float):
    if out(x, y):
        return 0
    if x == 0 and y == h/2:
        return E_0
    return 0


def initialize():
    global E, E_prev
    for i in range(Ny):
        for j in range(Nx):
            E[i, j] = starting_position(j * dx, i * dy)
    E_prev = E.copy()


def update_field():
    global E, E_prev, E_next

    for i in range(1, Ny - 1):
        for j in range(1, Nx - 1):
            if out(j * dx, i * dy):
                E[i, j] = 0
                continue
            d2E_dx2 = 1 / (dx * dx) * (E[i, j + 1] + E[i, j - 1] - 2 * E[i, j])
            d2E_dy2 = 1 / (dy * dy) * (E[i + 1, j] + E[i - 1, j] - 2 * E[i, j])
            E_next[i, j] = (c * c * dt * dt) * (d2E_dx2 + d2E_dy2) + 2 * E[i, j] - E_prev[i, j]
            '''if flag == 1 and d2E_dx2 != 0 and d2E_dy2 != 0 and i % 10 == 0 and j % 10 == 0:
                print(i, j, d2E_dx2, d2E_dy2, E_next[i, j], sep='  ,  ')'''

    E_prev = E.copy()
    E = E_next.copy()



