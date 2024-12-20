import numpy as np

frames = []

a = 1.0 # длина резонатора
h = 1.0 # высота резонатора
d = 0.2 # толщина волновода
l = 1.0 # длина волновода
E_0 = 0.1 # амплитуда поля
c = 3e8 # скорость света
omega = c * np.pi / a

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
    elif l < x < l + a and (x - l - a/2)**2 / (a/2)**2 + (y - h/2)**2 / ((h - d) / 2)**2 >= 1:
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


def update_field(t: float):
    global E, E_prev, E_next

    for i in range(1, Ny - 1):
        for j in range(1, Nx - 1):
            if out(j * dx, i * dy):
                E_next[i, j] = 0
                continue
            d2E_dx2 = 1 / (dx * dx) * (E[i, j + 1] + E[i, j - 1] - 2 * E[i, j])
            d2E_dy2 = 1 / (dy * dy) * (E[i + 1, j] + E[i - 1, j] - 2 * E[i, j])
            E_next[i, j] = (c * c * dt * dt) * (d2E_dx2 + d2E_dy2) + 2 * E[i, j] - E_prev[i, j]
    # E_next[Ny//2, 0] = E_0 * np.cos(omega * t)
    y = (h - d)/2
    while y < (h + d)/2:
        E_next[int(y / dy), 0] = E_next[int(y / dy), 1]
        E_next[int(y / dy), Nx - 1] = E_next[int(y / dy), Nx - 2]
        y += dy

    E_prev = E.copy()
    E = E_next.copy()

def calculate_intensity():
    x = 0
    y = (h - d)/2
    sum_intensity = 0.0
    while x < 2 * l + a:
        while y < (h + d)/2:
            sum_intensity += E[int(y / dy), int(x / dx)]**2
            y += dy
        y = (h - d) / 2
        x += dx
        if l + a > x >= l:
            x = l + a
    return sum_intensity

def calculate_intensity_in_resonator():
    x = l
    sum_intensity = 0
    while x < l + a:
        y = 0
        while y < h:
            sum_intensity += E[int(y / dy), int(x / dx)]**2
            y += dy
        x += dx
    return sum_intensity