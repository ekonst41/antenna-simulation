import numpy as np
import matplotlib.pyplot as plt
import imageio.v2 as imageio
import numpy as np
from tqdm import tqdm

frames = []

a = 1.0 # длина резонатора
h = 1.0 # высота резонатора
d = 0.2 # толщина волновода
l = 1.0 # длина волновода
E_0 = 0.01 # амплитуда поля
c = 3e8 # скорость света

omega = 10e10
k = omega / c # волновое число

steps_pl = 100

Nx = int(a * steps_pl + 2 * l * steps_pl)
Ny = int(h * steps_pl)

dx = 1.0 / steps_pl
dy = 1.0 / steps_pl
dt = 1e-11

E_x = np.zeros((Ny, Nx))
E_y = np.zeros((Ny, Nx))
E_x_prev = np.zeros((Ny, Nx))
E_y_prev = np.zeros((Ny, Nx))
E_x_next = np.zeros((Ny, Nx))
E_y_next = np.zeros((Nx, Nx))

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
    for i in range(Ny):
        for j in range(Nx):
          E_x[i, j] = starting_position(j * dx, i * dy)
          E_y[i, j] = starting_position(j * dx, i * dy)
          E[i, j] = np.sqrt(2) * starting_position(j * dx, i * dy)
    E_prev = E.copy()


def update_field():
    global E_x, E_y, E_x_prev, E_y_prev , E_x_next, E_y_next, E, E_prev, E_next

    for i in range(1, Ny - 1):
        for j in range(1, Nx - 1):
            if out(j * dx, i * dy):
                E_x[i, j] = 0
                E_y[i, j] = 0
                E[i, j] = 0
                continue
              
            d2E_x_dx2 = 1 / (dx ** 2) * (E_x[i, j + 1] + E_x[i, j - 1] - 2 * E_x[i, j])
            d2E_x_dy2 = 1 / (dy ** 2) * (E_x[i + 1, j] + E_x[i - 1, j] - 2 * E_x[i, j])
            E_x_next[i, j] = (d2E_x_dx2 + d2E_x_dy2 + k**2 * E_x[i, j]) / (k ** 2)
            
            d2E_y_dx2 = 1 / (dx ** 2) * (E_y[i, j + 1] + E_y[i, j - 1] - 2 * E_y[i, j])
            d2E_y_dy2 = 1 / (dy ** 2) * (E_y[i + 1, j] + E_y[i - 1, j] - 2 * E_y[i, j])
            E_y_next[i, j] = (d2E_y_dx2 + d2E_y_dy2 + k**2 * E_y[i, j]) / (k ** 2)
            
            E_next[i, j] = np.sqrt(E_x_next[i, j] ** 2 + E_y_next[i, j] ** 2)
            '''if flag == 1 and d2E_dx2 != 0 and d2E_dy2 != 0 and i % 10 == 0 and j % 10 == 0:
                print(i, j, d2E_dx2, d2E_dy2, E_next[i, j], sep='  ,  ')'''

    E_x_prev = E_x.copy()
    E_y_prev = E_y.copy()
    E_prev = E.copy()
    E_x = E_x_next.copy()
    E_y = E_y_next.copy()
    E = E_next.copy()



frames = []
steps = 10000
lines = np.array([(0, h/2 - d/2), (l, h/2 - d/2), (l, 0), (l + a, 0), (l + a, h/2 - d/2), (2*l + a, h/2 - d/2), (2*l + a, h/2 + d/2), (l + a, h/2 + d/2), (l + a, h), (l, h), (l, h/2 + d/2), (0, h/2 + d/2), (0, h/2 - d/2)])
lines = lines * steps_pl

initialize()
for t in tqdm(range(steps), desc="Итерации"):
    update_field()
    if t % 10 == 0:
        fig = plt.subplots(figsize=(9, 3))
        plt.imshow(E, cmap='hot', interpolation='bilinear', aspect='auto')  # Использование цветовой карты 'hot'
        x, y = zip(*lines)
        plt.plot(x, y, color='black')

        plt.colorbar(label='Значения')  # Добавление цветовой шкалы
        plt.title('Тепловая карта')  # Заголовок графика
        plt.xlabel('Ось X')  # Подпись оси X
        plt.gca().invert_yaxis()
        plt.ylabel('Ось Y')  # Подпись оси Y

        plt.savefig('temp.png')  # Сохраняем текущий график во временный файл
        plt.close()  # Закрываем текущее окно графика
        frames.append(imageio.imread('temp.png')) # Читаем изображение и добавляем его в список кадров

# Создание анимированного GIF из кадров
imageio.mimsave('heatmap3.gif', frames, duration=0.1)  # duration - время отображения каждого кадра в секундах


