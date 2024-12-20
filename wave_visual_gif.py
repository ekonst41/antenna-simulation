import wave_equation_in_resonator as cmp
import matplotlib.pyplot as plt
import imageio.v2 as imageio
import numpy as np
from tqdm import tqdm


frames = []
steps = 3000
dt = cmp.dt
a = cmp.a
h = cmp.h
d = cmp.d
l = cmp.l
lines = np.array([(0, h/2 - d/2), (l, h/2 - d/2), (l, 0), (l + a, 0), (l + a, h/2 - d/2), (2*l + a, h/2 - d/2), (2*l + a, h/2 + d/2), (l + a, h/2 + d/2), (l + a, h), (l, h), (l, h/2 + d/2), (0, h/2 + d/2), (0, h/2 - d/2)])
lines = lines * cmp.steps_pl

# X, Y = np.meshgrid(np.arange(cmp.Nx),np.arange(cmp.Ny))
Z = np.array([[int(cmp.out(i * cmp.dx, j * cmp.dy)) for i in range(cmp.Nx)] for j in range(cmp.Ny)])

intensity_guide = []
intensity_res = []
time = []

cmp.initialize()
for t in tqdm(range(steps), desc="Итерации"):
    cmp.update_field(t * dt)

    intensity_guide.append(cmp.calculate_intensity())
    intensity_res.append(cmp.calculate_intensity_in_resonator())
    time.append(t * dt)

    if t % 10 == 0:
        fig = plt.subplots(figsize=(9, 3))
        plt.imshow(cmp.E, cmap='hot', interpolation='bilinear', aspect='auto')  # Использование цветовой карты 'hot'
        plt.contour(Z, colors='white')
        # x, y = zip(*lines)
        #plt.plot(x, y, color='black')


        # plt.colorbar(label='Значения')  # Добавление цветовой шкалы
        plt.title('Тепловая карта')  # Заголовок графика
        plt.xlabel('Ось X')  # Подпись оси X
        plt.gca().invert_yaxis()
        plt.ylabel('Ось Y')  # Подпись оси Y

        plt.savefig('temp.png')  # Сохраняем текущий график во временный файл
        plt.close()  # Закрываем текущее окно графика
        frames.append(imageio.imread('temp.png')) # Читаем изображение и добавляем его в список кадров

# Создание анимированного GIF из кадров
imageio.mimsave('heatmap10.gif', frames, duration=0.2)  # duration - время отображения каждого кадра в секундах

plt.plot(time, intensity_guide)
plt.plot(time, intensity_res)
plt.savefig('intensity1.png')


