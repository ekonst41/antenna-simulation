import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

data = pd.read_csv('data/field_on_axe_different_focals_01-02.csv')
'''
focus = data['focus'].to_numpy()
real_coord = data['real_coord'].to_numpy()'''

focus = data['focus'].to_numpy()
real_coord = data['real_coord'].to_numpy()
#field = data['field'].to_numpy()
'''
time_max1 = time[np.argmax(field)]
time_max2 = time[np.argmax(field[field.size//2:]) + field.size//2]

fig, axs = plt.subplots(2, 1, figsize=(8, 6))
axs[0].plot(time, coord, linestyle='-', color='blue')
axs[0].axvline(x=time_max1, color='gray', linestyle='--')
axs[0].axvline(x=time_max2, color='gray', linestyle='--')
axs[0].axhline(y=0.4, color='red', linestyle='--')
axs[1].plot(time, field, linestyle='-', color='blue')
axs[1].axhline(y=np.max(field[field.size//2:]), color='red', linestyle='--')
axs[1].axhline(y=np.max(field), color='red', linestyle='--')
axs[1].axvline(x=time_max1, color='gray', linestyle='--')
axs[1].axvline(x=time_max2, color='gray', linestyle='--')
axs[0].set_title('Координата максимума от времени')
axs[0].set_xlabel('t')
axs[0].set_ylabel('coord')
axs[0].grid()

axs[1].set_title('Величина максимума от времени')
axs[1].set_xlabel('t')
axs[1].set_ylabel('field')
axs[1].grid()

plt.tight_layout()
plt.savefig('visuals/field_on_axe_040_focus_length.png')'''

plt.figure(figsize=(10, 6))
#plt.plot(time, field, linestyle='-', color='b')
plt.scatter(focus, focus, color='lightgray', label='геометрический фокус')
plt.errorbar(focus, focus, yerr=0.005, fmt='.', color='lightgray')
plt.scatter(focus, real_coord, label='реальный фокус')
#plt.plot(time, coord)
#plt.axhline(y=np.max(field), color='red', linestyle='--', label='y = 0.25')
plt.title('Координата реального фокуса в зависимости от фокального расстояния антенны')
plt.xlabel(r'координата фокуса')
plt.ylabel(r'реальная координата фокуса')
plt.grid(True)
plt.legend()
plt.savefig('visuals/field_on_axe_real_focus_01-02.png')