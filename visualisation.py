# import pandas as pd
# import matplotlib.pyplot as plt

# import os

# # Загрузка данных из CSV файла
# file_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data/resonator_frequencies.txt')
# data = pd.read_csv(file_path)

# frequency = data['frequency1'].to_numpy()
# b_value = data['b_value'].to_numpy()
# freq_1 = frequency[3::16]
# freq_2 = frequency[6::16]
# freq_3 = frequency[9::16]
# freq_4 = frequency[12::16]
# freq_5 = frequency[15::16]
# b_value = b_value[::16]

# # Проверка загруженных данных
# print(data.head())

# # Построение графиков
# plt.figure(figsize=(10, 8))

# # График зависимости частоты от ширины
# plt.subplot(2, 1, 1)
# plt.plot(b_value, freq_1, label='m = 4, n = 1')
# plt.plot(b_value, freq_2, label='m = 3, n = 2')
# plt.plot(b_value, freq_3, label='m = 2, n = 3')
# plt.plot(b_value, freq_4, label='m = 1, n = 4')
# plt.plot(b_value, freq_5, label='m = 4, n = 4')

# frequency = data['frequency2'].to_numpy()
# b_value = data['b_value'].to_numpy()
# freq_1 = frequency[3::16]
# freq_2 = frequency[6::16]
# freq_3 = frequency[9::16]
# freq_4 = frequency[12::16]
# freq_5 = frequency[15::16]
# b_value = b_value[::16]
# plt.plot(b_value, freq_1, label='m = 4, n = 1')
# plt.plot(b_value, freq_2, label='m = 3, n = 2')
# plt.plot(b_value, freq_3, label='m = 2, n = 3')
# plt.plot(b_value, freq_4, label='m = 1, n = 4')
# plt.plot(b_value, freq_5, label='m = 4, n = 4')

# plt.title('Зависимость $\omega^2/(a^2c^2)$ от ширины (b)')
# plt.xlabel('Ширина (b), м')
# plt.ylabel('Частота (Гц)')
# plt.grid()
# plt.legend()

# # График зависимости мнимой части от ширины
# plt.subplot(2, 1, 2)
# plt.plot(data['b_value'], data['imaginary_part'], color='red', label='Мнимая часть')
# plt.title('Зависимость мнимой части от ширины (b)')
# plt.xlabel('Ширина (b), м')
# plt.ylabel('Мнимая часть')
# plt.grid()
# plt.legend()

# # Показать графики
# plt.tight_layout()
# output_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'avoided_crossing.png')
# plt.savefig(output_path)



import matplotlib.pyplot as plt
import pandas as pd
import os

# Чтение данных из файла
file_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data/resonator_frequencies.txt')
data = pd.read_csv(file_path)

# График
plt.figure(figsize=(10, 6))
plt.scatter(data["b_value"], data["frequency1"], label="Hybrid Mode 1", marker="o")
plt.scatter(data["b_value"], data["frequency2"], label="Hybrid Mode 2", marker="o")
plt.xlabel("Ширина резонатора (b)")
plt.ylabel("Частота (Гц)")
plt.title("Avoided Crossing in Resonator Frequencies")
plt.legend()
plt.grid()
plt.tight_layout()
output_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'avoided_crossing.png')
plt.savefig(output_path)

