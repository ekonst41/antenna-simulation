import pandas as pd
import matplotlib.pyplot as plt

# Загрузка данных из CSV файла
file_path = 'cmake-build-debug/resonator_frequencies.txt'  # Укажите путь к вашему файлу
data = pd.read_csv(file_path)

frequency = data['frequency'].to_numpy()
b_value = data['b_value'].to_numpy()
freq_1 = frequency[3::16]
freq_2 = frequency[6::16]
freq_3 = frequency[9::16]
freq_4 = frequency[12::16]
freq_5 = frequency[15::16]
b_value = b_value[::16]

# Проверка загруженных данных
print(data.head())

# Построение графиков
plt.figure(figsize=(10, 8))

# График зависимости частоты от ширины
plt.subplot(2, 1, 1)
plt.plot(b_value, freq_1, color='blue', label='m = 4, n = 1')
plt.plot(b_value, freq_2, color='blue', label='m = 3, n = 2')
plt.plot(b_value, freq_3, color='blue', label='m = 2, n = 3')
plt.plot(b_value, freq_4, color='blue', label='m = 1, n = 4')
plt.plot(b_value, freq_5, color='blue', label='m = 4, n = 4')
plt.title('Зависимость $\omega^2/(a^2c^2)$ от ширины (b)')
plt.xlabel('Ширина (b), м')
plt.ylabel('Частота (Гц)')
plt.grid()
plt.legend()

# График зависимости мнимой части от ширины
plt.subplot(2, 1, 2)
plt.plot(data['b_value'], data['imaginary_part'], color='red', label='Мнимая часть')
plt.title('Зависимость мнимой части от ширины (b)')
plt.xlabel('Ширина (b), м')
plt.ylabel('Мнимая часть')
plt.grid()
plt.legend()

# Показать графики
plt.tight_layout()
plt.savefig('output.png')
# plt.show()
