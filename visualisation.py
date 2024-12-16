import pandas as pd
import matplotlib.pyplot as plt

# Загрузка данных из CSV файла
file_path = 'resonator_frequencies.txt'  # Укажите путь к вашему файлу
data = pd.read_csv(file_path)

# Проверка загруженных данных
print(data.head())

# Построение графиков
plt.figure(figsize=(10, 8))

# График зависимости частоты от ширины
plt.subplot(2, 1, 1)
plt.plot(data['b_value'], data['frequency'], color='blue', label='Частота')
plt.title('Зависимость частоты от ширины (b)')
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
plt.show()
