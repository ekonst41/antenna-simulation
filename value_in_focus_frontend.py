import pandas as pd
import matplotlib.pyplot as plt

data = pd.read_csv('data/field_in_focus_005-040.csv')

focus_lengths = data['focus_lengths'].to_numpy()
field_values = data['field_values'].to_numpy()

plt.figure(figsize=(10, 6))
plt.plot(focus_lengths, field_values, marker='o', linestyle='-', color='b')
plt.title('Зависимость поля в фокусе от фокусного расстояния')
plt.xlabel(r'$ f $')
plt.ylabel(r'$ E $')
plt.grid(True)
plt.savefig('visuals/field_in_focus_005-040.png')
