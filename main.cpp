#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <iomanip>

const double c = 3e8; // Скорость света, м/с

// Функция для расчета резонансной частоты
double resonanceFrequency(double a, double b, int m, int n) {
    return (c / 2) * sqrt((m / a) * (m / a) + (n / b) * (n / b));
}

// Функция для расчета мнимой части резонанса
double imaginaryPart(double f) {
    return f * 0.1; // Примерная модель мнимой части
}

int main() {
    double a; // Высота резонатора
    std::cout << "Введите высоту резонатора (a): ";
    std::cin >> a;

    std::vector<double> b_values; // Ширина резонатора
    std::vector<double> frequencies; // Частоты резонатора
    std::vector<double> imaginary_parts; // Мнимые части резонансов

    for (double b = 0.01; b <= 1.0; b += 0.01) { // Изменение ширины от 0.01 до 1.0 м
        double f = resonanceFrequency(a, b, 1, 1); // Основной режим (m=1, n=1)
        b_values.push_back(b);
        frequencies.push_back(f);
        imaginary_parts.push_back(imaginaryPart(f));
    }

    // Запись результатов в файл
    std::ofstream out("resonator_frequencies.txt");
    out << std::fixed << std::setprecision(6);
    out << "b_value,frequency,imaginary_part\n";

    for (size_t i = 0; i < b_values.size(); ++i) {
        out << b_values[i] << "," << frequencies[i] << "," << imaginary_parts[i] << "\n";
    }

    out.close();

    std::cout << "Данные записаны в файл resonator_frequencies.txt." << std::endl;

    return 0;
}
