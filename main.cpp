#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <filesystem>

const double c = 3e8; // Скорость света, м/с
const double pi = 3.141592653589793;
double a = 4.0;

// Функция для расчета резонансной частоты
double resonanceFrequency(double a, double b, int m, int n) {
    return c * pi * sqrt((m / a) * (m / a) + (n / b) * (n / b));
}


int main(int argc, char *argv[]) {
    if (argc < 2) {
        std::cerr << "Использование: " << argv[0] << " <путь_к_корневой_директории>" << std::endl;
        return 1;
    }
     std::filesystem::path root_dir = argv[1];

    std::cout << "Calculating frequencies..." << std::endl;

    std::vector<double> b_values; // Ширина резонатора
    std::vector<double> frequencies; // Частоты резонатора
    std::vector<double> imaginary_parts; // Мнимые части резонансов

    for (double b = 3.0; b <= 5.0; b += 0.05) { // Изменение ширины от 0.01 до 1.0 м
        for (int n = 1; n < 5; n++) {
            for (int m = 1; m < 5; m++) {
                double f = resonanceFrequency(a, b, m, n);
                b_values.push_back(b);
                frequencies.push_back(f * f / (a * a * c * c));
                imaginary_parts.push_back(0);
            }
        }
    }

    std::filesystem::path data_dir = root_dir / "data";

    if (!std::filesystem::exists(data_dir)) {
        std::filesystem::create_directory(data_dir);
        std::cout << "Data directory has been just created!" << std::endl;
    }

    // Запись результатов в файл
    std::ofstream out(data_dir / "resonator_frequencies.txt");
    out << std::fixed << std::setprecision(6);
    out << "b_value,frequency,imaginary_part\n";

    for (size_t i = 0; i < b_values.size(); ++i) {
        out << b_values[i] << "," << frequencies[i] << "," << imaginary_parts[i] << "\n";
    }

    out.close();

    std::cout << "Данные записаны в файл resonator_frequencies.txt." << std::endl;

    return 0;
}
