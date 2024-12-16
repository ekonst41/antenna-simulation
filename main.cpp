#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <filesystem>

const double c = 3e8; // Скорость света, м/с
const double pi = 3.141592653589793;
double a = 4.0; // Длина резонатора
double kappa = 0.3; // Коэффициент взаимодействия между модами

// Функция для расчета резонансной частоты
double resonanceFrequency(double a, double b, int m, int n) {
    return c * pi * sqrt((m / a) * (m / a) + (n / b) * (n / b));
}

// Функция для расчета гибридных частот с учетом взаимодействия
std::vector<double> hybridFrequencies(double f1, double f2, double kappa) {
    double delta = f1 - f2;
    double sum = f1 + f2;
    double detuning = sqrt(delta * delta + 4 * kappa * kappa);

    double f_hybrid1 = (sum - detuning) / 2.0; // Первая гибридная частота
    double f_hybrid2 = (sum + detuning) / 2.0; // Вторая гибридная частота

    return {f_hybrid1, f_hybrid2};
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        std::cerr << "Использование: " << argv[0] << " <путь_к_корневой_директории>" << std::endl;
        return 1;
    }
    std::filesystem::path root_dir = argv[1];

    std::cout << "Calculating frequencies with avoided crossing..." << std::endl;

    std::vector<double> b_values; // Ширина резонатора
    std::vector<double> frequencies1; // Частоты первой моды
    std::vector<double> frequencies2; // Частоты второй моды
    std::vector<double> imaginary_parts; // Мнимые части резонансов

    for (double b = 3.0; b <= 5.0; b += 0.05) {
        for (int n = 4; n < 5; ++n) { // Пока только для двух частот, для большего числа нужно что-то с оператором Гамильтона делать (это квантфиз...)
            for (int m = 4; m < 5; ++m) {
                double f1 = resonanceFrequency(a, b, m + 2, n);
                double f2 = resonanceFrequency(a, b, m, n + 2);

                // Расчет гибридных частот
                std::vector<double> hybrid_freqs = hybridFrequencies(f1, f2, kappa);

                b_values.push_back(b);
                frequencies1.push_back(hybrid_freqs[0]);
                frequencies2.push_back(hybrid_freqs[1]);
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
    out << "b_value,frequency1,frequency2,imaginary_part\n";

    for (size_t i = 0; i < b_values.size(); ++i) {
        out << b_values[i] << "," << frequencies1[i] << "," << frequencies2[i] << "," << imaginary_parts[i] << "\n";
    }

    out.close();

    std::cout << "Данные записаны в файл resonator_frequencies.txt." << std::endl;

    return 0;
}