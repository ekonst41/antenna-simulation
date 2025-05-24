import numpy as np
import cmath
import math
from numpy import pi
import scipy.optimize as opt
import matplotlib.pyplot as plt
import numericalunits as nu
from typing import Callable, List, Tuple, Optional
from tmm.tmm import bc_matrix
from physics.modes import find_kzs
from system.optical_system import OpticalSystem

from tqdm import tqdm

INF = float('inf')

class ModeFinder:
    def __init__(self) -> None:
        """Инициализация класса для поиска мод в слоистых структурах"""
        self.nu = nu  # numerical units

    def inverse_fn(self, fn: Callable[[complex], complex], z: complex) -> float:
        """
        Вычисляет 1/fn(z) с обработкой деления на ноль
        
        Аргументы:
            fn: Функция комплексного переменного
            z: Точка в комплексной плоскости
            
        Возвращает:
            Обратное значение функции или INF при делении на ноль
        """
        f = fn(z)
        return INF if f == 0 else 1/f

    def create_search_grid(self, 
                         min_re: float, 
                         max_re: float, 
                         min_im: float, 
                         max_im: float, 
                         grid_points: int) -> Tuple[np.ndarray, np.ndarray, float, float]:
        """
        Создает 2D сетку для поиска в комплексной плоскости
        
        Аргументы:
            min_re: Минимальное значение действительной части
            max_re: Максимальное значение действительной части
            min_im: Минимальное значение мнимой части
            max_im: Максимальное значение мнимой части
            grid_points: Количество точек на сетке
            
        Возвращает:
            Кортеж (re_list, im_list, re_step, im_step) - списки значений и шаги
        """
        re_list, re_step = np.linspace(min_re, max_re, num=grid_points, retstep=True)
        im_list, im_step = np.linspace(min_im, max_im, num=grid_points, retstep=True)
        return re_list, im_list, re_step, im_step

    def calculate_contour_integral(self, 
                                 fn: Callable[[complex], complex], 
                                 z: complex, 
                                 d_re: float, 
                                 d_im: float) -> complex:
        """
        Аппроксимирует контурный интеграл 1/fn вокруг точки z
        
        Аргументы:
            fn: Функция комплексного переменного
            z: Центральная точка
            d_re: Полуширина прямоугольника по действительной оси
            d_im: Полувысота прямоугольника по мнимой оси
            
        Возвращает:
            Значение контурного интеграла
        """
        below = self.inverse_fn(fn, z - 1j * d_im)
        above = self.inverse_fn(fn, z + 1j * d_im)
        left = self.inverse_fn(fn, z - d_re)
        right = self.inverse_fn(fn, z + d_re)
        return (below * (2 * d_re) + right * (2j * d_im) + 
                above * (-2 * d_re) + left * (-2j * d_im))

    def plot_search_regions(self, 
                          fn: Callable[[complex], complex], 
                          min_re: float, 
                          max_re: float, 
                          min_im: float, 
                          max_im: float) -> Tuple[plt.Axes, plt.Axes]:
        """
        Генерирует диагностические графики для области поиска
        
        Аргументы:
            fn: Функция для анализа
            min_re, max_re: Границы по действительной оси
            min_im, max_im: Границы по мнимой оси
            
        Возвращает:
            Кортеж объектов осей matplotlib
        """
        res, ims, re_step, im_step = self.create_search_grid(
            min_re, max_re, min_im, max_im, 100)
        
        # График 1: log(|fn(z)|)
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        data = [[math.log10(abs(fn(re + 1j * im))) for re in res] for im in ims]
        ax1.imshow(data, extent=(min_re * self.nu.um, max_re * self.nu.um,
                                min_im * self.nu.um, max_im * self.nu.um),
                  origin='lower')
        ax1.set_xlabel('Re(kx) [rad/um]')
        ax1.set_ylabel('Im(kx) [rad/um]')
        ax1.set_title('log(|fn(z)|) -- Поиск минимумов (синие области)')

        # График 2: Контурные интегралы
        fig = plt.figure()
        ax2 = fig.add_subplot(111)
        data = [[-math.log10(abs(self.calculate_contour_integral(
            fn, re + 1j * im, re_step, im_step))) for re in res] for im in ims]
        ax2.imshow(data, extent=(min_re * self.nu.um, max_re * self.nu.um,
                                min_im * self.nu.um, max_im * self.nu.um),
                  origin='lower')
        ax2.set_xlabel('Re(kx) [rad/um]')
        ax2.set_ylabel('Im(kx) [rad/um]')
        ax2.set_title('Контурные интегралы выделяют нули функции fn(z)')

        return ax1, ax2

    def find_local_minima(self, 
                        fn: Callable[[complex], complex], 
                        region: List[float], 
                        grid_points: int, 
                        iteration_number: int, 
                        show_progress: bool) -> List[complex]:
        """
        Находит локальные минимумы в заданной области
        
        Аргументы:
            fn: Функция для анализа
            region: Границы области [min_re, max_re, min_im, max_im]
            grid_points: Количество точек на сетке
            iteration_number: Номер текущей итерации
            show_progress: Флаг вывода прогресса
            
        Возвращает:
            Список комплексных чисел - кандидатов в минимумы
        """
        min_re, max_re, min_im, max_im = region
        re_list, im_list, re_step, im_step = self.create_search_grid(
            min_re, max_re, min_im, max_im, grid_points)
        
        results_grid = np.array([[abs(fn(re + 1j * im)) for im in im_list] 
                               for re in re_list])
        
        local_mins = []
        for i in range(grid_points):
            for j in range(grid_points):
                if self.is_local_minimum(results_grid, i, j, grid_points):
                    local_mins.append((i, j))
        
        return self.filter_edge_points(local_mins, re_list, im_list, fn, 
                                     iteration_number, grid_points, show_progress)

    def is_local_minimum(self, 
                       grid: np.ndarray, 
                       i: int, 
                       j: int, 
                       size: int) -> bool:
        """
        Проверяет, является ли точка (i,j) локальным минимумом
        
        Аргументы:
            grid: Массив значений функции
            i, j: Индексы проверяемой точки
            size: Размер сетки
            
        Возвращает:
            True если точка является локальным минимумом
        """
        return all(grid[i2, j2] >= grid[i,j]
                  for i2 in [i-1, i, i+1]
                  for j2 in [j-1, j, j+1]
                  if (0 <= i2 < size and 0 <= j2 < size))

    def filter_edge_points(self, 
                         local_mins: List[Tuple[int, int]], 
                         re_list: np.ndarray, 
                         im_list: np.ndarray, 
                         fn: Callable[[complex], complex], 
                         iteration_number: int, 
                         grid_points: int, 
                         show_progress: bool) -> List[complex]:
        """
        Фильтрует точки на границах после начальных итераций
        
        Аргументы:
            local_mins: Список найденных минимумов (индексы)
            re_list, im_list: Списки значений координат
            fn: Анализируемая функция
            iteration_number: Номер текущей итерации
            grid_points: Размер сетки
            show_progress: Флаг вывода прогресса
            
        Возвращает:
            Отфильтрованный список минимумов (комплексные числа)
        """
        valid_mins = []
        for i, j in local_mins:
            z = re_list[i] + 1j * im_list[j]
            if iteration_number >= 2 and (i == 0 or j == 0 or 
                                         i == grid_points-1 or j == grid_points-1):
                if show_progress:
                    print('----')
                    print(f'Удаление граничной точки: (i,j)=({i},{j}), kx={z/self.nu.um**-1}, fn(z)={fn(z)}')
            else:
                valid_mins.append(z)
        return valid_mins

    def find_all_zeros(self, 
                      min_re: float, 
                      max_re: float, 
                      min_im: float, 
                      max_im: float, 
                      fn: Callable[[complex], complex],
                      grid_points: int = 20, 
                      iterations: int = 9, 
                      reduction_factor: int = 9,
                      plot_full_region: bool = True, 
                      show_progress: bool = False) -> List[complex]:
        """
        Основная функция для поиска всех нулей функции в комплексной плоскости
        
        Аргументы:
            min_re, max_re: Границы по действительной оси
            min_im, max_im: Границы по мнимой оси
            fn: Функция для анализа
            grid_points: Количество точек на сетке (по умолчанию 20)
            iterations: Количество итераций (по умолчанию 9)
            reduction_factor: Коэффициент уменьшения области (по умолчанию 9)
            plot_full_region: Флаг построения графиков (по умолчанию True)
            show_progress: Флаг вывода прогресса (по умолчанию False)
            
        Возвращает:
            Список найденных нулей (комплексные числа)
        """
        # Проверка входных данных
        assert reduction_factor > 1 and max_re > min_re and max_im > min_im
        assert grid_points > 2 * reduction_factor

        if plot_full_region:
            self.plot_search_regions(fn, min_re, max_re, min_im, max_im)

        regions = [[min_re, max_re, min_im, max_im]]
        region_width_re = max_re - min_re
        region_width_im = max_im - min_im

        all_zeros = []
        for iteration in tqdm(range(iterations)):
            current_zeros = []
            for region in regions:
                current_zeros.extend(self.find_local_minima(
                    fn, region, grid_points, iteration, show_progress))
            
            # Удаление дубликатов
            current_zeros = self.remove_duplicate_zeros(current_zeros, 
                                                      region_width_re/grid_points,
                                                      region_width_im/grid_points)
            
            if show_progress:
                print(f'Итерация {iteration}: Найдено кандидатов {len(current_zeros)}')

            # Подготовка к следующей итерации
            region_width_re /= reduction_factor
            region_width_im /= reduction_factor
            regions = [[z.real - region_width_re/2, z.real + region_width_re/2,
                       z.imag - region_width_im/2, z.imag + region_width_im/2]
                      for z in current_zeros]
            
            all_zeros = current_zeros

        return self.process_final_zeros(all_zeros)

    def remove_duplicate_zeros(self, 
                              zeros: List[complex], 
                              re_tol: float, 
                              im_tol: float) -> List[complex]:
        """
        Удаляет дубликаты нулей в пределах допуска
        
        Аргументы:
            zeros: Список нулей
            re_tol: Допуск по действительной оси
            im_tol: Допуск по мнимой оси
            
        Возвращает:
            Список уникальных нулей
        """
        unique_zeros = []
        for z in zeros:
            if not any(abs((z - uz).real) <= re_tol and 
                      abs((z - uz).imag) <= im_tol 
                      for uz in unique_zeros):
                unique_zeros.append(z)
        return unique_zeros

    def process_final_zeros(self, zeros: List[complex]) -> List[complex]:
        """
        Сортирует и фильтрует окончательные кандидаты в нули
        
        Аргументы:
            zeros: Список нулей
            
        Возвращает:
            Обработанный список нулей
        """
        zeros = sorted(zeros, key=lambda kx: abs(kx))
        
        # Удаление дубликатов с противоположными знаками
        i = 0
        while i < len(zeros) - 1:
            if abs(zeros[i] + zeros[i+1]) <= 1e-6 * (abs(zeros[i]) + abs(zeros[i+1])):
                zeros.pop(i)
            else:
                i += 1
        
        # Обеспечение правильной конвенции знаков
        return [(-z if (z.imag < 0 or (z.imag == 0 and z.real < 0)) else z)
                for z in zeros]

    def calculate_mode_error(self, 
                           input_params: OpticalSystem, 
                           kx: complex) -> float:
        """
        Вычисляет, насколько kx близко к валидной моде
        
        Аргументы:
            input_params: Параметры системы
            kx: Проверяемое волновое число
            
        Возвращает:
            Мера отклонения от условия моды
        """
        if kx == 0:
            return INF
        
        params = input_params.copy()
        params['kx'] = kx
        N = len(params['mu_list'])
        
        determinant = np.linalg.det(bc_matrix(find_kzs(params)))
        return determinant / kx**(N+1)  # Нормировка для удаления полюса в k=0

    def find_kx_modes(self, 
                     input_params: OpticalSystem, 
                     search_domain: Optional[List[float]] = None, 
                     show_progress: bool = False,
                     grid_points: int = 20, 
                     iterations: int = 9, 
                     reduction_factor: int = 9,
                     plot_full_region: bool = True) -> List[complex]:
        """
        Находит моды поверхностных плазмонов kx для заданных параметров
        
        Аргументы:
            input_params: Словарь параметров системы
            search_domain: Область поиска [min_re, max_re, min_im, max_im] (опционально)
            show_progress: Флаг вывода прогресса (по умолчанию False)
            grid_points: Количество точек на сетке (по умолчанию 20)
            iterations: Количество итераций (по умолчанию 9)
            reduction_factor: Коэффициент уменьшения области (по умолчанию 9)
            plot_full_region: Флаг построения графиков (по умолчанию True)
            
        Возвращает:
            Список найденных мод (комплексные числа)
        """
        w = input_params['w']
        d_list = input_params['d_list']
        ex_list = input_params['ex_list']
        ez_list = input_params['ez_list']
        mu_list = input_params['mu_list']
        
        N = len(mu_list)
        assert N == len(d_list) == len(ex_list) == len(ez_list)

        # Определение области поиска если не задана
        if search_domain is None:
            kx_re_max = max(
                max(abs((20 / (2 * pi * d_list[i])) * 
                    cmath.sqrt(ez_list[i] / ex_list[i])) 
                    for i in range(1, N)),
                3 * w / self.nu.c0)
            search_domain = [-kx_re_max, kx_re_max, 0, abs(kx_re_max)]

        # Поиск всех нулей функции ошибки
        error_fn = lambda kx: self.calculate_mode_error(input_params, kx)
        return self.find_all_zeros(
            *search_domain, error_fn,
            grid_points=grid_points,
            iterations=iterations,
            reduction_factor=reduction_factor,
            plot_full_region=plot_full_region,
            show_progress=show_progress)