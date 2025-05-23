import matplotlib.pyplot as plt
from copy import deepcopy
import numericalunits as nu
import numpy as np

from physics.fields import calculate_Ex, calculate_Ez

def _get_z_values(params, z_range=None):
    """
    Возвращает массив координат z для построения графика.

    Args:
        params (dict): параметры моды
        z_range (list or None): пользовательский диапазон [z_min, z_max]

    Returns:
        np.ndarray: точки z для построения
    """
    kz_list = params['kz_list']
    layer_bottom_list = params['layer_bottom_list']
    N = len(kz_list)

    if z_range is not None:
        return np.linspace(z_range[0], z_range[1], num=200)

    if N == 2:
        # Для двух слоёв показываем затухание с обеих сторон от границы
        z_max = min(4 / abs(kz.imag) for kz in kz_list)
        z_min = -z_max
    else:
        # Для многослойных систем — центральная область
        z_max = 1.5 * layer_bottom_list[-1]
        z_min = -0.5 * layer_bottom_list[-1]

    return np.linspace(z_min, z_max, num=200)

def plot_mode(params, filename_x=None, filename_z=None, z_range=None):
    """
    Рисует профиль поля Ex и Ez вдоль оси z.

    Args:
        params (dict): параметры моды, содержащий kz_list, layer_bottom_list и т.п.
        filename_x (str or None): имя файла для сохранения Ex.
        filename_z (str or None): имя файла для сохранения Ez.
        z_range (list or None): диапазон построения [z_min, z_max].

    Returns:
        None
    """
    N = len(params['kz_list'])
    zs = _get_z_values(params, z_range)

    # Расчёт полей
    Exs = np.array([calculate_Ex(z, params) for z in zs])
    Ezs = np.array([calculate_Ez(z, params) for z in zs])

    max_E = max(np.abs(Exs).max(), np.abs(Ezs).max())
    Exs /= max_E
    Ezs /= max_E

    # === График Ex ===
    plt.figure(figsize=(10, 6))
    plt.plot(zs / nu.nm, Exs.real, label='Re(Ex)')
    plt.plot(zs / nu.nm, Exs.imag, '--', label='Im(Ex)')
    for i in range(1, N):
        plt.axvline(x=params['layer_bottom_list'][i] / nu.nm, color='gray', linestyle='--')
    plt.title('Профиль Ex в момент времени 0 и через четверть цикла')
    plt.xlabel('Позиция z (нм)')
    plt.ylabel('E(x) (произвольные единицы)')
    plt.grid(True)
    plt.legend()
    if filename_x is not None:
        plt.savefig(filename_x)
    plt.close()

    # === График Ez ===
    plt.figure(figsize=(10, 6))
    plt.plot(zs / nu.nm, Ezs.real, label='Re(Ez)')
    plt.plot(zs / nu.nm, Ezs.imag, '--', label='Im(Ez)')
    for i in range(1, N):
        plt.axvline(x=params['layer_bottom_list'][i] / nu.nm, color='gray', linestyle='--')
    plt.title('Профиль Ez в момент времени 0 и через четверть цикла')
    plt.xlabel('Позиция z (нм)')
    plt.ylabel('E(z) (произвольные единицы)')
    plt.grid(True)
    plt.legend()
    if filename_z is not None:
        plt.savefig(filename_z)
    plt.close()

def rescale_fields(factor, params):
    """
    Масштабирует амплитуду поля и энергии.

    Args:
        factor (complex): множитель, на который умножаются все амплитуды
        params (dict): параметры моды

    Returns:
        dict: обновлённые параметры с новыми масштабированными значениями
    """
    new_params = deepcopy(params)
    N = len(new_params['d_list'])

    for name in ['H_up_list', 'H_down_list', 'Ex_up_list', 'Ex_down_list', 'Ez_up_list', 'Ez_down_list']:
        for i in range(N):
            new_params[name][i] *= factor

    for i in range(N):
        new_params['Sx_list'][i] *= abs(factor)**2
    new_params['Sx_total'] *= abs(factor)**2

    return new_params