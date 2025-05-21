import yaml
import numericalunits as nu
from math import pi
import re

def _parse_complex(value):
    """Парсит комплексное число из строки в формате 'a + bj' или 'a + bj'"""
    if isinstance(value, (float, int)):
        return complex(value)
    if isinstance(value, str):
        # Удаляем все пробелы и приводим к нижнему регистру
        s = value.replace(" ", "").lower()
        # Проверяем, есть ли мнимая часть
        if 'j' in s:
            # Разделяем на реальную и мнимую части
            parts = re.split(r'[+-]', s)
            parts = [p for p in parts if p]  # Удаляем пустые строки
            real_part = 0.0
            imag_part = 0.0
            for part in parts:
                if 'j' in part:
                    imag_part = float(part.replace('j', ''))
                else:
                    real_part = float(part)
            # Учитываем знаки
            if '-' in s[s.find(real_part):]:
                real_part *= -1
            if '-' in s[s.find(imag_part):]:
                imag_part *= -1
            return complex(real_part, imag_part)
        else:
            return complex(float(s))
    return complex(value)

def load_config(config_path):
    with open(config_path, 'r') as f:
        data = yaml.full_load(f)

    config = data['system']
    N = len(config['layers'])

    w_parts = config['w'].split()
    wavelength = float(w_parts[0])
    if len(w_parts) == 1 or w_parts[1].lower() == 'nm':
        w = 2 * pi * nu.c0 / (wavelength * nu.nm)
    elif w_parts[1].lower() == 'um':
        w = 2 * pi * nu.c0 / (wavelength * nu.um)
    else:
        raise ValueError("Поддерживаемые единицы для длины волны: nm, um")

    d_list = []
    ex_list = []
    ez_list = []
    mu_list = []

    for layer in config['layers']:
        thickness = layer.get('thickness', 'inf')
        d = float('inf') if isinstance(thickness, str) and thickness.lower() == 'inf' else thickness * nu.nm
        
        ex = _parse_complex(layer['ex'])
        ez = _parse_complex(layer.get('ez', ex))
        mu = layer.get('mu', 1)

        d_list.append(d)
        ex_list.append(ex)
        ez_list.append(ez)
        mu_list.append(mu)

    return {
        'w': w,
        'd_list': d_list,
        'ex_list': ex_list,
        'ez_list': ez_list,
        'mu_list': mu_list
    }