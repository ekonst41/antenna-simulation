# pip install Helmi-FEM
from helmi import Helmholtz
import skfem
import numpy as np

x_pts = np.linspace(0, 100, 101)
y_pts = np.linspace(-15, 15, 61)
mesh = skfem.MeshTri.init_tensor(x_pts, y_pts) # создание сетки
mesh = mesh.with_subdomains({'air': lambda x: x[0] < 50,
                             'plastic': lambda x: x[0] >= 50}) # границы волновода
mesh = mesh.with_boundaries({'bound_xmin': lambda x: np.isclose(x[0], x_pts[0]),
                             'bound_xmax': lambda x: np.isclose(x[0], x_pts[-1]),
                             'waveguide_min': lambda x: np.isclose(x[1], 5) & (x[0] <= 50),
                             'waveguide_max': lambda x: np.isclose(x[1], 10) & (x[0] <= 50),
                             'bound_ymin': lambda x: np.isclose(x[1], -15) & (x[0] >= 50),
                             'bound_ymax': lambda x: np.isclose(x[1], 15) & (x[0] >= 50),
                             'resonator_left_min': lambda x: np.isclose(x[0], 50) & (x[1] <= 5),
                             'resonator_left_max': lambda x: np.isclose(x[0], 50) & (x[1] >= 10)})

element = skfem.ElementTriP2()
fem = Helmholtz(mesh, element)

k0 = 0.1667
eps_air = 1
mu_air = 1
#eps_plastic = 2 - 0.1j
eps_plastic = 2
mu_plastic = 1
fem.assemble_subdomains(alpha={'air': 1 / mu_air,
                               'plastic': 1 / mu_plastic},
                        beta={'air': -1 * k0 ** 2 * eps_air,
                              'plastic': -1 * k0 ** 2 * eps_plastic},
                        f={'air': 0,
                           'plastic': 0}) # TM волна (см README.md)

fem.assemble_boundaries_dirichlet(value={'bound_ymin': 0,
                                         'bound_ymax': 0,
                                         'waveguide_min': 0,
                                         'waveguide_max': 0,
                                         'resonator_left_min': 0,
                                         'resonator_left_max': 0})

fem.assemble_boundaries_3rd(gamma={'bound_xmin': 1 / mu_air * 1j * k0,
                                   'bound_xmax': 1 / mu_plastic * 1j * k0},
                            q={'bound_xmin': 1 / mu_air * 2j * k0,
                               'bound_xmax': 0})

fem.solve()

from skfem.visuals.matplotlib import plot
import matplotlib.pyplot as plt

fig, ax = plt.subplots(2, 1,figsize=(10, 4))
plot(fem.basis, fem.phi_re, ax=ax[0])
plot(fem.basis, fem.phi_im, ax=ax[1])
ax[0].set_aspect(1)
ax[1].set_aspect(1)
ax[0].set_title('Real Part')
ax[1].set_title('Imaginary Part')
plt.tight_layout()
plt.savefig('./waveguide1.png')
plt.close()