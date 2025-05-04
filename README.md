Слои параллельны плоскости $x-y$. "Вверх" определяется как направление увеличения $z$. Я рассчитываю поверхностные плазмон-поляритоны (SPP), распространяющиеся вдоль оси $x$, и однородные в направлении $y$ ($k_y = 0$). Общее количество слоёв — $N$, они нумеруются от $0$ до $N-1$, где слои $0$ и $N-1$ имеют бесконечную толщину. Слой $0$ находится внизу ($z \ll 0$), а слой $N-1$ — сверху ($z \gg 0$).

Слой $m$ имеет (AC) диэлектрическую проницаемость $\varepsilon_{xm}$ в направлении $x$ и $\varepsilon_{zm}$ в направлении $z$, а также магнитную проницаемость $\mu_{ym}$. Хотя мы не предполагаем изотропность диэлектрической или магнитной проницаемости, мы предполагаем, что недиагональные элементы, такие как $\varepsilon_{xz}$, равны нулю.

Весь текст представлен в системе СИ. Диэлектрическая и магнитная проницаемости безразмерны (по сравнению с $\varepsilon_0$ или $\mu_0$), за исключением самих $\varepsilon_0$ и $\mu_0$. Все волновые числа и волновые векторы — угловые.

## Формулы

$k_x$ — комплексное волновое число в плоскости. Мы заранее не знаем его значение; его нужно определить.

Волновые числа $k_{zm}$ выражаются через $k_x$ следующим образом:

$$
k_{zm} = \pm \sqrt{\omega^2 \mu_{ym} \varepsilon_{xm} / c^2 - (\varepsilon_{xm} / \varepsilon_{zm}) k_x^2}
$$
(выберите корень с неотрицательной мнимой частью)

> Если $k_{zm}$ вещественное, то можно выбрать любой знак, это не имеет значения. Единственное место, где это может быть важно — полубесконечные слои, но в этом случае $k_{zm}$ никогда не будет вещественным, иначе волна не будет локализована.

Когда я записываю формулу для $\vec{E}(z)$ или $\vec{H}(z)$, подразумевается, что её нужно умножить на $e^{i k_x x - i \omega t}$ и взять действительную часть.

Обсуждение основано главным образом на $H$-поле, поскольку оно скалярное (направлено по оси $y$), в отличие от электрического поля, которое имеет две компоненты. Иногда я могу использовать обозначение $H(z)$ вместо $H_y(z)$. Для слоя $m$:

$$
H_y(z) = H_{m\uparrow} e^{i k_{zm} (z-z_{\text{низ слоя } m})} + H_{m\downarrow} e^{i k_{zm} (z_{\text{верх слоя } m} - z)}
$$

$$
E_x(z) = E_{xm\uparrow} e^{i k_{zm} (z-z_{\text{низ слоя } m})} + E_{xm\downarrow} e^{i k_{zm} (z_{\text{верх слоя } m} - z)}
$$

$$
E_z(z) = E_{zm\uparrow} e^{i k_{zm} (z-z_{\text{низ слоя } m})} + E_{zm\downarrow} e^{i k_{zm} (z_{\text{верх слоя } m} - z)}
$$

где

$$
E_{xm\uparrow} = \frac{H_{m\uparrow} k_{zm}}{\omega \varepsilon_{xm} \varepsilon_0}, \quad
E_{xm\downarrow} = -\frac{H_{m\downarrow} k_{zm}}{\omega \varepsilon_{xm} \varepsilon_0}
$$

$$
E_{zm\uparrow} = -\frac{H_{m\uparrow} k_x}{\omega \varepsilon_{zm} \varepsilon_0}, \quad
E_{zm\downarrow} = -\frac{H_{m\downarrow} k_x}{\omega \varepsilon_{zm} \varepsilon_0}
$$

(при этом $X_{m\uparrow}=0$ в слое $0$ (нет нижней границы) и $X_{m\downarrow}=0$ в слое $N-1$ (нет верхней границы)). $X_{m\uparrow}$ описывает компоненту, затухающую при увеличении $z$, а $X_{m\downarrow}$ — при уменьшении $z$.

### Проверка закона Фарадея:

$$
-i\omega (\mu_y \mu_0 H) = \partial_t B \stackrel{?}{=} -\nabla \times E \;\; \rightarrow \;\; H \stackrel{?}{=} (-i/\mu_0\mu_y\omega)\nabla \times E
$$

$$
H_y \stackrel{?}{=} \frac{-i}{\mu_0 \mu_{ym} \omega} \left(\partial_z E_x - \partial_x E_z\right)
= \frac{H_{m\uparrow} e^{i k_{zm} (z-z_{\text{низ слоя } m})}}{\mu_0 \mu_{ym} \omega^2 \varepsilon_0} \left(\frac{k_{zm}^2}{\varepsilon_{xm}} + \frac{k_x^2}{\varepsilon_{zm}}\right)
+ \frac{H_{m\downarrow} e^{i k_{zm} (z_{\text{верх слоя } m} - z)}}{\mu_0 \mu_{ym} \omega^2 \varepsilon_0} \left(\frac{k_{zm}^2}{\varepsilon_{xm}} + \frac{k_x^2}{\varepsilon_{zm}}\right)
$$

✅ Работает!

### Проверка закона Ампера:

$$
\nabla \times H \stackrel{?}{=} -i \omega (\varepsilon \varepsilon_0) E \;\; \rightarrow \;\; E \stackrel{?}{=} \frac{i}{\omega \varepsilon \varepsilon_0} \nabla \times H
$$

$$
E_z \stackrel{?}{=} \frac{i}{\omega \varepsilon_{zm} \varepsilon_0} \partial_x H_y
= \frac{i}{\omega \varepsilon_{zm} \varepsilon_0} \left( i k_x H_{m\uparrow} e^{i k_{zm} (z-z_{\text{низ слоя } m})} + i k_x H_{m\downarrow} e^{i k_{zm} (z_{\text{верх слоя } m} - z)}\right)
$$

$$
E_x \stackrel{?}{=} \frac{-i}{\omega \varepsilon_{xm} \varepsilon_0} \partial_z H_y
= \frac{-i}{\omega \varepsilon_{xm} \varepsilon_0} \left( i k_{zm} H_{m\uparrow} e^{i k_{zm} (z-z_{\text{низ слоя } m})} - i k_{zm} H_{m\downarrow} e^{i k_{zm} (z_{\text{верх слоя } m} - z)}\right)
$$

✅ Работает!

### Проверка закона Гаусса:

$$
\nabla \cdot \vec{D} \stackrel{?}{=} 0 \;\; \rightarrow \;\; \varepsilon_x \partial_x E_x + \varepsilon_z \partial_z E_z \stackrel{?}{=} 0
$$

$$
0 \stackrel{?}{=} \left( i \varepsilon_{xm} k_x E_{xm\uparrow} e^{i k_{zm} (z-z_{\text{низ слоя } m})} + i \varepsilon_{xm} k_x E_{xm\downarrow} e^{i k_{zm} (z_{\text{верх слоя } m} - z)} \right)
+ \left( i \varepsilon_{zm} k_{zm} E_{zm\uparrow} e^{i k_{zm} (z-z_{\text{низ слоя } m})} - i \varepsilon_{zm} k_{zm} E_{zm\downarrow} e^{i k_{zm} (z_{\text{верх слоя } m} - z)} \right)
$$

✅ Работает!

## Стратегия решения

Я *угадываю* $k_x$. Затем вычисляю все $k_{zm}$. У меня есть $2N-2$ неизвестных (все $H_{m\uparrow}, H_{m\downarrow}$, кроме $H_{0\downarrow}$ и $H_{N-1,\uparrow}$) и $(N-1)$ границ, каждая из которых даёт два уравнения непрерывности ($E_x$ непрерывно и $\varepsilon_z E_z$ непрерывно). Я думаю, что $H_y$ тоже непрерывно, но это избыточно по сравнению с другими двумя. Таким образом, это система линейных уравнений с нетривиальным решением. Существует связанная матрица, определитель которой должен быть равен нулю. Я могу вычислить этот определитель для каждого возможного $k_x$ и использовать его как меру качества для поиска реального решения.

### Непрерывность $E_x$:

$$
E_{x0\downarrow} = E_{x1\uparrow} + E_{x1\downarrow} e^{i k_{z1} d_1}
$$

$$
E_{x1\uparrow} e^{i k_{z1} d_1} + E_{x1\downarrow} = E_{x2\uparrow} + E_{x2\downarrow} e^{i k_{z2} d_2}
$$

$$
\cdots
$$

$$
E_{x(N-2)\uparrow} e^{i k_{z(N-2)} d_{(N-2)}} + E_{x(N-2)\downarrow} = E_{x(N-1)\uparrow}
$$

### Непрерывность $\varepsilon_z E_z$:

$$
\varepsilon_{z0} E_{z0\downarrow} = \varepsilon_{z1} E_{z1\uparrow} + \varepsilon_{z1} E_{z1\downarrow} e^{i k_{z1} d_1}
$$

$$
\varepsilon_{z1} E_{z1\uparrow} e^{i k_{z1} d_1} + \varepsilon_{z1} E_{z1\downarrow} = \varepsilon_{z2} E_{z2\uparrow} + \varepsilon_{z2} E_{z2\downarrow} e^{i k_{z2} d_2}
$$

$$
\cdots
$$

$$
\varepsilon_{z(N-1)} E_{z(N-2)\uparrow} e^{i k_{z(N-2)} d_{(N-2)}} + \varepsilon_{z(N-2)} E_{z(N-2)\downarrow} = \varepsilon_{z(N-1)} E_{z(N-1)\uparrow}
$$

Для краткости, пусть $\delta_m = e^{i k_{zm} d_m}$. Тогда матрица системы:

$$
\begin{pmatrix}
\frac{E_{x0\downarrow}}{H_{0\downarrow}} & -\frac{E_{x1\uparrow}}{H_{1\uparrow}} & -\frac{E_{x1\downarrow}}{H_{1\downarrow}} \delta_1 & 0 & 0 & 0 \\
0 & \frac{E_{x1\uparrow}}{H_{1\uparrow}} \delta_1 & \frac{E_{x1\downarrow}}{H_{1\downarrow}} & -\frac{E_{x2\uparrow}}{H_{2\uparrow}} & -\frac{E_{x2\downarrow}}{H_{2\downarrow}} \delta_2 & 0 \\
0 & 0 & 0 & \frac{E_{x2\uparrow}}{H_{2\uparrow}} \delta_2 & \frac{E_{x2\downarrow}}{H_{2\downarrow}} & -\frac{E_{x3\uparrow}}{H_{3\uparrow}} \\ \hline
\varepsilon_{z0}\frac{E_{z0\downarrow}}{H_{0\downarrow}} & -\varepsilon_{z1}\frac{E_{z1\uparrow}}{H_{1\uparrow}} & -\varepsilon_{z1}\frac{E_{z1\downarrow}}{H_{1\downarrow}} \delta_1 & 0 & 0 & 0 \\
0 & \varepsilon_{z1} \frac{E_{z1\uparrow}}{H_{1\uparrow}} \delta_1 & \varepsilon_{z1} \frac{E_{z1\downarrow}}{H_{1\downarrow}} & -\varepsilon_{z2} \frac{E_{z2\uparrow}}{H_{2\uparrow}} & -\varepsilon_{z2} \frac{E_{z2\downarrow}}{H_{2\downarrow}} \delta_2 & 0 \\
0 & 0 & 0 & \varepsilon_{z2} \frac{E_{z2\uparrow}}{H_{2\uparrow}} \delta_2 & \varepsilon_{z2} \frac{E_{z2\downarrow}}{H_{2\downarrow}} & -\varepsilon_{z3} \frac{E_{z3\uparrow}}{H_{3\uparrow}}
\end{pmatrix}
\begin{pmatrix} H_{0\downarrow} \\ H_{1\uparrow} \\ H_{1\downarrow} \\ H_{2\uparrow} \\ H_{2\downarrow} \\ H_{3\uparrow} \end{pmatrix} = \begin{pmatrix} 0 \\ 0 \\ 0 \\ 0 \\ 0 \\ 0 \end{pmatrix}
$$

(Горизонтальная линия разделяет уравнения для $E_x$ и $E_z$.)

## Вектор Пойнтинга

Усреднённый по времени вектор Пойнтинга: $S = \frac{1}{2} E \times H^*$ (см. Jackson (6.132)). Действительная часть $S$ указывает средний поток мощности. Меня интересует только $x$-компонента $S$, $S_x = -(1/2)E_z H_y^*$.

$$
S_x = -(1/2)E_z H_y^*
= -(1/2) (\stack{E_{zm\uparrow} e^{i k_{zm} (z-z_{\text{низ слоя } m})} \; + \; }{\; + \; E_{zm\downarrow} e^{i k_{zm} (z_{\text{верх слоя } m} - z)}})(\stack{H_{m\uparrow}^* e^{-i k_{zm}^* (z-z_{\text{низ слоя } m})} \; + \; }{\; + \; H_{m\downarrow}^* e^{-i k_{zm}^* (z_{\text{верх слоя } m} - z)}})
$$

$$
= \frac{-E_{zm\uparrow}H_{m\uparrow}^*}{2}e^{-2 \Im(k_{zm})(z-z_{\text{низ слоя } m})} + \frac{-E_{zm\downarrow}H_{m\downarrow}^*}{2}e^{-2 \Im(k_{zm})(z_{\text{верх слоя } m}-z)}
$$

$$
+ \frac{-E_{zm\downarrow}H_{m\uparrow}^*}{2} e^{ik_{zm} d_m}e^{-2i\Re(k_{zm})(z-z_{\text{низ слоя } m})} + \frac{-E_{zm\uparrow}H_{m\downarrow}^*}{2} e^{ik_{zm} d_m}e^{-2i\Re(k_{zm})(z_{\text{верх слоя } m}-z)}
$$

Проинтегрируем:

$$
\int_{z_{\text{низ слоя } m}}^{z_{\text{верх слоя } m}} S_x = \frac{-E_{zm\uparrow}H_{m\uparrow}^*}{4 \Im(k_{zm})}(1-e^{-2\Im(k_{zm})d_m}) + \frac{-E_{zm\downarrow}H_{m\downarrow}^*}{4 \Im(k_{zm})}(1-e^{-2\Im(k_{zm})d_m})
$$

$$
+ \frac{-E_{zm\downarrow}H_{m\uparrow}^*}{4i\Re(k_{zm})} e^{i k_{zm} d_m} (1 - e^{-2i\Re(k_{zm})d_m}) + \frac{-E_{zm\uparrow}H_{m\downarrow}^*}{4i \Re(k_{zm})} e^{ik_{zm} d_m}(1 - e^{-2i\Re(k_{zm})d_m})
$$
