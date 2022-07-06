import numpy as np


class ElementMatrix:
    def __init__(self, l, I, E, ro, A):
        self.l = l            # длина конечного элемента
        self.I = I
        self.E = E
        self.ro = ro          # плотность
        self.A = A            # площадь сечения
        self.L = self.l * 6   # длина балки
        self.EI = self.E * self.I

    # Составим матрицу жесткости балочного элемента
    def forming_stiffness_matrix_el(self):
        K_el = np.array([[12 * self.EI / self.l ** 3, 6 * self.EI / self.l ** 2, -12 * self.EI / self.l ** 3,
                          6 * self.EI / self.l ** 2],
                         [6 * self.EI / self.l ** 2, 4 * self.EI / self.l, -6 * self.EI / self.l ** 2,
                          2 * self.EI / self.l],
                         [-12 * self.EI / self.l ** 3, -6 * self.EI / self.l ** 2, 12 * self.EI / self.l ** 3,
                          -6 * self.EI / self.l ** 2],
                         [6 * self.EI / self.l ** 2, 2 * self.EI / self.l, -6 * self.EI / self.l ** 2,
                          4 * self.EI / self.l]])
        return K_el

    # Составим матрицу инерции балочного элемента
    def forming_inertia_matrix_el(self):
        M_el = np.array([[156, 22 * self.l, 54, -13 * self.l],
                         [22 * self.l, 4 * self.l ** 2, 13 * self.l, -3 * self.l ** 2],
                         [54, 13 * self.l, 156, -22 * self.l],
                         [-13 * self.l, -3 * self.l ** 2, -22 * self.l, 4 * self.l ** 2]])  # матрица инерции элемента
        M_el = self.ro * self.A * self.l / 420 * M_el
        return M_el


class GlobalMatrix:
    def __init__(self, matrix_el, n_sh):
        self.matrix_el = matrix_el
        self.matrix_glob = np.zeros([14, 14])  # сначала заполняем нулями глобальную матрицу
        self.n_sh = n_sh

    def forming_global_matrix(self, c=0, n_c=''):
        # Сформируем глобальную матрицу жесткости
        self.matrix_glob = np.zeros([14, 14])  # сначала заполняем нулями глобальную матрицу
        # Затем добавляем матрицы элементов (всего их 6)
        self.matrix_glob[0:4, 0:4] = self.matrix_el
        self.matrix_glob[2:6, 2:6] = self.matrix_glob[2:6, 2:6] + self.matrix_el
        self.matrix_glob[4:8, 4:8] = self.matrix_glob[4:8, 4:8] + self.matrix_el
        self.matrix_glob[6:10, 6:10] = self.matrix_glob[6:10, 6:10] + self.matrix_el
        self.matrix_glob[8:12, 8:12] = self.matrix_glob[8:12, 8:12] + self.matrix_el
        self.matrix_glob[10:14, 10:14] = self.matrix_glob[10:14, 10:14] + self.matrix_el

        for i in range(len(n_c)):
            self.matrix_glob[n_c[i], n_c[i]] = self.matrix_glob[n_c[i], n_c[i]] + c

        # Учтем граничные условия
        self.matrix_glob = np.delete(self.matrix_glob, self.n_sh, 0)
        self.matrix_glob = np.delete(self.matrix_glob, self.n_sh, 1)
        return self.matrix_glob



