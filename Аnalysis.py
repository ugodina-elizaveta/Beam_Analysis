import scipy.sparse.linalg
import numpy as np


class Analysis:
    def __init__(self, K_glob, M_glob):
        self.K_glob = K_glob
        self.M_glob = M_glob
        # Решение частичной проблемы собственных значений
        self.X, self.D = scipy.sparse.linalg.eigsh(self.K_glob, 3, self.M_glob, sigma=0, which='LM')

    def finding_frequencies(self):
        for i in range(0, 3):
            W = np.sqrt(abs(self.X))   # i-ая собст. частота (рад/с)
            f = W / 2 / np.pi          # i-ая техн. частота  (Гц)
            print(f'{i+1}-ая собственная частота(рад/с):{W[i]}')
            print(f'{i+1}-ая техническая частота(Гц):{f[i]}')
            print()
        return W

    def formation_of_a_matrix_of_waveforms(self, n_sh):
        # Формирование матрицы форм колебаний с учетом граничных условий
        j = 1
        while j != len(n_sh)+1:
            A = np.insert(self.D, n_sh[j-1], [0, 0, 0], axis=0)
            self.D = A
            j += 1

        V1 = self.D[:, 0]      # вектор первой формы колебаний
        V2 = self.D[:, 1]      # -\\- второй
        V3 = self.D[:, 2]      # -\\- третьей
        return V1, V2, V3



