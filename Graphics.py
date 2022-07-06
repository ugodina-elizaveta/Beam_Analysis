import matplotlib.pyplot as plt
import numpy as np


class FormsOscillation:
    def __init__(self, l, V1, V2, V3):
        self.l = l     # длина элемента балки
        self.V1 = V1   # вектор первой формы колебаний
        self.V2 = V2   # -\\- второй
        self.V3 = V3   # -\\- третьей
        self.x = 0
        self.N1 = 0
        self.N2 = 0
        self.N3 = 0
        self.N4 = 0

    # Участки для прогибов
    def finding_areas_of_oscillation(self, i):
        self.x = np.arange(i*self.l, (i+1)*self.l+0.1, 0.1)
        return self.x

    # Функции Эрмита
    def finding_the_hermite_function(self):
        self.finding_areas_of_oscillation(0)
        self.N1 = 1-3*self.x**2/self.l**2+2*self.x**3/self.l**3
        self.N2 = self.x-2*self.x**2/self.l+self.x**3/self.l**2
        self.N3 = 3*self.x**2/self.l**2-2*self.x**3/self.l**3
        self.N4 = -self.x**2/self.l+self.x**3/self.l**2
        return self.N1, self.N2, self.N3, self.N4

    # Функции прогиба для локальных координат
    def deflection_functions_for_local_coordinates(self, i):

        self.finding_the_hermite_function()
        # первая форма колебаний
        w1 = self.N1*self.V1[i]+self.N2*self.V1[i+1]+self.N3*self.V1[i+2]+self.N4*self.V1[i+3]
        w2 = self.N1*self.V2[i]+self.N2*self.V2[i+1]+self.N3*self.V2[i+2]+self.N4*self.V2[i+3]
        w3 = self.N1*self.V3[i]+self.N2*self.V3[i+1]+self.N3*self.V3[i+2]+self.N4*self.V3[i+3]
        return w1, w2, w3


class Graphics:
    def __init__(self, x, w):
        self.x = x
        self.w = w

    # Построение графиков форм колебаний для каждой частоты
    def plotting(self):
        plt.plot(self.x[: 31], self.w)
        plt.grid()
        plt.xlabel('Длина балки L, м')
        plt.ylabel('Форма колебаний, м')



