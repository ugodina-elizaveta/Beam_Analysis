import copy
import logging
from Matrix import *
from Аnalysis import *
from Graphics import *
import scipy.sparse.linalg


logger = logging.getLogger('user')
logger.setLevel(logging.INFO)
handler = logging.FileHandler('user_log.txt')
logger_math = logging.getLogger('math')
logger_math.setLevel(logging.INFO)
handler_math = logging.FileHandler('math_log.txt')
handler.setFormatter(logging.Formatter(fmt='[%(asctime)s: %(levelname)s] %(message)s'))
handler_math.setFormatter(logging.Formatter(fmt='[%(asctime)s: %(levelname)s] %(message)s'))
logger.addHandler(handler)
logger_math.addHandler(handler_math)

logger.info('the program is running')

# Входные данные
l = 3               # длина конечного элемента
logger.info('the user entered the FE length: ' + str(l))
L = l*6             # длина балки
I = 572/(100**4)
logger.info('the user entered the moment of inertia: ' + str(I))
E = 2*10**11
logger.info('the user has entered the modulus of elasticity: ' + str(E))
EI = E*I
ro = 7.8*10**3      # плотность
logger.info('the user entered the density: ' + str(ro))
A = 17.4/(100**2)   # площадь сечения
logger.info('the user entered the cross-sectional area: ' + str(A))
c = 2353            # жесткость пружины
logger.info('the user entered the stiffness: ' + str(c))
n_c = (0,)           # номер элемента с пружиной
logger.info('the user entered the number of the spring element: ' + str(n_c))
n_sh = (6, 12)      # номера элементов с шарниром
logger.info('the user entered the number of the element with the hinge: ' + str(n_sh))
n_sh_ = copy.copy(n_sh)

n_sh = list(n_sh)
n_c = list(n_c)

if (l and I and E and ro and A and c) >= 0:
    for i in range(len(n_c)):
        if 0 < n_c[i] < 14 and type(n_c[i]) in [int, float]:
            logger_math.info('the user is input data is correct')
    for i in range(len(n_sh)):
        if 0 < n_sh[i] < 14 and type(n_sh[i]) in [int, float]:
            logger_math.info('the user is input data is correct')
else:
    logger_math.warning('the user is input data is incorrect')
    raise ValueError("input data is incorrect")


def graphics(l, V1, V2, V3, c):
    i = 0
    j = 0
    while i != 6:
        if i == 5:
            x = np.arange(i*l, (i+1)*l+0.01, 0.1)
        else:
            x = FormsOscillation(l, V1, V2, V3).finding_areas_of_oscillation(i)

        w1, w2, w3 = FormsOscillation(l, V1, V2, V3).deflection_functions_for_local_coordinates(j)
        plt.figure(1)
        Graphics(x, w1).plotting()
        plt.title(f'Первая форма колебаний при с={c}')
        plt.figure(2)
        Graphics(x, w2).plotting()
        plt.title(f'Вторая форма колебаний при с={c}')
        plt.figure(3)
        Graphics(x, w3).plotting()
        plt.title(f'Третья форма колебаний при с={c}')
        i += 1
        j += 2
    plt.show()


def checking_frequencies(W, c):
    flag = False
    for i in range(0, len(W)):
        if W[i] > 0:
            flag = True
    if flag:
        logger_math.info(f'frequencies at c = {c} found correctly')
    else:
        logger_math.info(f'frequencies at c = {c} found incorrectly')


K_el = ElementMatrix(l, I, E, ro, A).forming_stiffness_matrix_el()
K_glob = GlobalMatrix(K_el, n_sh).forming_global_matrix(c, n_c)

M_el = ElementMatrix(l, I, E, ro, A).forming_inertia_matrix_el()
M_glob = GlobalMatrix(M_el, n_sh).forming_global_matrix()

W_c_const = Analysis(K_glob, M_glob).finding_frequencies()
V1, V2, V3 = Analysis(K_glob, M_glob).formation_of_a_matrix_of_waveforms(n_sh)
checking_frequencies(W_c_const, 'const')
graphics(l, V1, V2, V3, 'const')


print('ПРЕДЕЛЬНЫЕ СЛУЧАИ')
print('Пусть с => 0')
M_c0 = M_glob
K_c0 = GlobalMatrix(K_el, n_sh).forming_global_matrix(0, n_c)

W_c_0 = Analysis(K_c0, M_c0).finding_frequencies()
V1_c0, V2_c0, V3_c0 = Analysis(K_c0, M_c0).formation_of_a_matrix_of_waveforms(n_sh)
checking_frequencies(W_c_0, '0')
graphics(l, V1_c0, V2_c0, V3_c0, '0')


print('Пусть с => inf')
for i in range(len(n_c)):
    n_sh.append(n_c[i])

n_sh_inf = sorted(n_sh)
M_c_inf = GlobalMatrix(M_el, n_sh_inf).forming_global_matrix()
K_c_inf = GlobalMatrix(K_el, n_sh_inf).forming_global_matrix(0, n_c)

W_c_inf = Analysis(K_c_inf, M_c_inf).finding_frequencies()
V1_c_inf, V2_c_inf, V3_c_inf = Analysis(K_c_inf, M_c_inf).formation_of_a_matrix_of_waveforms(n_sh_inf)
checking_frequencies(W_c_inf, 'inf')
graphics(l, V1_c_inf, V2_c_inf, V3_c_inf, 'inf')

# Исследуем влияние жёсткости опорной связи на первую собственную частоту изгибных колебаний.
# Построим график зависимости частоты от жёсткости.
WW = []
cc = 1000*np.array([0, 2.353, 4, 10, 20, 40, 60, 100, 150, 200, 300, 400, 600, 800, 1000])  # массив жесткостей пружины
for i in range(0, 15):
    K = GlobalMatrix(K_el, n_sh_).forming_global_matrix(cc[i], n_c)
    X, D = scipy.sparse.linalg.eigsh(K, 1, M_glob, sigma=0, which='LM')
    W = np.sqrt(abs(X))
    WW.append(W)
plt.figure(9)
asimptota_x = [0, 10**6]
asimptota_y = [W_c_inf[0], W_c_inf[0]]
plt.plot(cc, WW)
plt.plot(asimptota_x, asimptota_y)
plt.grid()
plt.xlabel('Жесткость, Н/м')
plt.ylabel('Собственная частота, рад/с')
plt.show()

logger.info('the program is completed')
