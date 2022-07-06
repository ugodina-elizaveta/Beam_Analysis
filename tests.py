import unittest
from Matrix import *
from Аnalysis import *


class Test(unittest.TestCase):

    def test_frequencies(self):
        K_el_1 = ElementMatrix(3, 350/(100**4), 2*10**11, 7.870*10**3, 14.7/(100**2)).forming_stiffness_matrix_el()
        K_glob_1 = GlobalMatrix(K_el_1, [0, 1, 6]).forming_global_matrix(12003, [12])
        M_el_1 = ElementMatrix(3, 350/(100**4), 2*10**11, 7.870*10**3, 14.7/(100**2)).forming_inertia_matrix_el()
        M_glob_1 = GlobalMatrix(M_el_1, [0, 1, 6]).forming_global_matrix()
        self.assertEqual(['%.2f' % elem for elem in Analysis(K_glob_1, M_glob_1).finding_frequencies().tolist()],
                         ['%.2f' % elem for elem in [19.6409, 49.5472, 69.6142]])

        K_el_2 = ElementMatrix(3, 350*10**-8, 2*10**11, 7.8*10**3, 14.7*10**-4).forming_stiffness_matrix_el()
        K_glob_2 = GlobalMatrix(K_el_2, [0, 1]).forming_global_matrix(1920, [6, 12])
        M_el_2 = ElementMatrix(3, 350*10**-8, 2*10**11, 7.8*10**3, 14.7*10**-4).forming_inertia_matrix_el()
        M_glob_2 = GlobalMatrix(M_el_2, [0, 1]).forming_global_matrix()
        self.assertEqual(['%.2f' % elem for elem in Analysis(K_glob_2, M_glob_2).finding_frequencies().tolist()],
                         ['%.2f' % elem for elem in [6.7232,  18.4485,  47.5442]])

        K_el_3 = ElementMatrix(3.5, 572/(100**4), 2*10**11, 7.8*10**3, 17.4/(100**2)).forming_stiffness_matrix_el()
        K_glob_3 = GlobalMatrix(K_el_3, [0, 12]).forming_global_matrix((20*2*10**11*572/(100**4))/(3.5*6), [13])
        M_el_3 = ElementMatrix(3.5, 572/(100**4), 2*10**11, 7.8*10**3, 17.4/(100**2)).forming_inertia_matrix_el()
        M_glob_3 = GlobalMatrix(M_el_3, [0, 12]).forming_global_matrix()
        self.assertEqual(['%.2f' % elem for elem in Analysis(K_glob_3, M_glob_3).finding_frequencies().tolist()],
                         ['%.2f' % elem for elem in [9.3544, 30.6517, 64.7084]])

        K_el_4 = ElementMatrix(3, 350/(100**4), 2*10**11, 7.870*10**3, 14.7/(100**2)).forming_stiffness_matrix_el()
        K_glob_4 = GlobalMatrix(K_el_4, [0, 1, 12, 13]).forming_global_matrix(1928.439, [6])
        M_el_4 = ElementMatrix(3, 350/(100**4), 2*10**11, 7.870*10**3, 14.7/(100**2)).forming_inertia_matrix_el()
        M_glob_4 = GlobalMatrix(M_el_4, [0, 1, 12, 13]).forming_global_matrix()
        self.assertEqual(['%.2f' % elem for elem in Analysis(K_glob_4, M_glob_4).finding_frequencies().tolist()],
                         ['%.2f' % elem for elem in [17.664, 46.9142, 92.5519]])
