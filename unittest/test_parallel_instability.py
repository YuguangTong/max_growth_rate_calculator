import sys, os
module_dir = os.path.expanduser('~/Google Drive/research/projects/ion_anisotropy/anisotropy_instabilities/')
sys.path.append(module_dir)

import numpy.testing as npt
import numpy as np
import unittest
from find_start_freq import find_freq

class Test_paralleel(unittest.TestCase):

    def test_find_firehose_1(self):
        """
        Test find_freq() for parallel firehose instability.
        The test case data are taken from a figure in
        BA Maruca's PhD thesis (p.82). Maruca used SP Gary's code to generate
        the plot.

        """
        target_kz = 0.5 * np.sqrt(2)
        target_kp = 0
        target_beta = 3.16
        target_t = [1., 1., 4., 1.]
        target_a = [0.5, 1., 1., 1.]
        target_n = [1., 0, 0.05, 1.1]
        target_v = [0, 0, 0, 0]
        target_aol =1/5000
        param = [target_kz, target_kp, target_beta, target_t,
                 target_a, target_n, target_v, target_aol]

        freq = find_freq(param, wave_mode='pfh')
        npt.assert_allclose(np.real(freq),0.370972, rtol = 1e-3)
        npt.assert_allclose(np.imag(freq), 0.0241541, rtol = 1e-3)        

    def test_find_firehose_2(self):
        """
        Test find_freq() for parallel firehose instability.
        The test case data are taken from a figure in
        BA Maruca's PhD thesis (p.82). Maruca used SP Gary's code to generate
        the plot.

        """
        target_kz = 0.4 * np.sqrt(2)
        target_kp = 0
        target_beta = 3.16
        target_t = [1., 1., 4., 1.]
        target_a = [0.5, 1., 1., 1.]
        target_n = [1., 0, 0.05, 1.1]
        target_v = [0, 0, 0, 0]
        target_aol =1/5000
        param = [target_kz, target_kp, target_beta, target_t,
                 target_a, target_n, target_v, target_aol]

        freq = find_freq(param, wave_mode='pfh')
        npt.assert_allclose(np.real(freq), 0.248512, rtol = 1e-3)
        npt.assert_allclose(np.imag(freq), 0.0116537, rtol = 1e-3)            


if __name__ == '__main__':
    unittest.main()
