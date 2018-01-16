import unittest, tempfile, os, shutil, glob, subprocess
from pyglib.basic.data import compare_array as compare_data


class KnowValues(unittest.TestCase):
    def test_r_matrix(self):
        self.assertTrue(compare_data('/IMPURITY_1/R'))

    def test_etot_matrix(self):
        self.assertTrue(compare_data('/etot_model'))

    def test_nphy_matrix(self):
        self.assertTrue(compare_data('/IMPURITY_1/NC_PHY'))



if __name__ == "__main__":
    print("test CyGutz results for cerium-soc only")
    cwd = os.getcwd()
    tempd = tempfile.mkdtemp(prefix='gtest_tmp_')
    for f in glob.glob('*.h5') + ['Ce.kgen']:
        shutil.copy(f, tempd)
    os.chdir(tempd)
    cmd = [os.environ['WIEN_GUTZ_ROOT2']+'/CyGutz']
    subprocess.call(cmd)
    unittest.main()
