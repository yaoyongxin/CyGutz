import unittest, tempfile, os, shutil, glob, h5py, numpy, subprocess


def compare_data(path):
    with h5py.File('GLOG.h5', 'r') as f:
        data1 = f[path][()]
    with h5py.File('GLOG_REF.h5', 'r') as f:
        data2 = f[path][()]
    return numpy.allclose(data1, data2)


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
