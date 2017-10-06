import unittest, tempfile, os, shutil, glob, subprocess
from pyglib.basic.data import compare_array as compare_data


class KnowValues(unittest.TestCase):
    def test_dm_matrix(self):
        self.assertTrue(compare_data('/DM', fname1='EMBED_HAMIL_RES_1.h5', \
                fname2='EMBED_HAMIL_RES_1_REF.h5'))

    def test_etot_matrix(self):
        self.assertTrue(compare_data('/emol', fname1='EMBED_HAMIL_RES_1.h5', \
                fname2='EMBED_HAMIL_RES_1_REF.h5'))



if __name__ == "__main__":
    print("test exe_spci_s2_mott using feo example.")
    cwd = os.getcwd()
    tempd = tempfile.mkdtemp(prefix='gtest_tmp_')
    for f in glob.glob('*.h5'):
        shutil.copy(f, tempd)
    os.chdir(tempd)
    cmd = [os.environ['WIEN_GUTZ_ROOT2']+'/exe_spci_s2_mott']
    subprocess.call(cmd)
    unittest.main()
