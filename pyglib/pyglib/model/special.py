import numpy

'''Special models'''


class semicirular(object):
    '''Semi-circular DOS.
    '''
    def __init__(self):
        self.dos = lambda e: 2./numpy.pi * numpy.sqrt(1-e**2)
        self.cdos = lambda e: (e*numpy.sqrt(1-e**2) \
                + numpy.arcsin(e)) / numpy.pi + 0.5

    def get_e_list_of_uniform_wt(self, nmesh=5000):
        '''Get the energy mesh for semi-circular DOS with uniform weight.
        '''
        cdos_list = numpy.linspace(0,1,nmesh+1)
        from scipy.optimize import bisect
        e_list = [bisect(lambda x: self.cdos(x)-a, -1 ,1) \
                for a in cdos_list]
        e_list = numpy.asarray(e_list)
        e_list = (e_list[1:] + e_list[0:-1])/2
        return e_list



if __name__=='__main__':
    import matplotlib.pyplot as plt

    # Uniform e-mesh with different weight.
    nmesh = 500
    e_list = numpy.linspace(-1,1,nmesh)
    sc = semicirular()
    plt.plot(e_list,sc.dos(e_list))
    plt.plot(e_list,sc.cdos(e_list))
    plt.show()

    # Non-uniform e-mesh with uniform weight.
    nmesh = 50
    e_list = sc.get_e_list_of_uniform_wt(nmesh=nmesh)
    cdos_list = numpy.linspace(0,1,nmesh)
    plt.plot(e_list,sc.cdos(e_list),'o')
    plt.plot(e_list,cdos_list,'-')
    plt.show()
