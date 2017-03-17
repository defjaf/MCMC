import cPickle
import numpy
import pylab

def test(write=True,asarray=False):
    
    a = numpy.linspace(-3,3,num=100)
    
    if write:
        f1 = file("a.cpkl", 'w')
        cPickle.dump(a, f1)
        f1.close()
        
    f1 = open("a.cpkl", 'r')
    a1 = cPickle.load(f1)
    f1.close()
            
    pylab.subplot(1,2,1)
    h = pylab.hist(a)

    if asarray:
        a1 = numpy.asarray(a1, dtype=numpy.float64)
        
    pylab.subplot(1,2,2)
    h1 = pylab.hist(a1)
    
    return a, a1

