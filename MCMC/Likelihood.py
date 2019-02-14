

from numpy import sum, log, exp, asarray
from math import pi

# to subclass, need to change
#    lnLike1(self, model_vals=None, data=None):
#    lnNorm1 [optional, rarely used]
#    setmodel_class, setModel (depending on how the model deals with parameters)

### AHJ April 2014
###  add method to calculate chi2 at a parameter value.

class Likelihood(object):
    """
    Represents a likelihood function for beam fitting.

    I.E., the specific case of additive gaussian errors, a
    parameterized model beamshape with the **amplitude as a nuisance
    parameter**
    
    initialize with a dataset and the name of a *class* (not instance)
    describing the model -- either in __init__ or the appropriate set*
    function.

    given those, the __call__ (=lnLike) operator takes the model
    parameters and returns the *parameter-dependent* terms of
    **logarithm** of the the likelihood value

    also calculate the paramater-independent terms of the ln-likelihood

    calculate the posterior (model must define lnPrior)
    problem: what if prior -> 0 for some params????

    now works for a sequence of datasets
     - store only sequence even if only one dataset
     - relies on a single not being a sequence
     
     snow also calculates the ML value of the amplitude parameter, 
     NOT DONE, BUT SEE CODE **and requires ML amplitude>0**
     
    """
    
    nDerived = 0
    
    def __init__(self, data=None, model=None):
        self.setdata(data)
        self.setmodel_class(model)

    def setmodel_class(self, model):
        self.model = model
        self.this_model = None
        self.model_vals = None

    def adddata(self, dataset):
        if self.data is None:
            self.setdata(data)
        elif isiterable(dataset):
            self.data.extend(dataset)
        else:
            self.data.append(dataset)

    def setdata(self, data):
        if isiterable(data):
            self.data = data
        else:
            self.data = [data]
        self.model_vals = None

    def setModel(self, params):
        """set the model for the given paramaters"""
        self.this_model = self.model(*params)
        self.model_vals = [self.this_model(d) for d in self.data]
        self.MLamplitude = None
        #self.model_vals = map(self.this_model,self.data)

    ## this doesn't need to be a class method, although can just use self.model_vals!
    def lnLike1(self, model_vals=None, data=None):
        """
        return ln p(data|params) up to a param-independent term,
        - for a single dataset
        - assuming params have already been set by self.setModel
        """
        ## should work for scalar or vector sig2
        ## nb. FNid = F^T N^-1 d, etc.

        FNiF = data.quadform(model_vals)     
        FNid = data.quadform(model_vals, data.d)
        if FNiF<=0:
            raise ZeroPosterior(FNiF)
        self.MLamplitude = FNid/FNiF
        
        # ## require positive amplitudes
        # if self.MLamplitude<0:
        #     raise ZeroPosterior(self.MLamplitude)
        
        return 0.5 * (FNid*self.MLamplitude - log(FNiF))
        
    def lnLike(self, params): 
        """
        calculate the sum of the lnLike for each of the datasets
        i.e., \sum_i ln p(data_i |params) up to a param-independent term,
        """
        self.setModel(params)
        return sum( [ self.lnLike1(v, d) for (v, d) in zip(self.model_vals, self.data) ] ) 

    def chi2_1(self, model_vals=None, data=None):
        """
        return chi2(data|params) 
        - for a single dataset
        - assuming params have already been set by self.setModel
        """

        return data.chi2(model_vals)

        
    def chi2(self, params): 
        """
        calculate the sum of the chi2_1 for each of the datasets
        i.e., \sum_i chi2(data_i |params) 
        """
        self.setModel(params)
        return sum( [ self.chi2_1(v, d) for (v, d) in zip(self.model_vals, self.data) ] ) 
        

    def __call__(self, params):
        return self.lnLike(params)
        
    def lnNorm1(self, data=None):
        """
        return the parameter-independent part of the ln-likelihood
        - for a single dataset
        """       
        
        norm = (data.n-1)*log(2*pi)
        norm2 = data.lnDetN + data.chi2()
    
        return -0.5*(norm + norm2)

    def lnNorm(self): 
        """
        return the parameter-independent part of the ln-likelihood
        """
        return sum(map(self.lnNorm1, self.data))

    def lnLikeFull(self, params):
        """
        return the full ln-likelihood, including param-independent parts
        """
        return self.lnLike(params) + self.lnNorm()

    def lnPost(self, params):
        """
        return the parameter-independent part of the ln-posterior
        """
        return self.lnLike(params) + self.model.lnPrior(params)

    def lnPostFull(self):
        """
        return the full ln-posterior (possibly up to constant terms in
        the prior)
        """
        return self.lnPost(params) + self.lnNorm()
    
# should check for mutable sequence instead?
def isiterable(obj): 
    try: iter(obj) 
    except TypeError: return False 
    else: return True

class ZeroPosterior(Exception):
    def __init__(self, *value):
        self.value = value
    def __str__(self):
        return repr(self.value)
