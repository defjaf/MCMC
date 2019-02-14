
from .. import Likelihood

class simpleLikelihood(Likelihood.Likelihood):


    def lnLike1(self, model_vals, data):
        return 0.0 #fill this in???

    def lnNorm1():
        return 0.0
    

    def setModel(self, params):
        """set the model for the given parameters"""
        self.this_model = self.model(params)
        self.model_vals = self.this_model()

    def lnLike(self, params): 
        """
        calculate the sum of the lnLike for each of the datasets
        i.e., \sum_i ln p(data_i |params) up to a param-independent term,
        """
        self.setModel(params)
        return sum( [ self.lnLike1(self.model_vals, d) for d in self.data ] )
