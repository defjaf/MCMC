from __future__ import division
import BeamData
import BeamModel
import Likelihood

from math import pi

def main():
    beam = BeamModel.GaussianBeamModel2D((0,0), (0.05, 0.05), pi/4)
    dat = BeamData.BeamSim(beam, 1000, 1.0, 5.0)
    like = Likelihood.Likelihood(data=dat, beam_model=BeamModel.GaussianBeamModel2D)

    return like
