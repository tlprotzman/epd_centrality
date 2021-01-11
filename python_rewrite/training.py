#
# \brief training.py contains the functions to fit data from the
#        Event Plane Detector to either refmult or the impact
#        parameter, dependent on if the data is real or simulated
# 
# \author Tristan Protzman
# \date 
# \email tlprotzman@gmail.com
# \affiliation Lehigh University
# 
#


from os import error
import typing

import ROOT
import numpy as np
# import matplotlib.pyplot as plt

class centrality_model():
    
    def __init__(self, simulated_data : bool):
        self.simulation : bool
        self.simulation = simulated_data
        
    def Import_data(self, data_input : str) -> None:
        """Loads the data from the specified root file into memory,
           allowing the model to be created

        Args:
            data_path (str): The path to the root file containing the summary
                             of the ring values and the impact parameter/refmult
        """

        
        inFile = ROOT.TFile(data_input)
        
        ring_data = inFile.Get("ring_sums")    # Open and read in ring data from the ROOT file
        self.num_events = ring_data.GetNcols()
        if ring_data.GetNrows() != 16:
            raise ValueError("16 rows were not found importing training data")

        self.training_input = np.empty((self.num_events, 16))

        for i in range(16):     # TODO There must be a more efficient way to copy all the data into a numpy array
            for j in range(self.num_events):
                self.training_input[j][i] = ring_data[i][j]

        target_data = None
        if self.simulation:
            target_data = inFile.Get("impact_parameter")
        else:
            target_data = inFile.Get("tpc_multiplicity")

        if target_data.GetNrows() != self.num_events:
            raise ValueError("Number of events are different for ring data and b/refmult")


        

        inFile.Close()

    
    def build_model(self):
        pass

    def apply_model(self):
        pass

    def evaluate_model(self):
        pass
