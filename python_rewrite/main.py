# \brief Handles running the analysis of the supplied input data as
#        well as applying it to actual data sets
# 
# \author Tristan Protzman
# \date 01/10/2021
# \email tlprotzman@gmail.com
# \affiliation Lehigh University
# 
#

import sys

import training



def main(args):
    # Training Phase
    model = training.centrality_model(True)
    model.Import_data("/home/tristan/rhic/epd_centrality/data/simulated_data.root")
    model.build_model()

    # Application Phase
    model.apply_model()
    model.evaluate_model()



if __name__ == "__main__":
    main(sys.argv)