"""
This main script generates a tissue network based on \nu
and performs a simulation. 
Create causal graphs for various \nu values (i.e change in cross connections which leads to fibrillations)

The output files are is saved to to ./output/
"""
import argparse
import multiprocessing
import os
import time
import numpy as np
import matplotlib.pyplot as plt
import CAAF.model as model
from CAAF.utils import generateCausalNetwork

parser = argparse.ArgumentParser()
parser.add_argument(
    "-n",
    "--num_iters",
    help="number of iterations to run the model",
    type=int,
    default=1000,
)
parser.add_argument(
    "-p",
    "--plot",
    help="saves animation of simulation in a .mp4 file",
    type=int,
    default=1,
)
args = parser.parse_args()

num_iters = args.num_iters
plot = args.plot

nuvals = [0.1,0.2,0.5,0.6,0.8,0.9]

def processnu(exptn,nu):   
    parameters = {
        "row_size": 200,
        "col_size": 200,
        "refractory_period": 50,
        "driving_period": 220,
        "prob_not_fire": 0.05,
        "prob_con": nu, # 0.09 less than 1, fibrillation after 0.14
        "prob_def": 0.05,
        "filenameprefix":f'./output/SubNx_{exptn}'
    }
    try:
        np.random.seed(0)
        start_time = time.time()
        heart = model.Heart(**parameters)
        heart.update(num_iters=num_iters, plot=plot)
        #Generate the Causal Network and specify sampling grid size and the alpha level (1e-3) for choosing the causal connections among grod nodes
        generateCausalNetwork(heart,f'{parameters["filenameprefix"]}.paj',10,10,1e-3)
        # Number of active cells with time
        fig, ax = plt.subplots(figsize=[5, 5])
        ax.plot(heart.num_active_cells)
        fig.savefig(f"./output/num_active_cells_{exptn}.png")
    except:
        return f"{exptn} failed nu value was {nu}"
    return f"{exptn} solved nu value was {nu}"


import multiprocessing
from multiprocessing import Pool
if __name__ == "__main__":
    
    os.makedirs('output', exist_ok=True)
    pool = multiprocessing.Pool(processes=len(nuvals))
    results = [pool.apply_async(processnu, args=(c,nuvals[c])) for c in range(len(nuvals))] # maps function to iterator
    output = [p.get() for p in results]        
