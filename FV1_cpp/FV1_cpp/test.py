import os
import sys
import subprocess
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.pylab  as pylab

from mpl_toolkits.mplot3d import Axes3D

def EXIT_HELP():
    help_message = (
        "Use as:\n\n" +
        "    python test.py test <MODE> <NUM_CELLS>\n\n" +
        "        MODE : [debug,release]\n\n"
        "    python test.py run <MODE> <TEST_CASE> <NUM_CELLS>\n\n" +
        "        MODE      : [debug,release]\n" +
        "        TEST_CASE : [1,2,3,4,5,6]\n\n" +
        "    Available test cases:\n" +
        "        1. Wet dam break\n" +
        "        2. Dry dam break\n" +
        "        3. Dry dam break with friction\n" +
        "        4. Wet c property\n" +
        "        5. Wet dry c property\n" +
        "        6. Building overtopping"
    )
    
    sys.exit(help_message)

test_names = [
    "wet-dam-break",
    "dry-dam-break",
    "dry-dam-break-fric",
    "wet-c-prop",
    "wet-dry-c-prop",
    "building-overtopping",
]

def plot_soln(
        x,
        y,
        quantity,
        ylabel,
        test_num,
        test_name
    ):
        filename = str(test_num) + "-" + test_name + "-" + quantity
        
        xlim = [ np.amin(x), np.amax(x) ]
        
        ylim = [np.amin(y) * 0.98, np.amax(y) * 1.02]
        
        plt.plot(x, y, color='b', linewidth="1.5")    
            
        plt.xlim(xlim)
        plt.ylim(ylim)
        plt.xlabel("$x \, (m)$")
        plt.ylabel(ylabel)
        plt.savefig(filename, bbox_inches="tight")
        plt.clf()

class Solution:
    def __init__(
        self,
        mode,
        test_num=0,
        test_name="ad-hoc"
    ):
        print("Initialising solution...")
        
        if mode != "debug" and mode != "release":
            EXIT_HELP()
        
        self.test_num  = test_num
        self.test_name = test_name
        
        dataframe = pd.read_csv("solution_data.csv")
        
        self.x   = dataframe["x"].values
        self.q   = dataframe["q"].values
        self.z   = dataframe["z"].values
        self.eta = dataframe["eta"].values
        
    def plot_soln(
        self
    ):
        print("Plotting solution...")
        
        plot_soln(self.x, self.q,   "q",   "$q \, (m^2s^{-1})$", self.test_num, self.test_name)
        plot_soln(self.x, self.eta, "eta", "$\eta \, (m)$", self.test_num, self.test_name)
        plot_soln(self.x, self.z,   "z",   "$z \, (m)$", self.test_num, self.test_name)
    
def run():
    print("Attempting to run...")
    
    if len(sys.argv) > 4:
        dummy, action, mode, test, num_cells = sys.argv
        
        if   mode == "debug":
            solver_file = os.path.join(os.path.dirname(__file__), "..", "x64", "Debug", "FV1_cpp.exe")
        elif mode == "release":
            solver_file = os.path.join(os.path.dirname(__file__), "..", "x64", "Release", "FV1_cpp.exe")
        else:
            EXIT_HELP()
        
        subprocess.run( [solver_file, str(test), num_cells] )
        Solution(mode, test, test_names[int(test) - 1]).plot_soln()
    else:
        EXIT_HELP()
        
def run_tests():
    print("Attempting to run all tests...")
    
    if len(sys.argv) > 3:
        dummy, action, mode, num_cells = sys.argv
        
        if   mode == "debug":
            solver_file = os.path.join(os.path.dirname(__file__), "..", "x64", "Debug", "FV1_cpp.exe")
        elif mode == "release":
            solver_file = os.path.join(os.path.dirname(__file__), "..", "x64", "Release", "FV1_cpp.exe")
        else:
            EXIT_HELP()
        
        tests = [1, 2, 3, 4, 5, 6]
        
        for i, test in enumerate(tests):
            subprocess.run( [solver_file, str(test), num_cells] )
            Solution(mode, test, test_names[i]).plot_soln()
    else:
        EXIT_HELP()

params = {
    "legend.fontsize" : "xx-large",
    "axes.labelsize"  : "xx-large",
    "axes.titlesize"  : "xx-large",
    "xtick.labelsize" : "xx-large",
    "ytick.labelsize" : "xx-large"
}

pylab.rcParams.update(params)

if len(sys.argv) > 1:
    action = sys.argv[1]
    
    if action == "run":
        run()
    elif action == "test":
        run_tests()
    else:
        EXIT_HELP()
else:
    EXIT_HELP()