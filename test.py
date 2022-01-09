import os
import sys
import imageio
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
        "        MODE : [debug,release,linux]\n\n"
        "    python test.py run <MODE> <TEST_CASE> <NUM_CELLS> <SAVE_INT>\n\n" +
        "        MODE      : [debug,release,linux]\n" +
        "        TEST_CASE : [1,2,3,4,5,6]\n" +
        "        SAVE_INT  : interval in seconds that the solution data are saved\n\n"
        "    Available test cases:\n" +
        "        1. Wet dam break\n" +
        "        2. Dry dam break\n" +
        "        3. Dry dam break with friction\n" +
        "        4. Wet c property\n" +
        "        5. Wet dry c property\n" +
        "        6. Building overtopping"
    )
    
    sys.exit(help_message)

####################################
# NATURAL SORTING, READ UP ON THIS #
####################################

import re

def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    '''
    return [ atoi(c) for c in re.split(r'(\d+)', text) ]

def get_filenames_natural_order():
    path = os.path.dirname( os.path.abspath(__file__) )
    
    filenames_natural_order = os.listdir(path)
    
    filenames_natural_order.sort(key=natural_keys)
    
    return filenames_natural_order

####################################
####################################
####################################

def set_path(mode):
    if mode == "debug":
        path = os.path.join(os.path.dirname( os.path.abspath(__file__) ), "out", "build", "x64-Debug")
    elif mode == "release":
        path = os.path.join(os.path.dirname( os.path.abspath(__file__) ), "out", "build", "x64-Release")
    elif mode == "linux":
        path = os.path.join(os.path.dirname( os.path.abspath(__file__) ), "build")
    else:
        EXIT_HELP()
        
    return path

def clear_jpg_files():
    path = os.path.dirname( os.path.abspath(__file__) )
    
    [ os.remove(filename) for filename in os.listdir(path) if filename.endswith(".jpg") ]

def clear_csv_files():
    path = os.path.dirname( os.path.abspath(__file__) )
    
    [ os.remove(filename) for filename in os.listdir(path) if filename.endswith(".csv") ]

test_names = [
    "wet-dam-break",
    "dry-dam-break",
    "dry-dam-break-fric",
    "wet-c-prop",
    "wet-dry-c-prop",
    "building-overtopping",
    "triange-dam-break"
]

sim_times = [
    2.5,
    1.3,
    1.3,
    0.5,
    0.5,
    10
]

class Limits:
    def __init__(
            self,
            intervals
        ):
            eta = []
            q   = []
            z   = []
            
            for interval in range(intervals):
                filename  = "solution_data-" + str(interval) + ".csv"
                
                dataframe = pd.read_csv(filename)
                
                eta += dataframe["eta"].tolist()
                q   += dataframe["q"].tolist()
                z   += dataframe["z"].tolist()
                
            self.etamax = np.max(eta)
            self.qmax   = np.max(q)
            self.zmax   = np.max(z)
            
            self.etamin = np.min(eta)
            self.qmin   = np.min(q)
            self.zmin   = np.min(z)

def plot_soln(
        x,
        y,
        topo,
        ylim,
        topolim,
        quantity,
        interval,
        ylabel,
        test_num,
        test_name
    ):
        filename = test_name + "-" + str(interval) + "-" + quantity + ".jpg"
        
        xlim = [ np.amin(x), np.amax(x) ]
        
        topo_scale_factor = ylim[1] / topolim[1]
        topo_scaled       = ( topo * topo_scale_factor + ylim[0] ) / 2
        
        plt.fill_between(x, y,           color=(65 / 255, 136 / 255, 239 / 255))    
        plt.fill_between(x, topo_scaled, color=(159 / 255, 159 / 255, 159 / 255))    
        
        plt.plot(x, y          , color='k', linewidth=0.5)    
        plt.plot(x, topo_scaled, color='k', linewidth=0.5)    
            
        plt.xlim(xlim)
        plt.ylim( 0.98 * ylim[0], 1.02 * ylim[1] )
        plt.xlabel("$x \, (m)$")
        plt.ylabel(ylabel)
        plt.savefig(filename, bbox_inches="tight")
        plt.clf()

class Solution:
    def __init__(
        self,
        mode,
        interval,
        test_num=0,
        test_name="ad-hoc"
    ):
        print("Initialising solution...")
        
        if mode != "debug" and mode != "release" and mode != "linux":
            EXIT_HELP()
        
        self.interval  = interval
        self.test_num  = test_num
        self.test_name = test_name
        
        filename = "solution_data-" + str(self.interval) + ".csv"
        
        dataframe = pd.read_csv(filename)
        
        self.x   = dataframe["x"].values
        self.q   = dataframe["q"].values
        self.z   = dataframe["z"].values
        self.eta = dataframe["eta"].values
        
    def plot_soln(
        self,
        limits
    ):
        print("Plotting solution interval " + str(self.interval) + "...")
        
        qlim   = (limits.qmin,   limits.qmax)
        etalim = (limits.etamin, limits.etamax)
        zlim   = (limits.zmin,   limits.zmax)
        
        plot_soln(self.x, self.q,   self.z, qlim,   zlim,   "q", self.interval, "$q \, (m^2s^{-1})$", self.test_num, self.test_name)
        plot_soln(self.x, self.eta, self.z, etalim, zlim, "eta", self.interval, "$\eta \, (m)$",      self.test_num, self.test_name)
        plot_soln(self.x, self.z,   self.z, zlim,   zlim,   "z", self.interval, "$z \, (m)$",         self.test_num, self.test_name)

def animate():
    images = []
    
    filenames_natural_order = get_filenames_natural_order()
    
    vars = ["eta", "q"]
    
    for test_name in test_names:
        for var in vars:
            suffix = var + ".jpg"
            for filename in filenames_natural_order:
                if filename.startswith(test_name) and filename.endswith(suffix):
                    image = imageio.imread(filename)
                    
                    images.append(image)
            
            if images:
                imageio.mimsave(test_name + "-" + var + ".gif", images)
                images = []

def run():
    print("Attempting to run...")
    
    if len(sys.argv) > 5:
        dummy, action, mode, test, num_cells, saveint = sys.argv
        
        path = set_path(mode)
        
        print("Trying to find executable in " + path)
        
        if mode == "debug" or mode == "release" or mode == "linux":
            solver_file = os.path.join(path, "fv1-cpp" if mode == "linux" else "fv1-cpp.exe")
        else:
            EXIT_HELP()
        
        subprocess.run( [solver_file, str(test), num_cells, saveint] )
        
        intervals = int( sim_times[int(test) - 1] / float(saveint) )
        
        [ Solution(mode, interval, test, test_names[int(test) - 1]).plot_soln( Limits(intervals) ) for interval in range(intervals) ]
    else:
        EXIT_HELP()
        
def run_tests():
    print("Attempting to run all tests...")
    
    if len(sys.argv) > 3:
        dummy, action, mode, num_cells = sys.argv
        
        path = set_path(mode)
        
        print("Trying to find executable in " + path)
        
        if mode == "debug" or mode == "release" or mode == "linux":
            solver_file = os.path.join(path, "fv1-cpp" if mode == "linux" else "fv1-cpp.exe")
        else:
            EXIT_HELP()
        
        tests = [1, 2, 3, 4, 5, 6]
        
        for i, test in enumerate(tests):
            subprocess.run( [ solver_file, str(test), num_cells, str( sim_times[i] ) ] )
            
            Solution( mode, 1, test, test_names[i] ).plot_soln( Limits(1) )
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
    clear_jpg_files()
    
    action = sys.argv[1]
    
    if action == "run":
        run()
        animate()
        clear_jpg_files()
        clear_csv_files()
    elif action == "test":
        run_tests()
    else:
        EXIT_HELP()
    
else:
    EXIT_HELP()