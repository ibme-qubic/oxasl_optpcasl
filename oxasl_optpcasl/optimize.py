"""
OXASL_OPTPCASL - Functions for optimizing the PLDs of a multi-pld pCASL protocol

Copyright 2019 University of Oxford
"""
import sys
import copy
import itertools

import numpy as np

from ._version import __version__

class OptimizationOutput(object):
    """
    Stores the output of the optimization process as attributes
    """

    def __init__(self):
        self.cost_history = []
        self.param_history = []
        self.best_cost = 1e99

class Optimizer(object):
    """
    Optimization base class
    
    Defines the basic algorithm for optimization. The L-optimal and D-optimal subclasses
    implement the particular cost functions for these methods
    """

    def __init__(self, scantype, log=sys.stdout):
        """
        :param params: ASLParams instance
        :param scan: ASLScan instance giving details of scan to optimize for
        :param att_dist: ATTDist instance giving ATT time distribution and weighting
        :param pld_lims: Limits instance, giving pld limiting values (max/min)
        :param log: Stream-like object for logging output
        """
        self.log = log
        self.scantype = scantype

    def gridsearch(self, gridpts=1e6):
        self.log.write("Doing initial grid search for best parameters\n")
        param_bounds = self.scantype.param_bounds()
        nparams = len(param_bounds)
        grid_npts =  max(3, int(gridpts**(1.0/nparams)))
        self.log.write(" - %i parameters, %i grid points\n" % (nparams, grid_npts))
        
        
        best_cost = 1e99
        best_params = []
        tol = 0.01
        it=1
        while 1:
            param_ranges = [np.linspace(bounds[0], bounds[1], grid_npts) for bounds in param_bounds]
            param_space = itertools.product(*param_ranges)
            #print(param_ranges)

            batch_size = 1000
            finish = True
            
            n_batches = int(np.ceil(grid_npts ** self.scantype.scan_params.npld / batch_size))
            for batch in range(n_batches):
            #for idx, params in enumerate(param_space):
                #self.log.write(" - batch %i/%i %f %s\n" % (batch+1, n_batches, best_cost, best_params))
                #try:
                batch = [next(param_space, None) for idx in range(batch_size)]
                #except StopIteration:
                #    IndexError:
                #    batch = param_space[batch*batch_size:, ...]
                batch = np.array([t for t in batch if t is not None])
                #if idx % 1000 == 0: print(idx, best_cost, best_params)
                #print(batch.shape)
                cost = self.scantype.cost(batch)
                #cost = self.scantype.cost(np.array(params))
                min_cost = np.min(cost)
                if min_cost < best_cost:
                    diff = abs(best_cost - min_cost)
                    best_params = batch[np.argmin(cost), :]
                    #best_params = np.array(params)
                    best_cost = min_cost
                    if diff > tol:
                        finish = False

            self.log.write(" - Iteration %i: best cost: %f params: %s\n" % (it, best_cost, best_params))
            if finish: break

            it += 1
            new_param_bounds = []
            for param, bounds in zip(best_params, param_bounds):
                width = (bounds[1] - bounds[0]) / (grid_npts-1)
                #print(param, bounds, width)
                new_bounds = [param - width, param + width]
                new_bounds[0] = max(new_bounds[0], bounds[0])
                new_bounds[1] = min(new_bounds[1], bounds[1])
                new_param_bounds.append(new_bounds)
                #print("new boudns", new_bounds)
            param_bounds = new_param_bounds

        best_params = sorted(best_params)
        self.log.write(" - Initializing with parameters: %s\n" % best_params)
        self.log.write("DONE")
        return np.array(best_params)

    def optimize(self, initial_params=None, reps=1):
        self.log.write("Optimizing PLDs for: %s\n" % self.scantype)
        self.log.write("PLD search limits: %s\n" % self.scantype.pld_lims)
        self.log.write("Optimizing for %i PLDs\n" % self.scantype.scan_params.npld)
        self.log.write("%s\n" % str(self.scantype.att_dist))
        #self.log.write("Optimization method: %s\n" % self.name)
        self.log.write("\n")

        best_output, best_cost = None, 1e99
        for rep in range(reps):
            self.log.write("Optimization %i/%i... " % (rep+1, reps))
            output = self._optimize_once(initial_params)
            self.log.write("DONE - Optimized parameters: %s (cost: %f)\n" % (output["params"], output["best_cost"]))
            if output["best_cost"] < best_cost:
                    best_cost = output["best_cost"]
                    best_output = output

        self.log.write("Final parameters: %s (cost: %f)\n" % (best_output["params"], best_output["best_cost"]))
        return best_output

    def _optimize_once(self, initial_params, rand_order=True):
        iters = 0
        current_params = np.copy(initial_params)
        cost_history, param_history = [], []
        while 1:
            old_params = np.copy(current_params)

            param_indices = range(len(current_params))
            if rand_order:
                param_indices = np.random.permutation(param_indices)

            for idx in param_indices:
                #print("Optimizing parameter %i" % idx)
                trial_params = self.scantype.trial_params(current_params, idx)
                #print("Trials: %s" % trial_params[:, idx])
                cost = self.scantype.cost(trial_params)
                #print("Cost: %s" % cost)
                if len(cost) == 0:
                    continue
                 
                #print(cost)
                #import sys
                #sys.exit(1)
                min_cost = np.min(cost)
                min_cost_idx = np.argmin(cost)
                current_params = trial_params[min_cost_idx]
                #print("New PLDs: %s (cost: %f)" % (current_params, min_cost))
            
                # Save the results from each step for visualising
                cost_history.append(min_cost)
                param_history.append(current_params)

                iters += 1

            if np.allclose(current_params, old_params):
                break
        

        num_av, total_tr = self.scantype.repeats_total_tr(current_params)
        output = {
            "params" : np.array(current_params),
            "best_cost" : min_cost,
            "cost_history" : cost_history,
            "param_history" : param_history,
            "num_iters" : len(cost_history),
            "num_av" : num_av,
            "total_tr" : total_tr,
            "scan_time": total_tr * num_av,
        }
        output.update(self.scantype.name_params(current_params))
        return output
