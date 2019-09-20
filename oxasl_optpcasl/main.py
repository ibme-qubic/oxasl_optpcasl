"""
OXASL_OPTPCASL - Entry point for command line program

Copyright 2019 University of Oxford
"""
import sys
import argparse

import oxasl_optpcasl as opt

USAGE = """oxasl_optpcasl <options>"""

class OptPcaslArgumentParser(argparse.ArgumentParser):
    """
    ArgumentParser for OPT_PCASL options
    """

    def __init__(self, **kwargs):
        argparse.ArgumentParser.__init__(self, prog="oxasl_optpcasl", usage=USAGE, **kwargs)

        group = self.add_argument_group("Main Options")
        group.add_argument("--asltype", help="ASL data type", required=True,
                           choices=["var_multi_pCASL", "var_te_pCASL", "look_locker", "var_te_pCASL_nPLD"])
        group.add_argument("-f", help="", type=float, default=50.0/6000)
        group.add_argument("--opt-type", help="Optimization type - L-optimal or D-optimal", choices=["L", "D"], default="D")
        group.add_argument("--debug", help="Debug mode", action="store_true")
        
        group = self.add_argument_group("ATT distribution")
        group.add_argument("--att-start", help="Starting value for ATT distribution (s)", type=float, default=0.2)
        group.add_argument("--att-end", help="Ending value for ATT distribution (s)", type=float, default=2.1)
        group.add_argument("--att-step", help="Step value for ATT distribution (s)", type=float, default=0.001)
        group.add_argument("--att-taper", help="Length of taper for ATT distribution (s)", type=float, default=0.3)

        group = self.add_argument_group("Scan to optimize for")
        group.add_argument("--scan-duration", help="Desired scan duration (s)", type=float, default=300)
        group.add_argument("--scan-npld", help="Number of PLDs", type=int, default=6)
        group.add_argument("--scan-readout", help="Scan readout time", type=float, default=0.5)
        
        group = self.add_argument_group("PLD limits")
        group.add_argument("--pld-min", help="Minimum PLD (s)", type=float, default=0.1)
        group.add_argument("--pld-max", help="Maximum PLD (s)", type=float, default=3.0)
        group.add_argument("--pld-step", help="Step to search for optimal PLDs (s)", type=float, default=0.025)        

def main():
    try:
        arg_parser = OptPcaslArgumentParser()
        options = arg_parser.parse_args()

        welcome = "OXASL - PCASL Optimizer %s" % opt.__version__
        print(welcome)
        print("=" * len(welcome))
        
        # Define the ASL parameters
        params = opt.ASLParams(options.asltype, options.f)
        
        # ATT (BAT) distribution
        att_dist = opt.BATDist(options.att_start, options.att_end, options.att_step, options.att_taper)

        # Details of the desired scan to optimize for
        scan = opt.Scan(duration=options.scan_duration, npld=options.scan_npld, readout=options.scan_readout)
        
        # PLD limits and step size to search over
        lims = opt.Limits(options.pld_min, options.pld_max, options.pld_step)
        
        # Type of optimisation
        # Note: the output best_min_variance is not comparable between D-optimal and L-optimal
        if options.opt_type == "L":
            opttype = opt.LOptimal([[1, 0],  [0, 0]])
        else:
            opttype = opt.DOptimal()

        # Run the optimisation
        best_plds, num_av, best_min_variance = opt.optimize(opttype, params, att_dist, scan, lims)

        print("DONE")
    except (RuntimeError, ValueError) as exc:
        sys.stderr.write("ERROR: %s\n" % str(exc))
        if "--debug" in sys.argv:
            import traceback
            traceback.print_exc()
