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
        group.add_argument("--optimize", help="Optimization target", choices=["CBF", "ATT", "both"], default="both")
        group.add_argument("--debug", help="Debug mode", action="store_true")
        
        group = self.add_argument_group("ATT distribution")
        group.add_argument("--att-start", help="Starting value for ATT distribution (s)", type=float, default=0.2)
        group.add_argument("--att-end", help="Ending value for ATT distribution (s)", type=float, default=2.1)
        group.add_argument("--att-step", help="Step value for ATT distribution (s)", type=float, default=0.001)
        group.add_argument("--att-taper", help="Length of taper for ATT distribution (s)", type=float, default=0.3)

        group = self.add_argument_group("Scan to optimize for")
        group.add_argument("--asltype", help="ASL data type", choices=["var_multi_pCASL", "look_locker"], default="var_multi_pCASL")
        group.add_argument("--scan-duration", help="Desired scan duration (s)", type=float, default=300)
        group.add_argument("--scan-npld", help="Number of PLDs", type=int, default=6)
        group.add_argument("--scan-readout", help="Scan readout (non-ASL) time (s)", type=float, default=0.5)
        group.add_argument("--scan-ld", "--scan-tau", "--scan-bolusdur", help="Labelling duration (s)", type=float, default=1.4)
        group.add_argument("--scan-noise", help="Additive noise std.dev. relative to M0", type=float, default=0.002)
        group.add_argument("-f", help="CBF value to optimize for", type=float, default=50.0/6000)
        
        group = self.add_argument_group("PLD limits")
        group.add_argument("--pld-min", help="Minimum PLD (s)", type=float, default=0.1)
        group.add_argument("--pld-max", help="Maximum PLD (s)", type=float, default=3.0)
        group.add_argument("--pld-step", help="Step to search for optimal PLDs (s)", type=float, default=0.025)        

        group = self.add_argument_group("Labelling duration limits")
        group.add_argument("--optimize-ld", help="Optimize over varying labelling durations", action="store_true", default=False)
        group.add_argument("--ld-min", help="Minimum LD (s)", type=float, default=0.1)
        group.add_argument("--ld-max", help="Maximum LD (s)", type=float, default=1.8)
        group.add_argument("--ld-step", help="Step to search for optimal LD (s)", type=float, default=0.025)        

def main():
    try:
        arg_parser = OptPcaslArgumentParser()
        options = arg_parser.parse_args()

        welcome = "OXASL - PCASL Optimizer %s" % opt.__version__
        print(welcome)
        print("=" * len(welcome))
        
        # Define the ASL parameters
        params = opt.ASLParams(**vars(options))
        
        # ATT (BAT) distribution
        att_dist = opt.ATTDist(options.att_start, options.att_end, options.att_step, options.att_taper)

        # Details of the desired scan to optimize for
        scan = opt.ASLScan(options.asltype, duration=options.scan_duration, npld=options.scan_npld, 
                           readout=options.scan_readout, ld=options.scan_ld, noise=options.scan_noise)
        
        # PLD limits and step size to search over
        pld_lims = opt.Limits(options.pld_min, options.pld_max, options.pld_step, name="PLD")
        
        # LD limits and step size to search over
        if options.optimize_ld:
            ld_lims = opt.Limits(options.ld_min, options.ld_max, options.ld_step, name="LD")
        else:
            ld_lims = None

        # Type of optimisation
        # Note: the output best_min_variance is not comparable between D-optimal and L-optimal
        if options.optimize == "CBF":
            optimizer = opt.LOptimal([[1, 0],  [0, 0]], params, scan, att_dist, pld_lims)
        elif options.optimize == "ATT":
            optimizer = opt.LOptimal([[0, 0],  [0, 1]], params, scan, att_dist, pld_lims)
        else:
            optimizer = opt.DOptimal(params, scan, att_dist, pld_lims)

        # Run the optimisation
        output = optimizer.optimize()
        print("Optimal PLDs: %s" % output.plds)
        print("Number of repeats: %i" % output.num_av)
        print("Scan time: %.1fs" % output.scan_time)
        print("DONE")
    except (RuntimeError, ValueError) as exc:
        sys.stderr.write("ERROR: %s\n" % str(exc))
        if "--debug" in sys.argv:
            import traceback
            traceback.print_exc()
