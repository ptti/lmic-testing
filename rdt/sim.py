import argparse
import kappy
import numpy as np
import pandas as pd
from pyabc import Distribution, RV, ABCSMC

class Model(object):
    __name__ = "Kappa RDT"
    def __init__(self, model_files, fixed={}, stepsize=1, tmax=365):
        self.model_files = model_files
        self.stepsize = stepsize
        self.tmax = tmax
        self.fixed = fixed
        self.fixed_vars = "\n".join("%var: {} {}".format(k,v) for k,v in fixed.items())

    def __call__(self, params):
        ## assemble the model
        kasim = kappy.KappaStd()
        ## variable declarations for those that are fixed for this simulation
        kasim.add_model_string(self.fixed_vars, 1, file_id="fixed")
        ## variable declarations for those that we are trying to fit
        fit_vars = "\n".join("%var: {} {}".format(k,v) for k,v in params.items())
        kasim.add_model_string(fit_vars, 2, file_id="fit")
        ## model files
        for i, f in enumerate(self.model_files):
            kasim.add_model_file(f, i+3)
        kasim.project_parse()

        ## conduct the simulation
        meta = kappy.SimulationParameter(self.stepsize, "[T] > {}".format(self.tmax))
        kasim.simulation_start(meta)
        kasim.wait_for_simulation_stop()
        plot = kasim.simulation_plot()
        ## make the result usable
        df = pd.DataFrame(np.array(plot["series"])[::-1,:], columns=plot["legend"])

        ## compute R(t) according to S9.3 of
        ## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6002118
        N = self.fixed["N"]
        t = df["[T]"]
        Sn = df["Sn"]
        In = df["In"]
        I = df["In"] + df["Iy"]
        c = df["contact"]
        n = len(df)

        X = np.zeros(len(I))
        np.true_divide(Sn*In, I, out=X, where=I != 0)
        ker = np.exp(-0.14*t)
        bcs = params["beta"]*c*X

        R = []
        for i, tau in enumerate(t):
            s = np.pad(bcs, (n-i-1,0), mode="edge")
            R.append(np.trapz(s[:n]*ker[::-1]/N, t))
        df["R"] = np.array(R)

        ## done
        return df

def distance_target_R(x, y):
    return np.absolute(x["R"] - y["R"]).sum()

def command():
    parser = argparse.ArgumentParser()
    parser.add_argument("model", nargs="+", help="Kappa model files to simulate")
    parser.add_argument("-fit", nargs="*", default="", help="Parameters to fit")
    parser.add_argument("-fix", nargs="*", default="", help="Parameters to fix")
    parser.add_argument("-N", type=int, default="10000", help="Population size")
    parser.add_argument("-I", type=int, default="10", help="Initial infected")
    parser.add_argument("-tmax", default=365, type=int, help="Simulation max time")
    parser.add_argument("-db", default="sqlite:///abc.db", help="Database for ABC MCMC results")
    parser.add_argument("-R", default=1.0, type=float, help="Target R(t)")
    args = parser.parse_args()

    fixed = dict((k,float(v)) for k,v in map(lambda s: s.split("="), args.fix))
    fixed["N"] = args.N
    fixed["INIT_I"] = args.I
    m = Model(args.model, fixed=fixed, tmax=args.tmax)

    priors = dict((n, RV("uniform", float(lb),float(ub)))
                  for (n, lb, ub)
                  in map(lambda v: v.split(":"), args.fit))
    prior = Distribution(priors)

    abc = ABCSMC(m, prior, distance_target_R)
    abc_id = abc.new(args.db, { "R": args.R })
    history = abc.run(max_nr_populations=15)

    df, w = history.get_distribution()
    best = np.argmax(w)
    print(df.iloc[best])

if __name__ == '__main__':
    command()
