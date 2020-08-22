#	Compute the function R(t) for the reproduction number according to the
#       provided time-series for the susceptible population. Taken from S9.3 of
#       https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6002118/

import numpy as np
import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-b", "--beta", type=float, default=0.036)
parser.add_argument("-g", "--gamma", type=float, default=7**-1)
parser.add_argument("-i", "--input", type=str, required=True)
parser.add_argument("-o", "--output", type=str, required=True)

args = parser.parse_args()

beta = args.beta

traj = pd.read_csv(args.input, delimiter="\t")

N = np.sum([traj.iloc[0][i] for i in [1,2,3,4,5,6,7,8]])

t  = traj["t"]
Sn = traj["Sn"]
In = traj["In"]
I  = traj["In"] + traj["Iy"]
c  = traj["c"]
n = len(t)

X = np.zeros(len(I))
np.true_divide(Sn*In, I, out=X, where=I != 0)

ker = np.exp(-args.gamma*t)
bcs = beta*c*X

Rs = []
for i, tau in enumerate(t):
  s = np.pad(bcs, (n-i-1, 0), mode="edge")
  Rs.append(np.trapz(s[:n]*ker[::-1]/N, t))

head = "t\tR"
np.savetxt(args.output, np.vstack([t, np.array(Rs)]).T, delimiter="\t", encoding="utf-8", comments="", header=head)
