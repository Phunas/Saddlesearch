from ase.optimize import FIRE, LBFGS, LBFGSLineSearch
from ase.optimize.precon import PreconLBFGS
from ase.build import bulk
from atomistica import Tersoff
from ode import ODE12rOptimizer

a0 = bulk('C', cubic=True) * 5
a0.rattle(0.05)
del a0[0]


for OPT in [ODE12rOptimizer, FIRE, LBFGS, LBFGSLineSearch, PreconLBFGS]:
    a = a0.copy()
    a.set_calculator(Tersoff())
    opt = OPT(a)
    opt.run(fmax=1e-3)
    print()
