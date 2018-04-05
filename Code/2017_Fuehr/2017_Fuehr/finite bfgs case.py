import sovi_bib as sovi
import numpy as np
from fenics import *
from dolfin import *

#l-bfgs test fuer endlichdimensionalen fall ohne shape

def f(x): return (x[0]-20.)**4 + (x[1]+50.)**4 + (x[0]-20)**2*(x[1]+50)**3
def f_grad(x): return np.array([4*(x[0]-20.)**3 + 2*(x[0]-20)*(x[1]+50)**3, 4*(x[1]+50.)**3+ 3*(x[0]-20)**2 * (x[1]+50)**2])
def bfgs_step(memories):

    if isinstance(memories, sovi.bfgs_memory): pass
    else: raise SystemExit("bfgs_step benoetigt Objekt der  BFGS-Memory-Klasse als Input!")

    b     = memories.gradient[0]
    alpha = np.zeros(memories.length)

    if(memories.step_nr + 1 >= memories.length):
        # bei voll besetzter memory werden alle Eintraege verwendet

        for i in range(memories.length-1):
            # Vorwaertsschleife
            i         = i+1
            diff_grad = memories.gradient[i-1] - memories.gradient[i]
            alpha[i]              = sum(memories.deformation[i-1]*b) / sum(diff_grad*memories.deformation[i-1])
            b = b - float(alpha[i])*diff_grad

        # Reskalierung von q
        first_diff_grad = memories.gradient[0] - memories.gradient[1]
        gamma           = np.inner(first_diff_grad, memories.deformation[0]) / sum(first_diff_grad*first_diff_grad)
        b               = np.inner(gamma,b)

        for i in range(memories.length-1):
            # Rueckwaertsschleife
            i                     = i+1
            diff_grad = memories.gradient[-(i+1)] - memories.gradient[-i]
            beta                  = sum(diff_grad*b) / sum( diff_grad*memories.deformation[ -(i+1)])
            b         = b + (float(alpha[-i]) - beta)*memories.deformation[ -(i+1)]

    elif(memories.step_nr == 0):
            # der erste BFGS-Schritt ist ein Gradientenschritt
            return -b

    else:
        # bei voll besetzter memory werden alle Eintraege verwendet

        for i in range(memories.step_nr):
            # Vorwaertsschleife
            i         = i+1
            diff_grad = memories.gradient[i-1] - memories.gradient[i]
            alpha[i]              = sum(memories.deformation[i-1]*b) / sum(diff_grad*memories.deformation[i-1])
            b = b - float(alpha[i])*diff_grad

        # Reskalierung von q
        first_diff_grad = memories.gradient[0] - memories.gradient[1]
        gamma           = np.inner(first_diff_grad, memories.deformation[0]) / sum(first_diff_grad*first_diff_grad)
        b               = np.inner(gamma,b)

        for i in range(memories.step_nr):
            # Rueckwaertsschleife
            shift = (memories.length-1)-memories.step_nr
            i                     = i+1
            diff_grad = memories.gradient[-(i+1)-shift] - memories.gradient[-i-shift]
            beta                  = sum(diff_grad*b) / sum( diff_grad*memories.deformation[ -(i+1)-shift])
            b         = b + (float(alpha[-i-shift]) - beta)*memories.deformation[ -(i+1)-shift]

    return -b

pointinspace = np.zeros(2)
length = 8
memo = sovi.bfgs_memory(np.zeros([length,2]), np.zeros([length,2]), length, 0)
for i in range(40):
    memo.update_grad(f_grad(pointinspace))
    upd = bfgs_step(memo)
    memo.update_defo(upd)
    memo.step_nr = memo.step_nr+1
    pointinspace = pointinspace + upd
    print(pointinspace)




