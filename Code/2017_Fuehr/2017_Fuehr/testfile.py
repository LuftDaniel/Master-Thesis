import sovi_bib as sovi
import numpy as np
from fenics import *
from dolfin import *


grad = [1.0]
defo = [3.0, 2.0, 1.0]

zero = np.zeros([3,3])

mem = sovi.bfgs_memory(zero, zero, 3, 1)
shit = range(100)

#stuff i used
"""
print(MeshData.mesh.num_vertices())
testerer = Function(VectorFunctionSpace(MeshData.mesh,"P",1, dim=2))
testerer.interpolate(Expression(("0.0", "0.0"), degree=1))
testerer.vector()[:] = np.zeros(2*3080)
print(testerer.vector().array())
print(sovi.bilin_a(MeshData, testerer, U, 2.0))
"""

# Test der Memoryklasse
"""
print(MeshData.mesh.num_vertices())
testerer = Function(VectorFunctionSpace(MeshData.mesh,"P",1, dim=2))
testerer.interpolate(Expression(("0.0", "0.0"), degree=1))
testerer.vector()[:] = np.zeros(2*3080)
tester2 = U-U
testclass = sovi.bfgs_memory(np.zeros([4,2*3080]), np.zeros([4,2*3080]), 4, 0)
testclass.update_defo(U.vector().array())
print(U.vector().array() == testclass.deformation[0])
testclass.update_defo(np.zeros(2*3080))
testfunc = testclass.initialize_defo(MeshData.mesh, 1)
print(sovi.bilin_a(MeshData, testfunc, U, 2.0))
print(sovi.bilin_a(MeshData, U, U, 2.0))
testfunc = testclass.initialize_defo(MeshData.mesh, 0)
print(sovi.bilin_a(MeshData, testfunc, U, 2.0))
"""


# test der BFGS_step Funktion
"""
mesh = RectangleMesh(Point(0,0), Point(1, 1), 10, 10, "right")
positions = mesh.coordinates()

V = VectorFunctionSpace(mesh, "P", 1, dim=2)
print(mesh.num_vertices())

u = Function(V)
v = Expression(("x[0]","x[1]"), degree=1)
ar = np.zeros(2*mesh.num_vertices())
for i in range(len(ar)): ar[i] = float(i)

u.vector()[:] = ar
print(2*u.vector().array())
grill = Function(V)
grill.vector()[:] = 2*u.vector()
print(grill.vector().array())
dx = Measure("dx", mesh)
a = inner(grill,u)*dx
b = inner(2*u,u)*dx
print(assemble(a))
print(assemble(b))

grill.vector()[[2]] = 0.
print(grill.vector().array())
"""







