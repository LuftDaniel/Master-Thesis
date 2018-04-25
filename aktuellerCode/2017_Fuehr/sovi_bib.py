from fenics import *
import os
import errno
from datetime import datetime
from copy import copy
import numpy as np

__this_files_path = os.path.realpath(__file__)
__this_files_dir  = os.path.dirname(__this_files_path)

# Pfad zum Ordner mit den Gitter Dateien
DATA_DIR = os.path.abspath(os.path.join(__this_files_dir, 'Meshes'))

def create_outputfolder():

    # Erstellen eines Output Ordners, falls nicht vorhanden
    try:
        os.mkdir('Output')
    except OSError as exc:
        if exc.errno != errno.EEXIST:
            raise exc
        pass

    # Erstelle Ordner fuer jeden Durchlauf nach Datum und Zeit
    outputfolder = os.path.join(__this_files_dir,
                                'Output',
                                datetime.now().strftime('%Y%m%d_%H%M%S'))
    os.mkdir(outputfolder)

    # Gibt den Pfad zum speziellen Order im Output Ordner zurueck
    return outputfolder


class MeshData:

    # Objekt mit allen Daten des Gitters
    def __init__(self, mesh, subdomains, boundaries, ind):

        # FEniCS Mesh
        self.mesh = mesh

        # FEniCS Subdomains
        self.subdomains = subdomains

        # FEniCS Boundaries
        self.boundaries = boundaries

        # Indizes der Knotenpunkte mit Traeger nicht am inneren Rand
        self.indNotIntBoundary = ind


class bfgs_memory:

    """
    Schreibe zwei Funktionen, die Gradienten und Deformationen bei verschobenen Meshes
    neu als Fenics Functions initialisieren, indem sie Meshes nehmen und die Nummer
    des Gradienten/Defo.
    """

    # Objekt mit gespeicherten Gradienten und Deformationen der letzten l Schritte in Form von Arrays
    def __init__(self, gradient, deformation, length, step_nr):

        # Liste von Gradientenvektorfeldern
        if (len(gradient) == length): self.gradient = gradient
        else: raise SystemExit("Fehler: Anzahl der Gradienten passt nicht zur Memorylaenge!")

        # Liste von Deformationsvektorfeldern
        if (len(deformation) == length): self.deformation = deformation
        else: raise SystemExit("Fehler: Anzahl der Deformationen passt nicht zur Memorylaenge!")

        # Anzahl der gespeicherten letzten Schritte
        self.length = length

        # Anzahl der bereits ausgefuehrten l-BFGS-Schritte
        if (step_nr >= 0 and isinstance(step_nr, int)): self.step_nr = step_nr
        else: raise SystemExit("Fehler: step_nr muss Integer groesser gleich 0 sein!")

    # macht ein Update der Memory; neueste Elemente bekommen Index 0
    def update_grad(self, upd_grad):

        for i in range(self.length-1): self.gradient[-(i+1)]    = self.gradient[-(i+2)]

        self.gradient[0] = upd_grad

    def update_defo(self, upd_defo):

        for i in range(self.length-1): self.deformation[-(i+1)] = self.deformation[-(i+2)]

        self.deformation[0] = upd_defo

    def initialize_grad(self, meshData, i):
        # erzeugt eine FEniCS-Funktion des Gradientenfeldes i auf einem Mesh aus den gespeicherten Arrays
        # ermoeglicht dadurch Transport; i entspricht Index in Speichermatrix self.gradient, aufsteigend im Alter

        if isinstance(meshData, MeshData): pass
        else: raise SystemExit("initialize_grad benoetigt Objekt der MeshData-Klasse als Input!")

        V = VectorFunctionSpace(meshData.mesh, "P", 1, dim=2)
        f = Function(V)
        f.vector()[:] = self.gradient[i]
        return f

    def initialize_defo(self, meshData, i):
        # erzeugt eine FEniCS-Funktion des Deformationsfeldes i auf einem Mesh aus den gespeicherten Arrays
        # ermoeglicht dadurch Transport; i entspricht Index in Speichermatrix self.deformation, aufsteigend im Alter

        if isinstance(meshData, MeshData): pass
        else: raise SystemExit("initialize_defo benoetigt Objekt der MeshData-Klasse als Input!")

        V = VectorFunctionSpace(meshData.mesh, "P", 1, dim=2)
        f = Function(V)
        f.vector()[:] = self.deformation[i]
        return f


def load_mesh(name):

    # Pfad zur speziellen Gitter Datei
    path_meshFile = os.path.join(DATA_DIR, name)

    # Erstelle FEniCS Mesh mit Subdomains und Boundaries
    mesh = Mesh(path_meshFile + ".xml")
    subdomains = MeshFunction("size_t", mesh,
                              path_meshFile + "_physical_region.xml")
    boundaries = MeshFunction("size_t", mesh,
                              path_meshFile + "_facet_region.xml")

    # Berechne Indizes mit Traeger nicht am inneren Rand
    ind = __get_index_not_interior_boundary(mesh, subdomains, boundaries)

    # Rueckgabe als MeshData Objekt
    return MeshData(mesh, subdomains, boundaries, ind)

def solve_targetfunction(meshData, deformation, y_z, fValues, nu):
    """
    Berechnet den Wert des Zielfunktionals nach Verschiebung
    """

    # Verschiebe Gitter
    ALE.move(meshData.mesh, deformation)

    # Berechne Zustand in verschobenem Gitter
    y = solve_state(meshData, fValues)

    # Assembliere Zielfunktional
    V = FunctionSpace(meshData.mesh, 'P', 1)
    z = project(y_z, V)

    j = 1./2.*norm(project(y-z, V), 'L2', meshData.mesh)**2

    ds = Measure('dS', subdomain_data=meshData.boundaries)
    ones = Function(V)
    ones.vector()[:] = 1.
    j_reg_integral = ones*ds(5) + ones*ds(6)
    j_reg = nu*assemble(j_reg_integral)

    J = j + j_reg

    return J

def solve_state(meshData, fValues):
    """
    Loest Zustandsgleichung ohne Variationsungleichung
    """

    # Funktionen Raum
    V = FunctionSpace(meshData.mesh, "P", 1)

    # Abbruchkriterium
    zero_function = Expression("0.0", degree=1)
    abort = Function(V)
    abort.interpolate(zero_function)

    # Randbedingungen
    y_out = Constant(0.0)
    bcs = [ DirichletBC(V, y_out, meshData.boundaries, i) for i in range(1, 5) ]

    # Problem definieren
    dx = Measure('dx',
                 domain=meshData.mesh,
                 subdomain_data=meshData.subdomains)

    f1 = Constant(fValues[0])
    f2 = Constant(fValues[1])

    y = TrialFunction(V)
    v = TestFunction(V)

    a = inner(grad(y), grad(v))*dx('everywhere')
    b = f1*v*dx(1) + f2*v*dx(2)('everywhere')

    # Loesung berechnen
    y = Function(V, name="state_sol")
    solve(a == b, y, bcs)

    # Rueckgabe der Loesungen
    y.rename("state_sol", "label")
    return y


def solve_adjoint(meshData, y, z):
    """
    Loest Adjungierte Gleichung ohne Variationsungleichung
    """

    # Funktionen Raum
    V = FunctionSpace(meshData.mesh, "P", 1)

    # Randbedingungen
    p_out = Constant(0.0)
    bcs = [ DirichletBC(V, p_out, meshData.boundaries, i) for i in range(1, 5) ]

    # Variationsproblem definieren
    dx = Measure('dx',
                 domain=meshData.mesh,
                 subdomain_data=meshData.subdomains)

    p = TrialFunction(V)
    v = TestFunction(V)

    a = inner(grad(p), grad(v))*dx
    l = -y*v*dx + z*v*dx

    # Loesung berechnen
    p = Function(V, name="adjoint_sol")
    solve(a == l, p, bcs)

    # Rueckgabe der Loesung
    return p


def calc_lame_par(meshData, mu_min_value, mu_max_value):
    """
    Berechnet die lokal variierenden Lame-Parameter mu_elas
    """

    # Funktionen Raum
    V = FunctionSpace(meshData.mesh, "P", 1)

    # Randbedingungen
    mu_min = Constant(mu_min_value)
    mu_max = Constant(mu_max_value)
    bcs = ([ DirichletBC(V, mu_min, meshData.boundaries, i) for i in range(1, 5) ]
           + [ DirichletBC(V, mu_max, meshData.boundaries, i) for i in range(5, 7) ])

    # Variationsproblem definieren
    dx = Measure('dx',
                 domain=meshData.mesh,
                 subdomain_data=meshData.subdomains)

    mu_elas = TrialFunction(V)
    v = TestFunction(V)

    f = Expression("0.0", degree=1)

    a = inner(grad(mu_elas), grad(v))*dx
    l = f*v*dx

    # Loesung berechnen
    mu_elas = Function(V, name="lame_par")
    solve(a == l, mu_elas, bcs)

    # Rueckgabe des Lame-Parameters
    return mu_elas


def solve_linelas(meshData, p, y, z, fValues, mu_elas, nu):
    """
    Loest lineare Elastizitaetsgleichung ohne Variationsungleichung
    """

    # Funktionenraum
    V = VectorFunctionSpace(meshData.mesh, "P", 1, dim=2)

    # Randbedingungen
    u_out = Constant((0.0, 0.0))
    bcs = [ DirichletBC(V, u_out, meshData.boundaries, i) for i in range (1, 5) ]

    # Variationsproblem definieren
    dx = Measure('dx',
                 domain=meshData.mesh,
                 subdomain_data=meshData.subdomains)

    dS = Measure('dS', subdomain_data=meshData.boundaries)

    f1 = Constant(fValues[0])
    f2 = Constant(fValues[1])

    U = TrialFunction(V)
    v = TestFunction(V)
    n = FacetNormal(meshData.mesh)

    epsilon_v = sym(nabla_grad(v))
    sigma_U   = 2.0*mu_elas*sym(nabla_grad(U))

    a   = inner(sigma_U, epsilon_v)*dx + Constant(0.0)*dot(U,v)*ds
    LHS = assemble(a)

    Dj = (-inner(grad(y), dot(epsilon_v*2, grad(p)))*dx
          + nabla_div(v)*(1/2*(y-z)**2 + inner(grad(y), grad(p)))*dx
          - nabla_div(v)*f1*p*dx(1) - nabla_div(v)*f2*p*dx(2))

    Dj_reg = nu*((nabla_div(v('+'))
                  - inner(dot(nabla_grad(v('+')), n('+')), n('+')))*dS(5)
                  + (nabla_div(v('+'))
                     - inner(dot(nabla_grad(v('+')), n('+')), n('+')))*dS(6))

    F_elas = assemble(Dj + Dj_reg)

    # Alle die keine Traeger am innerend Rand haben auf 0 setzen
    F_elas[meshData.indNotIntBoundary] = 0.0

    for bc in bcs:
        bc.apply(LHS)
        bc.apply(F_elas)

    # Norm der assemblierten rechten Seite wegen Abbruchbedingung
    nrm_f_elas = norm(F_elas, 'L2', meshData.mesh)

    # Berechne Loesung
    U = Function(V, name="deformation_vec")
    solve(LHS, U.vector(), -F_elas)

    # Rueckgabe des Deformationsvektorfeldes U und der Abbruchbedingung
    return U, nrm_f_elas



def __get_index_not_interior_boundary(mesh, subdomains, boundaries):
    """
    Berechnet Indizes der Elemente mit Traeger am inneren Rand
    """

    # Facetten Indizes des inneren Randes bestimmen
    ind_interior_boundary_facets = []
    for i in range(0,len(boundaries)):
        if boundaries[i] > 4:
            ind_interior_boundary_facets.append(i)

    # Knoten Indizes des inneren Randes bestimmen
    ind_interior_boundary_vertices = []
    for c in cells(mesh):
        for f in facets(c):
            if f.index() in ind_interior_boundary_facets:
                for v in vertices(f):
                    ind_interior_boundary_vertices.append(v.index())

    ind_interior_boundary_vertices = list(set(ind_interior_boundary_vertices))

    # Element Indizes des inneren Randes bestimmen
    ind_around_interior_boundary_cells = []
    for c in cells(mesh):
        ind = False
        for v in vertices(c):
            if v.index() in ind_interior_boundary_vertices:
                ind = True
        if ind:
            ind_around_interior_boundary_cells.append(c.index())

    # Als neue Subdomain definieren
    new_sub = MeshFunction("size_t", mesh, 2)
    new_sub.set_all(0)
    for i in ind_around_interior_boundary_cells:
        if subdomains[i] == 1:
            new_sub[i] = 1
        else:
            new_sub[i] = 2

    # Indizes berchenen mit Traeger nicht am inneren Rand ueber Testproblem
    V = VectorFunctionSpace(mesh, "P", 1, dim=2)

    dx_int = Measure('dx',
                     domain=mesh,
                     subdomain_data=new_sub)

    v = TestFunction(V)

    dummy_y = Constant((1.0, 1.0))

    f_elas_int_1 = inner(dummy_y, v)*dx_int(1)
    F_elas_int_1 = assemble(f_elas_int_1)
    f_elas_int_2 = inner(dummy_y, v)*dx_int(2)
    F_elas_int_2 = assemble(f_elas_int_2)

    # Indizes setzen durch alle Punkte mit Wert 0, die keinen Einfluss haben
    ind1 = (F_elas_int_1.array() == 0.0)
    ind2 = (F_elas_int_2.array() == 0.0)

    ind = ind1 | ind2

    return ind


def bilin_a(meshData, U, V, mu_elas):
    """
    Berechnet den Wert der Bilinearform (lin. El.) fuer gegebene Vektorfelder U, V
    Beide Vektorfelder muessen auf dem selben Mesh definiert sein
    Lame parameter lambda = 0
    """

    dx = Measure("dx", domain=meshData.mesh, subdomain_data=meshData.subdomains)

    epsilon_v = sym(nabla_grad(V))
    sigma_U   = 2.0*mu_elas*sym(nabla_grad(U))

    a     = inner(sigma_U, epsilon_v)*dx('everywhere')
    value = assemble(a)

    return value


def bfgs_step(meshData, memory, mu_elas):
    """
    berechnet aus einer BFGS-memory eine Mesh-Deformation q mittels double-loop-L-BFGS-Verfahren, welche zu memory.grad[0] gehoert
    benoetigt memory.grad[0] als aktuellen Gradienten, memory.deformation[0] als aktuell neueste Deformation
    Output q ist eine Fenics-Funktion der Art Function(V), V=VectorFunctionSpace(mesh, "P", 1, dim=2)
    """

    if isinstance(meshData, MeshData): pass
    else: raise SystemExit("bfgs_step benoetigt Objekt der MeshData-Klasse als Input!")

    if isinstance(memory, bfgs_memory): pass
    else: raise SystemExit("bfgs_step benoetigt Objekt der  BFGS-Memory-Klasse als Input!")

    q     = memory.initialize_grad(meshData, 0)
    alpha = np.zeros(memory.length)

    V = VectorFunctionSpace(meshData.mesh, "P", 1, dim=2)
    diff_grad = Function(V)
    first_diff_grad = Function(V)


    if(memory.step_nr + 1 >= memory.length):
        # bei voll besetzter memory werden alle Eintraege verwendet

        for i in range(memory.length-1):
            # Vorwaertsschleife
            i         = i+1
            diff_grad.vector()[:] = memory.initialize_grad(meshData, i-1).vector() - memory.initialize_grad(meshData, i).vector()
            alpha[i]              = bilin_a(meshData, memory.initialize_defo(meshData, i-1), q, mu_elas) / bilin_a(meshData, diff_grad, memory.initialize_defo(meshData, i-1), mu_elas)
            q.vector()[:]         = q.vector() - float(alpha[i])*diff_grad.vector()

        # Reskalierung von q
        first_diff_grad.vector()[:] = memory.initialize_grad(meshData, 0).vector() - memory.initialize_grad(meshData, 1).vector()
        gamma                       = bilin_a(meshData, first_diff_grad, memory.initialize_defo(meshData, 0), mu_elas) / bilin_a(meshData, first_diff_grad, first_diff_grad, mu_elas)
        q.vector()[:]               = gamma*q.vector()

        for i in range(memory.length-1):
            # Rueckwaertsschleife
            i                     = i+1
            diff_grad.vector()[:] = memory.initialize_grad(meshData, -(i+1)).vector() - memory.initialize_grad(meshData, -i).vector()
            beta                  = bilin_a(meshData, diff_grad, q, mu_elas) / bilin_a(meshData, diff_grad, memory.initialize_defo(meshData, -(i+1)), mu_elas)
            q.vector()[:]         = q.vector() + (float(alpha[-i]) - beta)*memory.initialize_defo(meshData, -(i+1)).vector()

    elif(memory.step_nr == 0):
            # der erste BFGS-Schritt ist ein Gradientenschritt
            q.vector()[:] = -1.*q.vector()
            return q

    else:
        # bei nicht voll besetzter memory werden lediglich die besetzten Eintraege verwendet

        for i in range(memory.step_nr):
            # Vorwaertsschleife
            i         = i+1
            diff_grad.vector()[:] = memory.initialize_grad(meshData, i-1).vector() - memory.initialize_grad(meshData, i).vector()
            alpha[i]              = bilin_a(meshData, memory.initialize_defo(meshData, i-1), q, mu_elas) / bilin_a(meshData, diff_grad, memory.initialize_defo(meshData, i-1), mu_elas)
            q.vector()[:]         = q.vector() - float(alpha[i])*diff_grad.vector()

        # Reskalierung von q
        first_diff_grad.vector()[:] = memory.initialize_grad(meshData, 0).vector() - memory.initialize_grad(meshData, 1).vector()
        gamma                       = bilin_a(meshData, first_diff_grad, memory.initialize_defo(meshData, 0), mu_elas) / bilin_a(meshData, first_diff_grad, first_diff_grad, mu_elas)
        q.vector()[:]               = gamma*q.vector()

        for i in range(memory.step_nr):
            # Rueckwaertsschleife
            shift = (memory.length-1) - memory.step_nr
            i                     = i+1
            diff_grad.vector()[:] = memory.initialize_grad(meshData, -(i+1)-shift).vector() - memory.initialize_grad(meshData, -i-shift).vector()
            beta                  = bilin_a(meshData, diff_grad, q, mu_elas) / bilin_a(meshData, diff_grad, memory.initialize_defo(meshData, -(i+1)-shift), mu_elas)
            q.vector()[:]         = q.vector() + (float(alpha[-i-shift]) - beta)*memory.initialize_defo(meshData, -(i+1)-shift).vector()

    q.vector()[:] = -1. * q.vector()
    return q

