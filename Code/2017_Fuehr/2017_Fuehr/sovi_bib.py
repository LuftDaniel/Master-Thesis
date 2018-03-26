from fenics import *
import os
import errno
from datetime import datetime
from copy import copy

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

    # Objekt mit gespeicherten Gradienten und Deformationen der letzten l Schritte
    def __init__(self, gradient, deformation, length):

        # Liste von Gradientenvektorfeldern
        if (len(gradient) == length): self.gradient = gradient
        else: raise SystemExit("Fehler: Anzahl der Gradienten passt nicht zur Memorylaenge!")

        # Liste von Deformationsvektorfeldern
        if (len(deformation) == length): self.deformation = deformation
        else: raise SystemExit("Fehler: Anzahl der Deformationen passt nicht zur Memorylaenge!")

        # Anzahl der gespeicherten letzten Schritte
        self.length = length

    # macht ein Update der Memory; neueste Elemente bekommen Index 0
    def update(self, upd_grad, upd_defo):

        for i in range(self.length-1):

            self.gradient[-(i+1)]    = self.gradient[-(i+2)]
            self.deformation[-(i+1)] = self.deformation[-(i+2)]

        self.gradient[0]    = upd_grad
        self.deformation[0] = upd_defo


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


def solve_state_vi(meshData, fValues, psi_values, c, tol):
    """
    Loest Zustandsgleichung mit Variationsungleichung
    """

    # Funktionen Raum
    V = FunctionSpace(meshData.mesh, "P", 1)

    # Psi als Funktion beschreiben
    psi = Function(V)
    psi.interpolate(psi_values)

    # Lambda initialisieren
    zero_function = Expression("0.0", degree=1)
    lmbda_c       = Function(V)
    lmbda_c.interpolate(zero_function)

    # Abbruchkriterium
    abort = Function(V)
    abort.interpolate(zero_function)

    # Randbedingungen
    y_out = Constant(0.0)
    bcs = [ DirichletBC(V, y_out, meshData.boundaries, i) for i in range(1, 5) ]

    # Variationsproblem definieren
    dx = Measure('dx',
                 domain=meshData.mesh,
                 subdomain_data=meshData.subdomains)

    f1 = Constant(fValues[0])
    f2 = Constant(fValues[1])

    y = TrialFunction(V)
    v = TestFunction(V)

    a = inner(grad(y), grad(v))*dx('everywhere')
    b = f1*v*dx(1) + f2*v*dx(2) - lmbda_c*v*dx('everywhere')

    # Loesung berechnen
    y = Function(V, name="state_sol")
    solve(a == b, y, bcs)

    # Abbruchkriterium initialisieren
    nrm = 1.0
    while nrm > tol:

        # Maximumsfunktion definieren
        max_func = project(lmbda_c+c*(y-psi))
        ind = max_func.vector()<=DOLFIN_EPS
        max_func.vector()[ind] = 0.0

        # Charakteristische Funktion definieren
        chi_delta_y = Function(V)
        chi_delta_y.interpolate(Expression("1.0", degree=1))
        chi_delta_y.vector()[ind] = 0.0

        for i in range(0,2):
            delta_y = TrialFunction(V)
            lhs = (inner(grad(delta_y), grad(v))*dx('everywhere')
                   + 2*c*(lmbda_c+c*(y-psi))*delta_y*v*chi_delta_y*dx('everywhere'))
            rhs = (f1*v*dx(1) + f2*v*dx(2)
                   - inner(grad(y), grad(v))*dx('everywhere')
                   - (lmbda_c+c*(y-psi))**2*v*chi_delta_y*dx)

            # Loesung berechnen
            delta_y = Function(V)
            solve(lhs == rhs, delta_y, bcs)

            y = project(y + delta_y, V)

        # Lambda updaten und Abbruchkriterium berechnen
        lmbda_c = project(lmbda_c+c*(y-psi))
        ind = lmbda_c.vector()<=DOLFIN_EPS
        abort = project(y-psi)

        lmbda_c.vector()[ind] = 0.0
        abort.vector()[ind] = 0.0

        nrm = norm(abort, 'L2', meshData.mesh)

    # Lambda_c quadrieren
    lmbda_c_squared = interpolate(Expression("f*f",
                                             f = interpolate(lmbda_c, V),
                                             degree=1),V)
    # Rueckgabe der Loesungen
    y.rename("state_sol", "label")
    return y, lmbda_c, lmbda_c_squared


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


def solve_adjoint_vi(meshData, y, z, lmbda_c_root, c):
    """
    Loest Adjungierte Gleichung im Fall einer Variationsungleichung
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

    a = inner(grad(p), grad(v))*dx + 2*c*lmbda_c_root*p*v*dx
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

    # Funktionen Raum
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


def solve_linelas_vi(meshData, p, y, lmbda_c, z, fValues, mu_elas, nu):
    """
    Loest lineare Elastizitaetsgleichung mit Variationsungleichung
    """

    # Funktionen Raum
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
          + nabla_div(v)*(1/2*(y-z)**2 + inner(grad(y), grad(p))+lmbda_c*p)*dx
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
    berechnet aus einer BFGS-memory eine Mesh-Deformation q mittels double-loop-L-BFGS-Verfahren
    Transporte der Gradienten und Deformationen sind vereinfacht, lediglich Fußpunkte sind retracted
    benötigt memory.grad[0] als aktuellen Gradienten
    """
    if isinstance(memory, bfgs_memory): pass
    else: raise SystemExit("bfgs_step benoetigt eine BFGS-Memory als input!")

    q = memory.gradient[0]
    alpha = np.zeros(memory.length)

    for i in range(memory.length-1):
        # Vorwärtsschleife
        i = i+1
        diff_grad = memory.gradient[i-1] - memory.gradient[i]
        alpha[i] = bilin_a(meshData, memory.deformation[i], q, mu_elas) / bilin_a(meshData, diff_grad, memory.deformation[i], mu_elas)
        q = q - alpha[i]*diff_grad

    # Reskalierung von q
    first_diff_grad = memory.gradient[0] - memory.gradient[1]
    gamma = bilin_a(meshData, first_diff_grad, memory.deformation[1], mu_elas) / bilin_a(meshData, first_diff_grad, first_diff_grad, mu_elas)
    q = gamma*q

    for i in range(memory.length-1):
        # Rückwärtsschleife
        i = i+1
        diff_grad = memory.gradient[-i-1] - memory.gradient[-i]
        beta = bilin_a(meshData, diff_grad, q, mu_elas) / bilin_a(meshData, diff_grad, memory.deformation[-i], mu_elas)
        q = q + (alpha[-i] - beta)*diff_grad

    return q

