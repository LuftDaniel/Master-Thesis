from fenics import *
import sovi_bib as sovi
import time
import argparse
import os
import numpy as np

parameters['allow_extrapolation'] = True

# ----------------------- #
#    BENUTZER AUSWAHL     #
# ----------------------- #

# Waehle die konstanten Funktionswerte in den zwei Gebieten
f1 = -10.0
f2 = 100.0

# Definition des minimalen und maximalen Lameparameters
mu_min = 1.0
mu_max = 30.0

# ----------------------- #
#  BENUTZER AUSWAHL ENDE  #
# ----------------------- #

print("####################################################")
print("# Shape Optimization with Variational Inequalities #")
print("####################################################\n")

string_mesh_desription = """\nFunktionswerte f1 und f2 und
                            \nmu_min und mu_max f1 und f2 in Datei aendern!
                            \nKombinationsmoeglichkeiten:
                            \n [1] kleiner Kreis -> grosser Kreis
                            \n [2] grosser Kreis -> kleiner Kreis
                            \n [3] Form -> kleiner Kreis
                            \n [4] kleiner Kreis -> verschobener Kreis
                            \n [5] verschobener Kreis -> kleiner Kreis
                            """

# ----------------------- #
#      INPUT ABFRAGE      #
# ----------------------- #

parser = argparse.ArgumentParser(description=string_mesh_desription,
                                 formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('-m','--mesh', type=int,
                    help='Waehle Gitter-Kombinationen (1, 2, 3, 4, 5)',
                    required=True)
parser.add_argument('-p','--psi', type=float,
                    help='Setze konstante Beschraenkung psi (z.B. 0,7 oder 0.7)',
                    required=True)
args = parser.parse_args()

# ----------------------- #
#   INPUT ABFRAGE ENDE    #
# ----------------------- #

# Liste aller Gitter-Kombinationen
mesh_combination_list = [("mesh_smallercircle", "mesh_circle"),
                         ("mesh_circle",        "mesh_smallercircle"),
                         ("mesh_form",          "mesh_smallercircle"),
                         ("mesh_smallercircle", "mesh_leftbottom"),
                         ("mesh_leftbottom",    "mesh_smallercircle")]

# Anzahl aller Kombinationen
n = len(mesh_combination_list)

# Ueberpruefe ob gewaehlte Kombination vorhanden, sonst Ende
if args.mesh - 1 >= n or args.mesh <= 0:
    print("Warnung: Waehle Kombination zwischen 1 und {0:1d}!".format(n))
    quit()

# Falls Komma als Dezimaltrenner umwandeln in float
if type(args.psi) is tuple:
    if len(args.psi) == 2:
        psi_parsed = float(str(args.psi[0]) + "." + str(args.psi[1]))
    else:
        print("Falsche Eingabe")
        quit()
else:
    psi_parsed = float(args.psi)

# Setze Start- und Zielgitter je nach Benutzerwahl aus der Liste
(startMesh_file, targetMesh_file) = mesh_combination_list[args.mesh - 1]

# ----------------------- #
#       PARAMETER         #
# ----------------------- #

# Setze Funktionswerte als Array
f_values = [f1, f2]

# Setze psi als Funktion
psi_values = Constant(psi_parsed) # Expression(psi_string, degree=2)

# Parameter fuer den PDAS Algorithmus
c       = 5.0
tol_ssn = 1.e-3

# Parameter fuer L-BFGS Algorithmus
memory_length = 3

# Parameter fuer den Formoptimierungsalgorithmus
nu        = 0.01
tol_shopt = 2.e-4

# Keine Ausgabe von FEniCS Funktionen in der Konsole
set_log_active(False)

# Generiere Output Ordner falls nicht vorhanden
outputfolder = sovi.create_outputfolder()

# Speichere Ausgabedateien fuer ParaView im Output Ordner
file_mesh        = File(os.path.join(outputfolder,
                                     'Mesh',            'mesh.pvd'))
file_sol_y       = File(os.path.join(outputfolder,
                                     'StateSolution',   'sol_y.pvd'))
file_adj_p       = File(os.path.join(outputfolder,
                                     'AdjointSolution', 'adj_p.pvd'))
file_def_u       = File(os.path.join(outputfolder,
                                     'Deformation',     'def_u.pvd'))
file_data_sol    = File(os.path.join(outputfolder,
                                     'TargetData',      'solution.pvd'))
file_data_lambda = File(os.path.join(outputfolder,
                                     'TargetData',      'lambda.pvd'))
file_data_mesh   = File(os.path.join(outputfolder,
                                     'TargetData',      'mesh.pvd'))
file_data_bound  = File(os.path.join(outputfolder,
                                     'TargetData',      'bound.pvd'))
file_lame        = File(os.path.join(outputfolder,
                                     'LameParameter',   'lamepar.pvd'))
file_bound_end   = File(os.path.join(outputfolder,
                                     'BoundaryData',   'bound_end.pvd'))
file_bound_start = File(os.path.join(outputfolder,
                                     'BoundaryData',   'bound_start.pvd'))

# Zusammenfassung aller Parameter als Ausgabe in der Konsole
print("\n-----------------PARAMETER-----------------------\n")
print("Zielgitter Dateiname:     " + targetMesh_file)
print("Ausgangsgitter Dateiname: " + startMesh_file)
print("\npsi    = {0:3.1f}".format(psi_parsed))
print("c      = {0:3.1f}".format(c))
print("mu_min = {0:3.0f}".format(mu_min))
print("mu_max = {0:3.0f}".format(mu_max))
print("\ntol    = {0:1.0E}".format(tol_shopt))
print("\n-------------------------------------------------")

# Zusammenfassung aller Parameter als Ausgabe in csv-Datei
file_csv = open(os.path.join(outputfolder, 'output.csv'), 'a')
file_csv.write("Gitter; Start; "  + startMesh_file                    + "\n")
file_csv.write("; Ziel; "       + targetMesh_file                   + "\n")
file_csv.write("Parameter; psi; " + str(psi_parsed).replace(".", ",") + "\n")
file_csv.write("; c; "            + str(c).replace(".", ",")          + "\n")
file_csv.write("; mu_min; "       + str(mu_min).replace(".", ",")     + "\n")
file_csv.write("; mu_max; "       + str(mu_max).replace(".", ",")     + "\n")
file_csv.write("; tol; "          + str(tol_shopt).replace(".", ",")  + "\n")
file_csv.write("Iteration; ||f_elas||_L2;J=j+j_reg;||U||_L2"            + "\n")

# ----------------------- #
#   ZIELDATEN BERECHNEN   #
# ----------------------- #

print('\nSchritt 1/2: Zieldaten berechnen\n')

# Gitter der Zieldaten laden
targetMeshData = sovi.load_mesh(targetMesh_file)

# Loesen der Zustandsgleichung auf dem Zielgitter im Variationsfall
# y_z, _, lmbda_c_squared = sovi.solve_state_vi(targetMeshData, f_values,
#                                                  psi_values, c, tol_ssn)

# Loesen der Zustandsgleichung auf dem Zielgitter ohne Variationsungleichung
y_z = sovi.solve_state(targetMeshData, f_values)
lmbda_c_squared = y_z

# Parameter Normalverteilung fuer Stoerung
mu = Constant(0.0)

#y_z_min = y_z.vector().array().min()
#y_z_max = y_z.vector().array().max()
#sigma = Constant(max(abs(y_z_min), y_z_max)*5/100)
#
#y_z.vector()[:] = y_z.vector()+np.random.normal(mu, sigma,
#                                                y_z.vector().size())

# Loesung y, lambda_c_squared, und das Gitter in pvd-Datei speichern
file_data_sol << y_z
file_data_lambda << lmbda_c_squared
file_data_mesh << targetMeshData.mesh
file_data_bound << targetMeshData.boundaries

# ----------------------- #
#   FORM OPTIMIERUNG      #
# ----------------------- #

print('\nSchritt 2/2: Formoptimierung')

# Gitter laden und speichern
MeshData = sovi.load_mesh(startMesh_file)
file_bound_start << MeshData.boundaries

# file_mesh << MeshData.mesh

# BFGS-memory initialisieren
bfgs_memory   = sovi.bfgs_memory(np.zeros([memory_length, 2 * MeshData.mesh.num_vertices()]),
                                 np.zeros([memory_length, 2 * MeshData.mesh.num_vertices()]), memory_length, 0)

# Lame-Parameter berechnen und speichern
mu_elas = sovi.calc_lame_par(MeshData, mu_min, mu_max)

# Dummy um Schleife zu durchlaufen
nrm_f_elas = 1.0

# Zaehler der Iterationen
counter = 0

# Startzeit zum Zeitmessen
start_time = time.time()

# Start der Optimierungsschritte
print("\nIteration " + "  ||f_elas||_L2 " + "  J = j + j_reg " + "  ||U||_L2\n")
# Start der Optimierungsschritte
while nrm_f_elas > tol_shopt:

    # Zaehler erhoehen und ausgeben in Konsole
    counter += 1

    # ----------------------- #
    #      INTERPOLATION      #
    # ----------------------- #

    V = FunctionSpace(MeshData.mesh, "P", 1)
    z = project(y_z, V)
    mu_elas_projected = project(mu_elas, V)
    mu_elas_projected.rename("lame_par", "label")

    # ----------------------- #
    #     ZUSTANDSGLEICHUNG   #
    # ----------------------- #

    # Im Variationsfall
    # y, lmbda_c, lmbda_c_squared = sovi.solve_state_vi(MeshData, f_values,
    #                                                 psi_values, c, tol_ssn)

    # Ohne Variationsungleichung
    y = sovi.solve_state(MeshData, f_values)

    # ----------------------- #
    #      ADJUNGIERTE        #
    # ----------------------- #

    # Im Variationsfall
    # p = sovi.solve_adjoint_vi(MeshData, y, z, lmbda_c, c)

    # Ohne Variationsungleichung
    p = sovi.solve_adjoint(MeshData, y, z)

    # ----------------------- #
    #     LIN. ELAS. PDE      #
    # ----------------------- #

    # Im Variationsfall
    # U, nrm_f_elas = sovi.solve_linelas_vi(MeshData, p, y, lmbda_c_squared, z,
    #                                   f_values, mu_elas_projected, nu)

    # Ohne Variationsungleichung
    U , nrm_f_elas = sovi.solve_linelas(MeshData, p, y, z, f_values, mu_elas_projected, nu)

    # L-BFGS Schritt mit Memory-Updates
    bfgs_memory.update_grad(U.vector().get_local())
    S = sovi.bfgs_step(MeshData, bfgs_memory, mu_elas_projected)

    bfgs_memory.update_defo(S.vector().get_local())
    bfgs_memory.step_nr = bfgs_memory.step_nr + 1

    # Norm des Deformationsvektorfeldes berechnen
    V           = FunctionSpace(MeshData.mesh, 'P', 1)
    U_magnitude = sqrt(dot(U, U))
    U_magnitude = project(U_magnitude, V)
    nrm_U_mag   = norm(U_magnitude, 'L2', MeshData.mesh)

    # Zielfunktional berechnen
    j = 1./2.*norm(project(y-z,V), 'L2', MeshData.mesh)

    ds = Measure('dS', subdomain_data=MeshData.boundaries)
    ones = Function(V)
    ones.vector()[:] = 1.
    j_reg_integral = ones*ds(5) + ones*ds(6)
    j_reg = nu*assemble(j_reg_integral)

    J = j + j_reg

    # Iterationsfortschritt in Konsole ausgeben
    print("   {0:3d}   ".format(counter)    + " | " +
          "   {0:.2E}  ".format(nrm_f_elas) + " | " +
          "   {0:.2E}  ".format(J)   + " | " +
          " {0:.2E}".format(nrm_U_mag))

    # Iterationsschritt in csv-Datei speichern
    file_csv.write(str(counter)                       +";" +
                   str(nrm_f_elas).replace(".", ",")  + ";" +
                   str(J).replace(".", ",") + ";" +
                   str(nrm_U_mag).replace(".", ",")   + "\n")

    # ----------------------- #
    #      DEFORMATION        #
    # ----------------------- #

    ALE.move(MeshData.mesh, S)

# ----------------------- #
#          FERTIG         #
# ----------------------- #

    # ----------------------- #
    #         OUTPUT          #
    # ----------------------- #

    # Speichern in pvd-Datei fuer ParaView
file_mesh      << MeshData.mesh
file_sol_y     << y
file_adj_p     << p
file_def_u     << U
file_lame      << mu_elas_projected
file_bound_end << MeshData.boundaries



# Zeitmessung
elapsed_time = time.time() - start_time

# Dauer der Berechnung in csv-Datei speichern
file_csv.write("Dauer der Berechnung in Sekunden;" +
               str(elapsed_time).replace(".", ",") + ";" + "\n")

# Datei schliessen
file_csv.close()

# Abschluss Ausgabe in Konsole
print("\n\n-------------------Fertig------------------------\n")
print("Iterationen:            {0:0d}".format(counter))
print("Dauer der Berechnungen: {0:.0f} min".format(int(elapsed_time)/60))
print("Ergebnisse in:          '{0:s}/'".format(outputfolder))
print("\n-------------------------------------------------\n")
