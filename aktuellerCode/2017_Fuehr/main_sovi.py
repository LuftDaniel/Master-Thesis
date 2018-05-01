from fenics import *
import sovi_bib as sovi
import time
import argparse
import os
import numpy as np
import matplotlib.pyplot as plt

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

# Parameter fuer L-BFGS Algorithmus
memory_length = 3

# Parameter fuer den Formoptimierungsalgorithmus
nu        = 0.01
tol_shopt = 2.e-4


# Parameter fuer Backtracking-Linesearch
shrinkage = 0.5
c = 0.99
start_scale = 5.0

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
                            \n [6] kleiner Kreis hoch aufgeloest -> grosser Kreis hoch aufgeloest
                            \n [7] grosser Kreis hoch aufgeloest -> kleiner Kreis hoch aufgeloest                            
                            """

# ----------------------- #
#      INPUT ABFRAGE      #
# ----------------------- #

parser = argparse.ArgumentParser(description=string_mesh_desription,
                                 formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('-m','--mesh', type=int,
                    help='Waehle Gitter-Kombinationen (1, 2, 3, 4, 5, 6, 7)',
                    required=True)
args = parser.parse_args()

# ----------------------- #
#   INPUT ABFRAGE ENDE    #
# ----------------------- #

# Liste aller Gitter-Kombinationen
mesh_combination_list = [("mesh_smallercircle",      "mesh_circle"),
                         ("mesh_circle",             "mesh_smallercircle"),
                         ("mesh_form",               "mesh_smallercircle"),
                         ("mesh_smallercircle",      "mesh_leftbottom"),
                         ("mesh_leftbottom",         "mesh_smallercircle"),
                         ("mesh_fine_smallercircle", "mesh_fine_circle"),
                         ("mesh_fine_circle",        "mesh_fine_smallercircle")]

# Anzahl aller Kombinationen
n = len(mesh_combination_list)

# Ueberpruefe ob gewaehlte Kombination vorhanden, sonst Ende
if args.mesh - 1 >= n or args.mesh <= 0:
    print("Warnung: Waehle Kombination zwischen 1 und {0:1d}!".format(n))
    quit()

# Setze Start- und Zielgitter je nach Benutzerwahl aus der Liste
(startMesh_file, targetMesh_file) = mesh_combination_list[args.mesh - 1]

# ----------------------- #
#       PARAMETER         #
# ----------------------- #

# Setze Funktionswerte als Array
f_values = [f1, f2]

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
                                     'BoundaryData',      'bound.pvd'))
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
print("mu_min = {0:3.0f}".format(mu_min))
print("mu_max = {0:3.0f}".format(mu_max))
print("\ntol    = {0:1.0E}".format(tol_shopt))
print("\n-------------------------------------------------")

# Zusammenfassung aller Parameter als Ausgabe in csv-Datei
file_csv = open(os.path.join(outputfolder, 'output.csv'), 'a')
file_csv.write("Gitter; Start; "  + startMesh_file                    + "\n")
file_csv.write("; Ziel; "       + targetMesh_file                   + "\n")
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

# Loesen der Zustandsgleichung auf dem Zielgitter
y_z = sovi.solve_state(targetMeshData, f_values)

# Parameter Normalverteilung fuer Stoerung
mu = Constant(0.0)

y_z_min = y_z.vector().get_local().min()
y_z_max = y_z.vector().get_local().max()
sigma = Constant(max(abs(y_z_min), y_z_max)*5/100)

#y_z.vector()[:] = y_z.vector()+np.random.normal(mu, sigma,
#                                                y_z.vector().size())

# Loesung y, lambda_c_squared, und das Gitter in pvd-Datei speichern
file_data_sol << y_z
file_data_mesh << targetMeshData.mesh
file_data_bound << targetMeshData.boundaries

# ----------------------- #
#   FORM OPTIMIERUNG      #
# ----------------------- #

print('\nSchritt 2/2: Formoptimierung')

# Gitter laden und speichern
MeshData = sovi.load_mesh(startMesh_file)
file_bound_start << MeshData.boundaries, 0


print("initial meshdistance: {0:8e}".format(sovi.mesh_distance(MeshData, targetMeshData)))
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

# Outputgraph
file_output = open(os.path.join(outputfolder, 'outputdata.txt'),'a')

# Formableitungen fuer Curv. Cond; 0 aktuelles Gitter, 1 vorheriges
curv_cond = np.zeros(2)

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

    # Ohne Variationsungleichung
    y = sovi.solve_state(MeshData, f_values)

    # ----------------------- #
    #         ADJOINT         #
    # ----------------------- #

    # Ohne Variationsungleichung
    p = sovi.solve_adjoint(MeshData, y, z)

    # ----------------------- #
    #   GRADIENT CALCULATION  #
    # ----------------------- #

    # Loese lineare Elastizitaetsgleichung
    U , nrm_f_elas = sovi.solve_linelas(MeshData, p, y, z, f_values, mu_elas_projected, nu)

    # ----------------------- #
    #       L-BFGS STEP       #
    # ----------------------- #

    # Curvature Condition
    V_2 = VectorFunctionSpace(MeshData.mesh, "P", 1, dim=2)
    last_defo = Function(V_2)
    last_defo.vector()[:] = bfgs_memory.deformation[0]
    curv_cond[0] = sovi.shape_deriv(MeshData, p, y, z, f_values, nu, last_defo)

    if(counter > 1):
        curv_cond_val = curv_cond[1] - curv_cond[0]
        print("Curvature Condition value: {0:4e}".format(curv_cond_val))
        if(curv_cond_val <= 0.):
            print("Condition not fullfilled! Exiting process, last step was invalid!")
            break

    # L-BFGS Schritt
    bfgs_memory.update_grad(U.vector().get_local())
    S = sovi.bfgs_step(MeshData, bfgs_memory, mu_elas_projected)

    # L-BFGS-Memory Update
    bfgs_memory.update_defo(S.vector().get_local())
    bfgs_memory.step_nr = bfgs_memory.step_nr + 1

    # speichere Formableitung fuer Berechnung der curv. cond. nach der naechsten Deformation,
    # d.h. in der naechsten Schleife (deshalb curv. cond. der vorrigen Bedingung VOR BFGS-Schritt

    curv_cond[1] = sovi.shape_deriv(MeshData, p, y, z, f_values, nu, S)

    # TEST ZU SKALARPRODUKT UND ABLEITUNG; ECHTER UND FALSCHER GRADIENT
    U_real, nrm_real = sovi.solve_linelas(MeshData, p, y, z, f_values, mu_elas_projected, nu, zeroed = False)

    #loler = sovi.shape_deriv(MeshData, p, y, z, f_values, nu, S)
    #loler2 = sovi.shape_deriv(MeshData, p, y, z, f_values, nu, U_real)

    #print("shape derivative bfgs {0:.10e}".format(loler))
    #print("bilin_a {0:.10e}".format(sovi.bilin_a(MeshData,S,-U_real,mu_elas_projected)))
    #print("shape derivative gradient {0:.10e}".format(loler2))
    #print("bilin_a {0:.10e}".format(sovi.bilin_a(MeshData,U_real,-U_real,mu_elas_projected)))

            # ----------------------- #
            # BACKTRACKING-LINESEARCH # die Frage ist, ob ich die runterskalierten Felder weiterverwende
            # ----------------------- # oder die oben berechneten vollen

    #print("Achtung: Wert des deformierten Gitters {0:.2E}  ".format(Jvalue))
    #print("Achtung: Wert der Ableitung {0:.2E}  ".format(sovi.shape_deriv(MeshData, p, y, z, f_values, nu, S)))
    #print("Achtung: Wert der bilin {0:.2E}  ".format(sovi.bilin_a(MeshData, U, S, mu_elas_projected)))

    scale_parameter = start_scale

    zero_function = Function(VectorFunctionSpace(MeshData.mesh, "P", 1, dim=2))
    current_value = sovi.solve_targetfunction(MeshData, zero_function, y_z, f_values, nu)
    S.vector()[:] = start_scale * S.vector()

    counterer = 0
    #print(current_value)
    #print(current_deriv)
    U_real, trash = sovi.solve_linelas(MeshData, p, y, z , f_values, mu_elas_projected, nu, zeroed = False)
    #print(sovi.bilin_a(MeshData, U_real, U, mu_elas_projected))

    # Skaliert das Deformationsfeld bei jedem Schritt automatisch dauerhaft: NOCH OHNE AMIJO
    #while(sovi.solve_targetfunction(MeshData, S, y_z, f_values, nu) > current_value + c*scale_parameter*current_deriv):
    while(sovi.solve_targetfunction(MeshData, S, y_z, f_values, nu) > current_value):
    #    print("shape deriv: {0:3e}".format(sovi.shape_deriv(MeshData, p, y, z, f_values, nu, S)))
        #print("current Value: {0:3e}".format(current_value))
        #print("next Value: {0:3e}\n".format(sovi.solve_targetfunction(MeshData, S, y_z, f_values, nu)))
        #print("bilin deriv: {0:3e}".format(sovi.bilin_a(MeshData,S,-U_real,mu_elas_projected)))

        scale_parameter = shrinkage * scale_parameter
        S.vector()[:] = shrinkage * S.vector()

        counterer = counterer + 1
        if(counterer > 20):
            print("had to break")
            break


    # Norm des Deformationsvektorfeldes berechnen
    V           = FunctionSpace(MeshData.mesh, 'P', 1)
    U_magnitude = sqrt(dot(U, U))
    U_magnitude = project(U_magnitude, V)
    nrm_U_mag   = norm(U_magnitude, 'L2', MeshData.mesh)

    # Zielfunktional berechnen: KANN ERSETZT WERDEN AUS BIB
    j = 1./2.*norm(project(y-z,V), 'L2', MeshData.mesh)**2

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

    # Norm von f_elas plotten
    #file_output.write(str(counter)       +";"+
    #                  str(nrm_f_elas)  + "\n")

    # Abstand zur Zielform plotten
    file_output.write(str(counter)       +";"+
                      str(np.log(sovi.mesh_distance(MeshData, targetMeshData)))  + "\n")

    # ----------------------- #
    #      DEFORMATION        #
    # ----------------------- #

    # Gradientenverfahren
    #ALE.move(MeshData.mesh, U)

    # L-BFGS-Verfahren
    ALE.move(MeshData.mesh, S)

    # ----------------------- #
    #     GITTER SPEICHERN    #
    # ----------------------- #

    file_bound_increment = File(os.path.join(outputfolder,
                                         'BoundaryData', 'bound_{}.pvd'.format(counter)))
    #MeshData.boundaries.rename("bound_{}".format(counter), "bound_{}".format(counter))
    file_bound_increment << MeshData.boundaries, counter

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

file_output.close()

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

# Konvergenzgraph printen

with open(os.path.join(outputfolder, 'outputdata.txt')) as output:
    items  = (map(float, line.split(";")) for line in output)
    xs, ys = zip(*items)

plt.plot(xs, ys)
plt.show()
