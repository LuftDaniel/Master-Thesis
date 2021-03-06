from fenics import *
import shape_bib as bib
import time
import argparse
import os
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

parameters['allow_extrapolation'] = True

# ----------------------- #
#    BENUTZER AUSWAHL     #
# ----------------------- #


# Waehle die konstanten Funktionswerte in den zwei Gebieten
f1 = -10.0
f2 = 100.0

# Waehle, ob das Zielgitter pertubiert werden soll
Pertubation = False
sigma = Constant(0.01)

# Definition des minimalen und maximalen Lameparameters
mu_min = 1.0
mu_max = 30.0

# Parameter fuer L-BFGS Algorithmus
L_BFGS = False
memory_length = 60

# Parameter fuer Perimeter-Regularisierung und Abbruchkriterium fuer Optimierung
nu        = 0.00001
tol_shopt = 9.e-4

# Parameter fuer Backtracking-Linesearch
Linesearch = True
Resetcounter = 0
shrinkage = 0.5
c = 0.9999
start_scale = 5.0

# ----------------------- #
#  BENUTZER AUSWAHL ENDE  #
# ----------------------- #

print("\n#############################################")
print("# Shape Optimization with L-BFGS-Algorithms #")
print("#############################################")

string_mesh_desription = """\nParameter und Verfahren muessen in Datei geaendert werden!
                            \nFalls L-BFGS-Verfahren ausgeschaltet, so wird ein Gradientenverfahren verwendet!
                            \nKombinationsmoeglichkeiten:
                            \n [1] kleiner Kreis -> grosser Kreis
                            \n [2] grosser Kreis -> kleiner Kreis
                            \n [3] verschobener Kreis -> kleiner Kreis
                            \n [4] kleiner Kreis -> verschobener Kreis
                            \n [5] verschobener Kreis -> kleiner Kreis
                            \n [6] kleiner Kreis hoch aufgeloest -> grosser Kreis hoch aufgeloest
                            \n [7] grosser Kreis hoch aufgeloest -> kleiner Kreis hoch aufgeloest
                            \n [8] kleiner Kreis -> kaputter Donut       
                            \n [9] kleiner Kreis hoch aufgeloest -> kaputter Donut hoch aufgeloest                     
                            """

# ----------------------- #
#      INPUT ABFRAGE      #
# ----------------------- #

parser = argparse.ArgumentParser(description=string_mesh_desription,
                                 formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('-m','--mesh', type=int,
                    help='Waehle Gitter-Kombinationen (1, 2, 3, 4, 5, 6, 7, 8, 9)',
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
                         ("mesh_fine_circle",        "mesh_fine_smallercircle"),
                         ("mesh_smallercircle",      "mesh_broken_donut"),
                         ("mesh_fine_smallercircle", "mesh_fine_broken_donut")]

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
outputfolder = bib.create_outputfolder()

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
                                     'BoundaryData',      'bound_target.pvd'))
file_lame        = File(os.path.join(outputfolder,
                                     'LameParameter',   'lamepar.pvd'))
file_bound_end   = File(os.path.join(outputfolder,
                                     'BoundaryData',   'bound_end.pvd'))
file_bound_start = File(os.path.join(outputfolder,
                                     'BoundaryData',   'bound_start.pvd'))

# Zusammenfassung aller Parameter als Ausgabe in der Konsole
print("\n-----------------PARAMETER------------------------------\n")
print("Zielgitter Dateiname:                   " + targetMesh_file)
print("Ausgangsgitter Dateiname:               " + startMesh_file)
print("Perimeter Regularisierung:              " + str(nu))
print("minimaler/ maximaler Lame-Parameter:    " + str(mu_min) + "   " + str(mu_max))
print("Backtracking-Linesearch:                " + str(Linesearch))
print("L-BFGS (bei False Gradientenverfahren): " + str(L_BFGS))
print("Abbruch bei Toleranz:                   " + str(tol_shopt))
print("\n--------------------------------------------------------")

# Zusammenfassung aller Parameter als Ausgabe in csv-Datei, welche in OPENOFFICE CALC geoffnet wird
file_csv = open(os.path.join(outputfolder, 'output.csv'), 'a')
file_csv.write("Gitter; Start; "  + startMesh_file                    + "\n")
file_csv.write("; Ziel; "       + targetMesh_file                   + "\n")
file_csv.write("; mu_min; "       + str(mu_min).replace(".", ",")     + "\n")
file_csv.write("; mu_max; "       + str(mu_max).replace(".", ",")     + "\n")
file_csv.write("; tol; "          + str(tol_shopt).replace(".", ",")  + "\n")
file_csv.write("; Perimeter; "          + str(nu).replace(".", ",")  + "\n")
if(L_BFGS):            file_csv.write("Iteration; ||f_elas||_L2;J=j+j_reg;||U||_L2; Curv. Cond.; Meshdistance" + "\n")
elif(L_BFGS == False): file_csv.write("Iteration; ||f_elas||_L2;J=j+j_reg;||U||_L2; Meshdistance" + "\n")


# -----------------------#
# TARGETDATA CALCULATION #
# -----------------------#

print('\nSchritt 1/2: Zieldaten berechnen\n')

# Gitterdaten laden
if(Pertubation):
    # laed Optimierungsproblem, bei welchem das Startgitter gestoert wird,
    # und das Zielgitter das ungestoerte Gitter ist

    MeshData       = bib.load_mesh(startMesh_file)
    targetMeshData = bib.load_mesh(startMesh_file)

    # Parameter Normalverteilung fuer Stoerung, sigma vorsichtig waehlen!
    mu = Constant(0.0)

    # finde Indizes der Knoten am Inneren Rand
    boundary_index_target = bib.__get_index_not_interior_boundary(MeshData.mesh, MeshData.subdomains,
                                                                  MeshData.boundaries, interior=False)

    # Stoere das Zielgitter!
    for i in boundary_index_target:
        MeshData.mesh.coordinates()[i,0] += np.random.normal(mu, sigma)
        MeshData.mesh.coordinates()[i,1] += np.random.normal(mu, sigma)

elif(Pertubation == False):
    # Problem wie in der Beschreibung wird geladen
    MeshData       = bib.load_mesh(startMesh_file)
    targetMeshData = bib.load_mesh(targetMesh_file)


# Loesen der Zustandsgleichung auf dem Zielgitter
y_z = bib.solve_state(targetMeshData, f_values)

# Loesung y, lambda_c_squared, und die Gitter in pvd-Datei speichern
file_data_sol << y_z
file_data_mesh << targetMeshData.mesh
file_data_bound << targetMeshData.boundaries
file_bound_start << MeshData.boundaries, 0

# ----------------------- #
#   FORM OPTIMIERUNG      #
# ----------------------- #

print('Schritt 2/2: Formoptimierung \n')

print("initial meshdistance: {0:8e} \n".format(bib.mesh_distance(MeshData, targetMeshData)))

# file_mesh << MeshData.mesh

# BFGS-memory initialisieren
bfgs_memory   = bib.bfgs_memory(np.zeros([memory_length, 2 * MeshData.mesh.num_vertices()]),
                                np.zeros([memory_length, 2 * MeshData.mesh.num_vertices()]), memory_length, 0)

# Lame-Parameter berechnen und speichern
mu_elas = bib.calc_lame_par(MeshData, mu_min, mu_max)

# Zaehler der Iterationen
counter = 0

# Startzeit zum Zeitmessen
start_time = time.time()

# Outputgraphendateien erstellen
file_output_grad = open(os.path.join(outputfolder, 'outputdata_grad.txt'),'a')
file_output_mshd = open(os.path.join(outputfolder, 'outputdata_mshd.txt'),'a')

# Formableitungen fuer Curv. Cond; 0 aktuelles Gitter, 1 vorheriges Gitter
curv_cond = np.zeros(2)

# Werte, welche aus voriger Iteration gespeichert sind; 0 aktuelles Gitter, 1 vorheriges Gitter
# Werte, ausser nrm_f_elas, werden in Schritt k nach printen von Schritt k-1 berechnet, aber
# erst in Schritt k+1 geprintet, deshalb keine Paare
nrm_f_elas = np.zeros(2)
J = 0.
nrm_U_mag = 0.
mesh_dist = 0.

# Berechnete Deformation s_k
last_defo = np.zeros([ 2 * MeshData.mesh.num_vertices()])

# Dummy um Schleife zu durchlaufen
nrm_f_elas[1] = 1.0


# Werte, welche fuer jeden Schritt geprintet werden
if(L_BFGS):
    print("Iteration " + "  ||f_elas||_L2 " + "  J = j + j_reg " + "  ||U||_L2  " + " Curv. Cond. " + " Meshdistance")
elif(L_BFGS == False):
    print("Iteration " + "  ||f_elas||_L2 " + "  J = j + j_reg " + "  ||U||_L2  " +  " Meshdistance")

# Start der Optimierungsschritte
while nrm_f_elas[1] > tol_shopt:

    # Zaehler erhoehen und ausgeben in Konsole
    counter += 1

    # ----------------------- #
    #      INTERPOLATION      #
    # ----------------------- #

    V = FunctionSpace(MeshData.mesh, "P", 1)
    z = project(y_z, V)
    mu_elas_projected = project(mu_elas, V)
    mu_elas_projected.rename("lame_par", "label")

    # --------------------------------------- #
    #  STATE, ADJOINT  & GRADIENT CALCULATION #
    # --------------------------------------- #

    # loese Zustandgleichung
    y = bib.solve_state(MeshData, f_values)

    # loese adjungierte Gleichung
    p = bib.solve_adjoint(MeshData, y, z)

    # Loese lineare Elastizitaetsgleichung, weise aktuelle nrm_f_elas zu
    U , nrm_f_elas[0] = bib.solve_linelas(MeshData, p, y, z, f_values, mu_elas_projected, nu)
    if(L_BFGS == False):
        # Ohne BFGS: Deformation ist negativer Gradient
        nrm_f_elas[1] = nrm_f_elas[0]
        S = Function(VectorFunctionSpace(MeshData.mesh, "P", 1, dim=2))
        S.vector()[:] = -1. * U.vector()

    #if(counter == 1): bfgs_memory.update_grad(U.vector().get_local())

    # ----------------------- #
    #       L-BFGS METHOD     #
    # ----------------------- #

    if(L_BFGS):

        # --------------------------------#
        #  CURVATURE CONDITION & PRINTING #
        # --------------------------------#

        # Curvature Condition berechnen
        V_k_1 = VectorFunctionSpace(MeshData.mesh, "P", 1, dim=2)
        last_defo_function = Function(V_k_1)

        if(bfgs_memory.step_nr == 0):
            last_defo_function.vector()[:] = U.vector().get_local()
            curv_cond[0] = bib.shape_deriv(MeshData, p, y, z, f_values, nu, last_defo_function)
        #last_defo_function.vector()[:] = bfgs_memory.deformation[0]
        if(bfgs_memory.step_nr > 0):
            last_defo_function.vector()[:] = last_defo
            curv_cond[0] = bib.shape_deriv(MeshData, p, y, z, f_values, nu, last_defo_function)

        # Curvature condition auswerten
        step_valid = True
        curv_cond_val = 0.
        if(bfgs_memory.step_nr > 0):
            curv_cond_val = curv_cond[1] - curv_cond[0]

        # Iterationsfortschritt in Konsole ausgeben, ein Zeitschritt versetzt, da curv_cond_val erst dann
        # berechnet werden kann
        if(counter > 1):
            print("   {0:3d}   ".format(counter-1) + " | " +
                  "   {0:.2E}  ".format(nrm_f_elas[1]) + " | " +
                  "   {0:.2E}  ".format(J) + " | " +
                  " {0:.2E}".format(nrm_U_mag) + " | " +
                  " {0:.2E}".format(curv_cond_val) + " | " +
                  " {0:.2E}".format(mesh_dist))

            # speichere in CSV-Datei (hier, da curv_cond_val berechnet werden musste) von vorigem Schritt
            file_csv.write(str(counter-1) + ";" +
                           str(nrm_f_elas[1]) + ";" +
                           str(J) + ";" +
                           str(nrm_U_mag) + ";" +
                           str(curv_cond_val) + ";" +
                           str(mesh_dist) + ";" + "\n")

        # aktualisiere, nachdem alte Werte geprintet wurden
        nrm_f_elas[1] = nrm_f_elas[0]

        # ----------------------- #
        #     NON-UPDATE STEP     #
        # ----------------------- #

        # Falls Curvature Condition nicht erfuellt, so update nicht und mache Schritt mit alter Information
        if(bfgs_memory.step_nr > 0):
            if(curv_cond_val <= 0.):
                print("Curvature Condition not fullfilled! Last step was invalid!"
                      "Programmiere das Mesh zurueckzusetzen? + break option")
                step_valid = False


        if(step_valid == False):
            print("break ausgeschaltet. Mache Schritt, aber kein Update!")
            #break
            # baue moeglichkeit zur auswahl, ob abbruch oder lieber weiterverwendung alter werte
            #bfgs_memory.update_grad(U.vector().get_local())
            if(bfgs_memory.step_nr == 0):
                # falls direkt Curvature condition verletzt, mache Gradientenschritt
                bfgs_memory.update_grad(U.vector().get_local())
                S = U
                S.vector()[:] = -1. * U.vector().get_local()
                bfgs_memory.update_defo(S.vector().get_local())
                last_defo = S.vector().get_local()

            if(bfgs_memory.step_nr > 0):

                # last_defo kann nicht bei aussetzen des updates zur berechnung der curv cond verwendet werden
                # ergo fuehre einen counter ein, s.d. man dann diesen nicht berechnet ( oder geht was besseres?)

                #bfgs_memory.update_grad(U.vector().get_local())
                #bfgs_memory.update_defo(last_defo)
                #S = bib.bfgs_step(MeshData, bfgs_memory, mu_elas_projected, bfgs_memory.gradient[0])

                S = bib.bfgs_step(MeshData, bfgs_memory, mu_elas_projected, U.vector().get_local())
                last_defo = S.vector().get_local()
                bfgs_memory.step_nr = bfgs_memory.step_nr - 1



        # ----------------------- #
        #       UPDATE STEP       #
        # ----------------------- #
        # L-BFGS Schritt
        elif(step_valid):
            bfgs_memory.update_grad(U.vector().get_local())
            if(bfgs_memory.step_nr == 0):
                S = U
                S.vector()[:] = -1. * U.vector().get_local()
                bfgs_memory.update_defo(S.vector().get_local())
                last_defo = S.vector().get_local()

            if(bfgs_memory.step_nr > 0):
                bfgs_memory.update_defo(last_defo)
                S = bib.bfgs_step(MeshData, bfgs_memory, mu_elas_projected, bfgs_memory.gradient[0])
                last_defo = S.vector().get_local()

        # L-BFGS Schritt erhoehen
        bfgs_memory.step_nr = bfgs_memory.step_nr + 1

        #curv_cond[1] = bib.shape_deriv(MeshData, p, y, z, f_values, nu, S)

        # TEST ZU SKALARPRODUKT UND ABLEITUNG; ECHTER UND FALSCHER GRADIENT
        #U_real, nrm_real = bib.solve_linelas(MeshData, p, y, z, f_values, mu_elas_projected, nu, zeroed = False)

        #loler = bib.shape_deriv(MeshData, p, y, z, f_values, nu, S)
        #loler2 = bib.shape_deriv(MeshData, p, y, z, f_values, nu, U_real)

        #print("shape derivative bfgs {0:.10e}".format(loler))
        #print("bilin_a {0:.10e}".format(bib.bilin_a(MeshData,S,-U_real,mu_elas_projected)))
        #print("shape derivative gradient {0:.10e}".format(loler2))
        #print("bilin_a {0:.10e}".format(bib.bilin_a(MeshData,U_real,-U_real,mu_elas_projected)))

    # ----------------------- #
    # BACKTRACKING-LINESEARCH #
    # ----------------------- #

    #print("Achtung: Wert des deformierten Gitters {0:.2E}  ".format(Jvalue))
    #print("Achtung: Wert der Ableitung {0:.2E}  ".format(bib.shape_deriv(MeshData, p, y, z, f_values, nu, S)))
    #print("Achtung: Wert der bilin {0:.2E}  ".format(bib.bilin_a(MeshData, U, S, mu_elas_projected)))
    if(Linesearch):

        scale_parameter = start_scale

        zero_function = Function(VectorFunctionSpace(MeshData.mesh, "P", 1, dim=2))
        current_value = bib.targetfunction(MeshData, zero_function, y_z, f_values, nu)
        S.vector()[:] = start_scale * S.vector()
        current_deriv = bib.shape_deriv(MeshData, p, y, z, f_values, nu, S)

        counterer = 0
        #print(current_value)
        #print(current_deriv)
        #U_real, trash = bib.solve_linelas(MeshData, p, y, z , f_values, mu_elas_projected, nu, zeroed = False)
        #print(bib.bilin_a(MeshData, U_real, U, mu_elas_projected))
        #print(Resetcounter)

        # Skaliert das Deformationsfeld bei jedem Schritt automatisch dauerhaft: NOCH OHNE AMIJO
        #while(bib.targetfunction(MeshData, S, y_z, f_values, nu) > current_value + c*scale_parameter*current_deriv):
        while(bib.targetfunction(MeshData, S, y_z, f_values, nu) >= current_value):
            print("shape deriv: {0:3e}".format(bib.shape_deriv(MeshData, p, y, z, f_values, nu, S)))
            #print("current Value: {0:3e}".format(current_value))
            #print("next Value: {0:3e}\n".format(bib.targetfunction(MeshData, S, y_z, f_values, nu)))
            #print("bilin deriv: {0:3e}".format(bib.bilin_a(MeshData,S,U_real,mu_elas_projected)))

            scale_parameter = shrinkage * scale_parameter
            S.vector()[:]   = shrinkage * S.vector()

            counterer = counterer + 1

            # ----------------------- #
            #      MEMORY REBOOT      #
            # ----------------------- #

            if(counterer >= 20):
                print("Had to break, restarting L-BFGS!")

                bfgs_memory.gradient = np.zeros([memory_length, 2 * MeshData.mesh.num_vertices()])
                bfgs_memory.deformation = np.zeros([memory_length, 2 * MeshData.mesh.num_vertices()])
                bfgs_memory.step_nr = 0
                Resetcounter = Resetcounter + 1
                S.vector()[:] = np.zeros(2*MeshData.mesh.num_vertices())
                break

        if(counterer < 20):
            #print("hurray")
            Resetcounter = 0

        # reskalierte Deformation als letzte Deformation uebergeben
        #bfgs_memory.deformation[0] = S.vector().get_local()
        last_defo = S.vector().get_local()


    #print("zweiter Resetcounter:")EIN
    #print(Resetcounter)

    if(Resetcounter >= 2):
        # Falls nach Restart immernoch keine gute Reskalierung gefunden, so stoppe das Verfahren
        print("Reboot didn't help, quitting L-BFGS optimization!")
        break

    if(L_BFGS):
        # setze neuen alte Formableitung zur Berechnung von curv_cond_val im naechsten Schritt
        curv_cond[1] = bib.shape_deriv(MeshData, p, y, z, f_values, nu, S)

    # ------------------------------ #
    #  MISC. CALCULATIONS & PRINTING #
    # ------------------------------ #

    # Norm des Deformationsvektorfeldes berechnen
    V           = FunctionSpace(MeshData.mesh, 'P', 1)
    U_magnitude = sqrt(dot(U, U))
    U_magnitude = project(U_magnitude, V)
    nrm_U_mag   = norm(U_magnitude, 'L2', MeshData.mesh)

    # Zielfunktional berechnen: KANN ERSETZT WERDEN AUS BIB
    zero_function = Function(VectorFunctionSpace(MeshData.mesh, "P", 1, dim=2))
    J = bib.targetfunction(MeshData, zero_function, y_z, f_values, nu)

    # Abstand zum Zielgitter berechnen
    mesh_dist = bib.mesh_distance(MeshData, targetMeshData)

    # Iterationsfortschritt in Konsole ausgeben, falls Gradientenverfahren genutzt wird
    # hier ist ein versetzen nicht noetig, da curv_cond_val nicht berechnet werden muss
    if(L_BFGS == False):
        print("   {0:3d}   ".format(counter) + " | " +
              "   {0:.2E}  ".format(nrm_f_elas[0]) + " | " +
              "   {0:.2E}  ".format(J) + " | " +
              " {0:.2E}".format(nrm_U_mag) + " | " +
              " {0:.2E}".format(mesh_dist))

        # Iterationsschritt in csv-Datei speichern
        file_csv.write(str(counter) + ";" +
                       str(nrm_f_elas[0]) + ";" +
                       str(J) + ";" +
                       str(nrm_U_mag) + ";" +
                       str(mesh_dist) + ";" + "\n")

    # Norm von aktuellem f_elas plotten
    file_output_grad.write(str(counter)       +";"+
                      str(np.log(nrm_f_elas[0]))  + "\n")

    # Abstand zur Zielform plotten
    file_output_mshd.write(str(counter)       +";"+
                      str(np.log(mesh_dist))  + "\n")

    # ----------------------- #
    #      DEFORMATION        #
    # ----------------------- #


    # KANN MAN EIGENTLICH ZU EINEM MACHEN, FALLS GRADIENTENVERFAHREN ZUSAMMEN MIT LINESEARCH GUT FUNKTIONIERT


    if(L_BFGS):
        # L-BFGS-Verfahren
        ALE.move(MeshData.mesh, S)

    elif(L_BFGS == False):
        # Gradientenverfahren
        #S = U
        #S.vector()[:] = -1.*U.vector().get_local()
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

if(L_BFGS):
    # hole bei L-BFGS das versetzte Printen des letzten Schrittes vor Ausstieg nach
    counter = counter + 1

    print("   {0:3d}   ".format(counter - 1) + " | " +
          "   {0:.2E}  ".format(nrm_f_elas[1]) + " | " +
          "   {0:.2E}  ".format(J) + " | " +
          " {0:.2E}".format(nrm_U_mag) + " | " +
          " {0:.2E}".format(curv_cond_val) + " | " +
          " {0:.2E}".format(mesh_dist))

    # speichere in CSV-Datei (hier, da curv_cond_val berechnet werden musste) von vorigem Schritt
    file_csv.write(str(counter - 1) + ";" +
                   str(nrm_f_elas[1]) + ";" +
                   str(J) + ";" +
                   str(nrm_U_mag) + ";" +
                   str(curv_cond_val) + ";" +
                   str(mesh_dist) + ";" + "\n")

    # Speichern in pvd-Datei fuer ParaView
file_mesh      << MeshData.mesh
file_sol_y     << y
file_adj_p     << p
file_def_u     << U
file_lame      << mu_elas_projected
file_bound_end << MeshData.boundaries

file_output_grad.close()
file_output_mshd.close()

# Zeitmessung
elapsed_time = time.time() - start_time

# Dauer der Berechnung in csv-Datei speichern
file_csv.write("Dauer der Berechnung in Sekunden;" +
               str(elapsed_time) + ";" + "\n")

# Datei schliessen
file_csv.close()

# Abschluss Ausgabe in Konsole
print("\n-------------------Fertig------------------------\n")
print("Iterationen:            {0:0d}".format(counter))
print("Dauer der Berechnungen: {0:.0f} min".format(int(elapsed_time)/60))
print("Ergebnisse in:          '{0:s}/'".format(outputfolder))
print("\n-------------------------------------------------\n")

# ----------------------- #
#     GRAPH PRINTING      #
# ----------------------- #


with open(os.path.join(outputfolder, 'outputdata_grad.txt')) as output:
    items  = (map(float, line.split(";")) for line in output)
    xs_grad, ys_grad = zip(*items)

plot1, = plt.plot(xs_grad, ys_grad, label="log Norm f_elas")

#legend_grad = plt.legend(handles=[plot1], loc=1)

with open(os.path.join(outputfolder, 'outputdata_mshd.txt')) as output:
    items  = (map(float, line.split(";")) for line in output)
    xs_mshd, ys_mshd = zip(*items)

plot2, = plt.plot(xs_mshd, ys_mshd, label = "log Meshdistance")

legend_mshd = plt.legend(handles=[plot1, plot2], loc=1)

plt.savefig(os.path.join(outputfolder, 'convergence_fig'))
plt.show()
