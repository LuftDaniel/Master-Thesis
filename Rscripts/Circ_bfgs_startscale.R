
###########################
# CIRCLE BFGS STARTSCALES #
###########################

circle_bfgs_target_fine <- c(0.0050020339613128795,
                        0.004366895286903954,
                        0.00018324673045683518,
                        7.60975106912989e-05,
                        2.3616979799604923e-05,
                        2.1998748836033196e-05,
                        2.191295116998402e-05
                        )

circle_bfgs_target_fine_highstartscale <- c(0.0050020339613128795,
                                       0.0026471180665470803,
                                       0.0015198920379566455,
                                       0.0002351076757080779,
                                       2.2048605819768815e-05,
                                       2.0493945816415698e-05,
                                       2.0421041476375052e-05
                                       )




plot(circle_bfgs_target_fine, type = "l", lty = 1, lwd = 1.2, xlab = "Schritt", ylab = "Zielfunktional J")
lines(circle_bfgs_target_fine_highstartscale, type = "l", lty = 2, lwd = 1.2, col = "black")
legend("topright", legend = c("L-BFGS; startscale 5",
                              "L-BFGS; startscale 20"),
       box.lty = 1, cex = 0.8, lty = c(1, 2), col = c("black", "black"), lwd = c(1.2, 1.2))

