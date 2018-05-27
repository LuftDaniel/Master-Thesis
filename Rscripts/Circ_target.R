
#########################################
# ERSTER PLOT; CIRCLE VERFAHREN TARGETS #
#########################################

circle_gradient_target <- c(0.0027392312518675593,
                            0.0024321737839198285,
                            0.00216609858365016,
                            0.0019363686339932188,
                            0.0017279192047640764,
                            0.0015451085650981839,
                            0.0013908157764862428
                            )

circle_bfgs_target <- c(0.0027392312518675593,
                        0.0024321737839198285,
                        0.0006183524058802135,
                        0.0006114389863089016,
                        0.0005624478119306424,
                        0.0005624478119306424
                        )

circle_bfgs_target_fine <- c(0.0050020339613128795,
                             0.004366895286903954,
                             0.00018324673045683518,
                             7.60975106912989e-05,
                             2.3616979799604923e-05,
                             2.1998748836033196e-05,
                             2.191295116998402e-05
)

circle_gradient_target_fine <- c(0.0050020339613128795,
                                 0.004366895286903954,
                                 0.0037771025504782854,
                                 0.003239521554555098,
                                 0.002751875265986985,
                                 0.002319260464510079,
                                 0.0019390320730832149
)


plot(circle_bfgs_target_fine, type = "l", lty = 1, lwd = 1.2, xlab = "Schritt", ylab = "Zielfunktional J")
lines(circle_gradient_target_fine, type = "l", lty = 2, lwd = 1.2, col = "black")
lines(circle_bfgs_target, type = "l", lty = 3, lwd = 1.2, col = "black")
lines(circle_gradient_target, type = "l", lty = 4, lwd = 1.2, col = "black")
legend("topright", legend = c("L-BFGS; fein",
                              "Gradientenverfahren; fein",
                              "L-BFGS; grob",
                              "Gradientenverfahren; grob"),
       box.lty = 1, cex = 0.8, lty = c(1, 2, 3, 4), col = c("black", "black", "black", "black"), lwd = c(1.2, 1.2, 1.2, 1.2))


