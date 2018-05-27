
##########################
# CIRCLE BFGS MEMLENGTHS #
##########################


circle_bfgs_target_fine_menlen3 <- c(0.0050020339613128795,
                                     0.004366895286903954,
                                     0.00018324673045683518,
                                     7.60975106912989e-05,
                                     2.3609103666971777e-05,
                                     2.1994310283399817e-05,
                                     2.1908675594888848e-05
                                     )

circle_bfgs_target_fine <- c(0.0050020339613128795,
                             0.004366895286903954,
                             0.00018324673045683518,
                             7.60975106912989e-05,
                             2.3616979799604923e-05,
                             2.1998748836033196e-05,
                             2.191295116998402e-05
)

circle_bfgs_target_fine_memlen2 <- c(0.0050020339613128795,
                                     0.004366895286903954,
                                     0.00018324673045683518,
                                     7.52830076495044e-05,
                                     2.3754789222989905e-05,
                                     2.221958362920098e-05,
                                     2.212663875917435e-05
                                     )


plot(circle_bfgs_target_fine, type = "l", lty = 1, lwd = 1.2, xlab = "Schritt", ylab = "Zielfunktional J")
lines(circle_bfgs_target_fine_menlen3, type = "l", lty = 2, lwd = 1.2, col = "black")
lines(circle_bfgs_target_fine_memlen2, type = "l", lty = 3, lwd = 1.2, col = "black")
legend("topright", legend = c("L-BFGS; mem.len. 60",
                              "L-BFGS; mem.len. 3",
                              "L-BFGS; mem.len. 2"),
       box.lty = 1, cex = 0.8, lty = c(1, 2, 3), col = c("black", "black", "black"), lwd = c(1.2, 1.2, 1.2))

