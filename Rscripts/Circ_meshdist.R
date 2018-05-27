
#########################################
# ZWEITER PLOT; LOG MESHDISTANCES CIRCLE#
#########################################

circle_gradient_meshdist <- log(c(0.16323917917634417,
                              0.15687729649269982,
                              0.15083731060173977,
                              0.14512292196290075,
                              0.1397343785294206,
                              0.1346895393893655,
                              0.12998461916900939
))

circle_bfgs_meshdist <- log(c(0.16323917917634417,
                          0.15687729649269982,
                          0.07610595521928065,
                          0.07687091383766395,
                          0.07583622747922938,
                          0.07583622747922938
))

circle_bfgs_meshdist_fine <- log(c(0.15768581328359355,
                                      0.14614827967277305,
                                      0.023594121289522662,
                                      0.026891398485430495,
                                      0.016825091174837197,
                                      0.018003451756615562,
                                      0.01767085353289591
))

circle_gradient_meshdist_fine <- log(c(0.15768581328359355,
                                   0.14614827967277305,
                                   0.13494986490476446,
                                   0.12418610256324535,
                                   0.11393433125129507,
                                   0.10426368410541428,
                                   0.09521880143459831
))


plot(circle_bfgs_meshdist_fine, type = "l", lty = 1, lwd = 1.2, xlab = "Schritt", ylab = "log(Meshdist to target)")
lines(circle_gradient_meshdist_fine, type = "l", lty = 2, lwd = 1.2, col = "black")
lines(circle_bfgs_meshdist, type = "l", lty = 3, lwd = 1.2, col = "black")
lines(circle_gradient_meshdist, type = "l", lty = 4, lwd = 1.2, col = "black")
legend("topright", legend = c("L-BFGS; fein",
                              "Gradientenverfahren; fein",
                              "L-BFGS; grob",
                              "Gradientenverfahren; grob"),
       box.lty = 1, cex = 0.8, lty = c(1, 2, 3, 4), col = c("black", "black", "black", "black"), lwd = c(1.2, 1.2, 1.2, 1.2))

