output$Va[g+1] <- varA(pop,simParam = SP)[1,1]

output$Vg[g+1] <- varG(pop)[1,1]

output$Vp[g+1] <- varP(pop)[1,1]

output$mean_gv[g+1] <- mean(pop@gv[,1])

output$genicVa[g+1] <- genicVarA(pop)[1]

output$genicVg[g+1] <- genicVarG(pop)[1]

output$h2[g+1] = varA(pop)[1,1]/varP(pop)[1,1]