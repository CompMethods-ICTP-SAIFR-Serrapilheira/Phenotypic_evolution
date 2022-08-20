# Scientific programming course ICTP-Serrapilheira
# First version: August 18th 2022
# Make the plots of manduca simulations

# Load the package
require("gplots")

# Load the outputs
output.i <-read.table("./output/output.i_1_0_1000_.Rtab",sep=";",header=F, nrow = 600*801)
output.j <-read.table("./output/output.j_1_0_1000_.Rtab",sep=";",header=F, nrow = 600*801)
colnames(output.i)<-c("generation","id","ecological.trait","op","sex","species","father","mother")
colnames(output.j)<-c("generation","id","ecological.trait","op","sex","species","father","mother")

output.i2 <-read.table("./output/output.i_1_0.1_1000_.Rtab",sep=";",header=F, nrow = 600*801)
output.j2 <-read.table("./output/output.j_1_0.1_1000_.Rtab",sep=";",header=F, nrow = 600*801)
colnames(output.i2)<-c("generation","id","ecological.trait","op","sex","species","father","mother")
colnames(output.j2)<-c("generation","id","ecological.trait","op","sex","species","father","mother")

output.i3 <-read.table("./output/output.i_1_1_1000_.Rtab",sep=";",header=F, nrow = 600*801)
output.j3 <-read.table("./output/output.j_1_1_1000_.Rtab",sep=";",header=F, nrow = 600*801)
colnames(output.i3)<-c("generation","id","ecological.trait","op","sex","species","father","mother")
colnames(output.j3)<-c("generation","id","ecological.trait","op","sex","species","father","mother")


  # Defining colors
  colors.i<- densCols(output.i$ecological.trait, colramp = colorRampPalette(rich.colors(600)[-(1:3)]))
  colors.j<- densCols(output.j$ecological.trait, colramp = colorRampPalette(rich.colors(600)[-(1:3)]))
  colors.op<- densCols(output.j$op, colramp = colorRampPalette(rich.colors(600)[-(1:3)]))

  colors.i2<- densCols(output.i2$ecological.trait, colramp = colorRampPalette(rich.colors(600)[-(1:3)]))
  colors.j2<- densCols(output.j2$ecological.trait, colramp = colorRampPalette(rich.colors(600)[-(1:3)]))
  colors.op2<- densCols(output.j2$op, colramp = colorRampPalette(rich.colors(600)[-(1:3)]))

  colors.i3<- densCols(output.i3$ecological.trait, colramp = colorRampPalette(rich.colors(600)[-(1:3)]))
  colors.j3<- densCols(output.j3$ecological.trait, colramp = colorRampPalette(rich.colors(600)[-(1:3)]))
  colors.op3<- densCols(output.j3$op, colramp = colorRampPalette(rich.colors(600)[-(1:3)]))

  # Making plots
  par(mfrow = c(3, 3))
  plot(output.i$ecological.trait~output.i$generation, col=colors.i, cex=0.007,
       ylab="Ecological trait (zA)", xlab="Generation", ylim = c(-7,10))

  plot(output.j$ecological.trait~output.j$generation, col=colors.j, cex=0.007,
       ylab="Ecological trait (zB)", xlab="Generation", ylim = c(-7,10))

  plot(output.j$op~output.j$generation, col=colors.op, cex=0.007,
       ylab="Sexual trait (oB)",xlab="Generation", ylim = c(-7,10))
 #------------------------------------------------------

  plot(output.i2$ecological.trait~output.i2$generation, col=colors.i2, cex=0.007,
       ylab="Ecological trait (zA)", xlab="Generation", ylim = c(-6,15))

  plot(output.j2$ecological.trait~output.j2$generation, col=colors.j2, cex=0.007,
       ylab="Ecological trait (zB)", xlab="Generation", ylim = c(-6,15))

  plot(output.j2$op~output.j2$generation, col=colors.op2, cex=0.007,
       ylab="Sexual trait (oB)",xlab="Generation", ylim = c(-6,15))
 #------------------------------------------------------

  plot(output.i3$ecological.trait~output.i3$generation, col=colors.i3, cex=0.007,
       ylab="Ecological trait (zA)", xlab="Generation", ylim = c(-5,50))

  plot(output.j3$ecological.trait~output.j3$generation, col=colors.j3, cex=0.007,
       ylab="Ecological trait (zB)", xlab="Generation", ylim = c(-5,50))

  plot(output.j3$op~output.j3$generation, col=colors.op3, cex=0.007,
       ylab="Sexual trait (oB)",xlab="Generation", ylim = c(-5,50))

  par(mfrow = c(1, 1))
