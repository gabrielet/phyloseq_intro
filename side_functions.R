# define plotting function using base R as vegan is mostly doing
# have a look here for some more ideas:
# https://fromthebottomoftheheap.net/2012/04/11/customising-vegans-ordination-plots/

# define function
plot_ordination <- function(data, metadata, var_fitting, scale_param, arrlen, maintitle) {
  
  # define shapes for depth legend
  shape_legend <- c(16, 17)
  # define shapes for depth dots. legend and dots shape must be in the same order
  shape <- c(21, 24)
  names(shape) <- c("topsoil", "deepsoil")
  
  # define colours for treatment
  treat_palette <- c("steelblue1", "firebrick1")
  names(treat_palette) <- c("C", "W")
  
  # adding species to triplot as mentioned here
  # https://github.com/vegandevs/vegan/issues/341
  
  # plot
  par(mar = c(2, 2, 2, 2) + 0.1)
  pl <- plot(data, type="none", cex.lab=1.5, cex.axis=1.5, cex.main=1.5, font=2, scaling=scale_param, main=maintitle)
  splen <- sqrt(rowSums(pl$species^2))
  #	with(metadata, text(pl, "species", select = splen > arrlen, arrow=TRUE, length=0.05))
  with(metadata, plot(var_fitting, add=T, col="black", cex=2, font=2, p.max=0.5))
  with(metadata, points(data, "sites", pch=shape[Depth], bg=treat_palette[Treatment], col="black", cex=3, scaling=scale_param))
  with(metadata, legend("topright", legend=levels(Treatment), fill=treat_palette, title="Treatment", cex=1.5))
  with(metadata, legend("bottomright", legend=levels(Depth), pch=shape_legend, title="Depth", cex=1.5))
  #with(metadata, text(data, cex=0.9, scaling=scale_param))
}
