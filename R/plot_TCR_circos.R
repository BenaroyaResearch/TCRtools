#' Make a circos plot of single cells, showing links between TCRs with identical junctions
#'
#' This function combines a number of steps for generating a circos plot from TCR sequences. The
#' purpose is to show clonotype sharing among cells within and between samples. Samples can optionally
#' be colored by grouping variables, and there are multiple options for determining TCR clonality.
#' 
#' @param tcr_cells data frame containing information at the cell level. Should include identifiers at minimum, plus optional columns for coloring the cells around the circos plot. Cells will be plotted in the order in which they are listed in this data frame.
#' @param tcr_links data frame with links to draw between cells. Should include identifiers for the two cells, plus optional columns for color and/or width of the links to be drawn. Often the output of \code{tabulate_shared_TCR_chains}.
#' @param ring_colors optional vector with number(s) or name(s) of the columns in \code{tcr_cells} to use for the colors of the ring. Current version can include up to 2 colors: list the inner ring first, outer ring second. 
#' @param link_colors optional number or name of the column in \code{tcr_links} to use for the colors of the links. Defaults to NULL, which yields black.
#' @param link_width optional number or name of the column in \code{tcr_links} to use for scaling the width of the links. Defaults to NULL, which results in all links having the same width.
#' @param filename the path to write the plot to. If provided, the function outputs a pdf of the plot, named "{filename}.pdf".
#' @param plottype character string, either "pdf" or "png". The type of plot object to output.
#' @param plotdims a numeric vector, the size (in inches) of the plotting object. Either the size of the pdf, or the size of the plotting window.
#' @export
#' 
#' @usage \code{
#' plot_TCR_circos(
#'      tcr_cells, tcr_links,
#'      ring_colors=NULL,
#'      link_colors=NULL, link_width=NULL,
#'      filename=NULL, plottype="pdf", plotdims=c(12,12))}
plot_TCR_circos <-
  function(
    tcr_cells, tcr_links,
    ring_colors=NULL,
    link_colors=NULL, link_width=NULL,
    filename=NULL, plottype="pdf", plotdims=c(12,12)) {
    
    if (!is.null(filename)) filename <- match.arg(filename, choices=c("pdf", "png"))
    
    # define limits for plotting
    n_cells <- nrow(tcr_cells)
    tcr_cells$xmin <- 0
    tcr_cells$xmax = n_cells
    
    #add sum values to tcr_cells, marking the x-position of the first links out (sum1) and in (sum2). Updated for further links in loop below.
    tcr_cells$sum <- length(unique(tcr_links$tcr1))
    
    
    if (!is.null(filename)) {
      if (plottype == "pdf") {
        pdf(file=filename, w=plotdims[1], h=plotdims[2])
      } else if (plottype == "png")
        png(file=filename, w=plotdims[1], h=plotdims[2], units="in")
      on.exit(dev.off(), add=TRUE) # close plotting device on exit
    } else dev.new(width=plotdims[1], height=plotdims[2])  ### open plotting window
    
    par(mar=rep(0,4))
    circos.clear()
    
    # basic circos graphic parameters
    circos.par(
      cell.padding=c(0,0,0,0), track.margin=c(0,0.15), start.degree = 90, gap.degree = 0)
    
    # sector details
    circos.initialize(factors = tcr_cells$tcr1, xlim = cbind(tcr_cells$xmin, tcr_cells$xmax))
    
    if (is.null(ring_colors)) {
      ring_colors <- "default"
      tcr_cells[,ring_colors] <- "black"
    }
      
    # plot sectors
    circos.trackPlotRegion(
      ylim = c(0, 1), factors = tcr_cells$tcr1, track.height=0.1,
      #panel.fun for each sector
      panel.fun = function(x, y) {
        #select details of current sector
        name = get.cell.meta.data("sector.index")
        i = get.cell.meta.data("sector.numeric.index")
        xlim = get.cell.meta.data("xlim")
        ylim = get.cell.meta.data("ylim")
        
        if (length(ring_colors)==1) {
          # single loop
          circos.rect(
            xleft=xlim[1], ybottom=ylim[1],
            xright=xlim[2], ytop=ylim[2],
            col = tcr_cells[i, ring_colors], border = tcr_cells[i, ring_colors])
        } else if (length(ring_colors)==2) {
          
          # colors for inner loop
          circos.rect(
            xleft=xlim[1], ybottom=ylim[1],
            xright=xlim[2], ytop=ylim[2]-0.5,
            col = tcr_cells[i, ring_colors[1]], border = tcr_cells[i, ring_colors[1]])
          
          # colors for outer loop
          circos.rect(
            xleft=xlim[1], ybottom=ylim[1]+0.5,
            xright=xlim[2], ytop=ylim[2],
            col = tcr_cells[i, ring_colors[2]], border = tcr_cells[i, ring_colors[2]])
        } else stop("Plotting with more than 2 loop colors is not currently supported.")
        
        # blank in part of main sector
        # circos.rect(xleft=xlim[1], ybottom=ylim[1], xright=xlim[2]-rowSums(m)[i], ytop=ylim[1]+0.3, col = "white", border = "white")
        
        #white line all the way around
        #circos.rect(xleft=xlim[1], ybottom=0.3, xright=xlim[2], ytop=0.32, col = "white", border = "white")
        
        #plot axis
        circos.axis(
          labels.cex=0.00000001, major.at=seq(from=0,to=floor(tcr_cells$xmax)[i],by=500), 
          labels.away.percentage = 0.15)
      })
    
    ### plot links
    if (is.null(link_colors)) {
      link_colors <- "default"
      tcr_links[,link_colors] <- "black"
    }
    
    if (is.null(link_width)) {
      link_width <- "default"
      tcr_links[,link_width] <- 1
    }
    for(k in 1:nrow(tcr_links)){
      # for(k in 1){
      # determine row of tcr_cells to use for coloring
      i <- match(tcr_links$tcr1[k], tcr_cells$tcr1)
      j <- match(tcr_links$tcr2[k], tcr_cells$tcr1)
      
      # draw links, colored by selected variable
      circos.link(
        sector.index1=tcr_cells$tcr1[i], point1=c(tcr_cells$sum[i]),
        sector.index2=tcr_cells$tcr1[j], point2=c(tcr_cells$sum[j]),
        col = tcr_cells[i, link_colors], rou1=0.75, rou2=0.75,
        lwd = tcr_links[k, link_width]*1.5)
    }
  }