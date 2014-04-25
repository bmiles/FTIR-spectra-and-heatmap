
require(ChemoSpec)
require(gplots)
require(RColorBrewer)
layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))

# Specify hooks for seperating data into groups
groupHooks <- c("A1", "B1", "C1","A2", "B2", "C2","A3", "B3", "C3")
groupColors <- brewer.pal(length(groupHooks), "Set1")

files2SpectraObject(gr.crit = groupHooks, gr.cols = groupColors,
                    freq.unit = "cm^-1",
                    int.unit = "AU",
                    descrip = "ATR-FTIR data",
                    format = "csv", 
                    alignTMS = FALSE,
                    out.file = "mydata",
                    debug = FALSE)

load("mydata.RData")
#baseline.irls(spectra, lambda1 = 5, lambda2 = 9, maxit = 200, wi = 0.05)
baseSpec <- baselineSpec(saveLoadReference, int=FALSE, retC=TRUE,lambda1 = 1, lambda2 = 9, maxit = 500, wi = 0.01, method="irls")
# remove noisy, uninteresting region:
newIR <- removeFreq(baseSpec, rem.freq =
                      baseSpec$freq > 1800 & baseSpec$freq < 2500)

irOffsetSpectra <- function () {
  plotSpectra(newIR, which=c(1,2,3,4,5,6,7,8,9), title="ATR-FTIR data for different sample coatings, 11-04-2014", offset= 0.2, yrange=c(-0.1,2), xaxt = "n")
  at <- seq(from = 0, to = max(newIR$freq), by = 100)
  
  # argument 'side': specifies which side of the plot the axis is to be drawn on.
  # 1 = x axis
  
  # argument 'labels'
  # "If labels is not specified, the numeric values supplied or calculated for 'at'
  # are converted to character strings as if they were a numeric vector printed
  # but note: "The code tries hard not to draw overlapping tick labels,
  # and so will omit labels where they would abut or overlap previously drawn labels"
  
  #axis(side = 1, at = at)
  
  # to make the labels fit, one can use the 'las' argument
  # "numeric in {0,1,2,3}; the style of axis labels"
  # where 2: always perpendicular to the axis
  
  #plot(Counts ~ Time, data = df, type = "l", xaxt = "n", xlab = "")
  axis(side = 1, at = at, las = 2, hadj = 0.9)
  mtext(text = "Absorbance", side = 2, line = 4)
  mtext(text = "wavenumber cm^-1", side = 1, line = 4)
}

# Heat Map, returns a heat map of intensities, takes a wavenumber or range (and calculates the mean)
irHeatmap <- function (wavenumber) {
  # Setup these parametric variables.
  pH <- c(3,7.3,12)
  collagenConc <- c(0,10,30)
  
  # Array to hold the intensities for a given wavenumber
  absVal <- c()

  
  # Loop gets intensities.
  for (i in c(1:length(newIR$groups))) {
    nextVal <-mean(c(newIR$data[i, wavenumber]))
    absVal <- c(absVal, nextVal)
  }
  
  # creates a matrix of the intensity values.
  vnames <-c()   
  for (i in c(1:length(newIR$names))) {
   nextName <-substring(c(newIR$names), 1, 2)
   vnames <- c(vnames, nextName)
  }
  mnames <- matrix(vnames, nrow = 3, ncol = 3, byrow = TRUE, dimnames = list(pH, collagenConc))

  mdat <- matrix(c(absVal), nrow = 3, ncol = 3, byrow = TRUE, dimnames = list(pH, collagenConc))
  colors <- brewer.pal(9, "YlGnBu")
  my_palette <- colorRampPalette(colors)(n = 299)
  heatmap.2(mdat, 
            main = "Collagen deposition dependence on pH and concentration", # heat map title
            cellnote = mnames,
            notecol="black",      # change font color of cell labels to black
            density.info="none",  # turns off density plot inside color legend
            trace="none",         # turns off trace lines inside the heat map
            margins =c(12,12),     # widens margins around plot
            col=my_palette,       # use on color palette defined earlier 
            #breaks=col_breaks,    # enable color transition at specified limits
            dendrogram="none",     # only draw a row dendrogram
            Colv="NA",            # turn off column clustering
            Rowv="NA",
            xlab = "Collagen Concentration v/v %",
            ylab = "pH",
            srtCol = 0,
            offsetRow = 2,
            offsetCol = 3)
}
irHeatmapRMS <- function (wavenumber) {
  # Setup these parametric variables.
  pH <- c(3,7.3,12)
  collagenConc <- c(0,10,30)
  
  # Array to hold the intensities for a given wavenumber
  absVal <- c()
  
  
  # Loop gets intensities.
  for (i in c(1:length(newIR$groups))) {
    y <- c(newIR$data[i, wavenumber])
    nextVal <-sqrt(sum(y^2)/length(y))
    absVal <- c(absVal, nextVal)
  }
  
  # creates a matrix of the intensity values.
  vnames <-c()   
  for (i in c(1:length(newIR$names))) {
    nextName <-substring(c(newIR$names), 1, 2)
    vnames <- c(vnames, nextName)
  }
  mnames <- matrix(vnames, nrow = 3, ncol = 3, byrow = TRUE, dimnames = list(pH, collagenConc))
  
  mdat <- matrix(c(absVal), nrow = 3, ncol = 3, byrow = TRUE, dimnames = list(pH, collagenConc))
  colors <- brewer.pal(9, "YlGnBu")
  my_palette <- colorRampPalette(colors)(n = 299)
  heatmap.2(mdat, 
            main = "Collagen deposition dependence on pH and concentration", # heat map title
            cellnote = mnames,
            notecol="black",      # change font color of cell labels to black
            density.info="none",  # turns off density plot inside color legend
            trace="none",         # turns off trace lines inside the heat map
            margins =c(12,12),     # widens margins around plot
            col=my_palette,       # use on color palette defined earlier 
            #breaks=col_breaks,    # enable color transition at specified limits
            dendrogram="none",     # only draw a row dendrogram
            Colv="NA",            # turn off column clustering
            Rowv="NA",
            xlab = "Collagen Concentration v/v %",
            ylab = "pH",
            srtCol = 0,
            offsetRow = 2,
            offsetCol = 3)
}
irOffsetSpectra()

