#dendro-arcs
#Assumes: 
#1) each variable represents a topic. meta-data should be included in the groupVars list - these variables are excluded from calculations
#2) the data includes a variable 'publications'. This is used to calculate proportions as shown in the colour-bands in the nodes
#3) nodes represent topics, links represent correlations

#other: The dendrogram is saved as a separate png file, which can be superimposed on the arcplot
#the function  'rolfsDendro2' sorts the variables. It is messy, and may need changing to suit your needs

Sys.setlocale("LC_CTYPE","russian")
load("sampleData.Rdata")
require(data.table)
require(ggdendro)
require(ggplot2)
#For installing arcdiagram
#require(devtools)
#install_github('arcdiagram', username='gastonstat')



#set up variables to be included/excluded
groupVars <- c("path","orientation","elections","date","publication","V1","predicted")
ignoreList=c()


dss <- sampleData

dendroArcs(dss,filename="out.pdf",nKeep=13,darker=F,titleTerm="'Egypt'")
#darker - if False correlations are relative to 1  (absolute correlation). If you see no or only weak links, try setting darker =T to emphasise weak connections. In this case links are relative to the highest observed correlation. 
#titleTerm - keyword linking texts for which correlations are observed. In the sample set texts about 'Egypt' are included. 
#filename= output file name
#Number of topics to include. Useful if you have hundreds of topics, but only want to visualise links between the most significant ones. 

zeroOne2 <- function(X) {(X - min(X))/diff(range(X))}


rolfsDendro2 <- function(dat,df3){
  require(ggdendro)
  vars <- as.hclust(dat)$order
  l <- as.hclust(dat)$label
  ttt <- colnames(df3)
  ttt <- ttt[!ttt %in% groupVars]
  #vars <- paste0("V",vars+1)
  vars <- ttt[vars]
  temp <- df3[,c(vars,"publication"),with=F]
  
  #assumes we are contrasting publications
  pub.names <- names(table(as.character(df3$publication)))
  m <- NULL
  for (i in pub.names){
    temp <- colMeans(df3[publication==i,vars,with=F])
    temp1 <- data.frame(temp)
    temp1$lab <- rownames(temp1)
    temp1$publication <- i
    m <- rbind(m,temp1)
    
  }
  temp <- m
  temp$lab <- factor(temp$lab, levels=unique(vars), ordered=TRUE)
  temp$temp <- (zeroOne2(temp$temp)/2.2) 

  dat=ggdendro:::dendrogram_data(dat)
  dat$segments$yend[dat$segments$yend<0.5] <- 0.5
  print(ggplot()+geom_segment(data = segment(dat), aes_string(x = "x",y = "y", xend = "xend", yend = "yend"))+
    theme_dendro()+ 
    coord_flip()+
    theme(axis.text.y = element_text( hjust = 1))+
    scale_x_discrete(labels = dat$labels$label)+
    geom_bar(data=temp,aes(lab,temp,fill=publication),stat="identity",position="dodge"))
  return(vars)
}

#Functions adapted from Gaston Sanchez:http://gastonsanchez.wordpress.com/2013/02/03/arc-diagrams-in-r-les-miserables/
dendroArcs <- function(df,filename,nKeep=500,darker=F,titleTerm="topic"){
  cols <- colnames(df)[!colnames(df) %in% groupVars]
  cs <- df[,lapply(.SD,sum),.SDcols=cols]
  cs <- sort(cs)
  keep <- tail(names(cs),nKeep)
  keep2 <- c(groupVars,keep)
  #colnames(df)[colnames(df) %in% keep2]
  df3 <- df[,c(colnames(df)[colnames(df) %in% keep2]),with=F]
  
  #correlations
  t <- cor(df3[,c(!colnames(df3) %in% groupVars),with=F])
  hc <- hclust(dist(t), "ward")
  vars <- colnames(df3)[!colnames(df3) %in% groupVars]

  hc2 <- as.dendrogram(hc)
  
  #reorder variables
  vars <- rolfsDendro2(hc2,df3)
  ggdendrogram(hc2)+ggsave(paste0(filename,".png"))
  #ADD ARCPLOT

  ct <- data.frame(t[vars,vars])
  labs <- vars
  #labs <- as.numeric(gsub("V","",vars))-1
  #labs <- topicsOld[labs]
  
  
  values=colSums(t)
  ord <- 1:length(vars)
  out <- NULL
  for (i in 1:nrow(ct)){
    out <- rbind(out,data.frame(rep(i,nrow(ct)),1:nrow(ct),ct[,i]))
  }
  out[out[,3]<0,3] <- 0
  out <- out[out[,3]>=0,]
  #out <- out[out[,3]>mean(out[,3])*2,]
  colnames(out) <- c("Source","Target","Weight")
  edgelist <- as.matrix(out[,1:2])
  
  #assumes we are contrasting newspapers, and that these are stored in a variable 'publications'
  pub.names <- names(table(as.character(df3$publication)))
  bands <- NULL
  for (i in pub.names){
      temp1<- sqrt(colSums(df3[publication==i,vars,with=F]))
      bands <- rbind(bands,temp1)
    }
  bands <- t(bands)
  colnames(bands) <- pub.names
  
  w2 <- out[,3]
  divideBy=max(w2)
  if (darker==T){
    w2[w2==1] <- 0
    divideBy=max(w2[w2<1])
  }
  print(w2/divideBy)
  cols = hsv(h=0, s=w2, v=0, alpha=0.5*zeroOne2(w2/divideBy))
  out[,3] <- zeroOne2(out[,3])
  require (RColorBrewer)
  col.bands = c("#4EA3CD", "#4E70CD", "#5E4ECD","#F2F2F2","#E5D9DA")
  col.bands=brewer.pal(5, "Set1")
  col.bands <- col.bands[1:length(pub.names)]
  
  pdf(filename, width=12, height=4)
  xmax <- arcBandBars(edgelist, bands,
              col.bands=col.bands, lwd=out[,3]*5, col=cols,
              col.terms="white", mar=c(1,1,3,1))
  title(c(paste0("Links between topics in texts about ",titleTerm), "Arc-diagram"), 
        cex.main=0.9, col.main="gray50")
  # add legend
  legend(x=xmax[2]-1.2, y=.65, title="Publication", text.col="gray25", cex=0.8,
         legend=c(paste0(pub.names)), pch=19, col=col.bands, bty="n")
  legend(x=xmax[2]-0.7, y=.65, text.col="gray25", cex=0.8,
         legend=c(paste0(1:length(labs)," ",labs)), bty="n")
  dev.off()
}

arcBandBars <- function(
  edgelist, bands, col.bands=NULL, sorted=TRUE, decreasing=FALSE,
  lwd=NULL, col=NULL, cex=NULL, col.nodes=NULL, cex.terms=NULL, col.terms=NULL,
  lend=1, ljoin=2, lmitre=1, bg=NULL, mar=c(4,1,3,1))
{
  # ARGUMENTS
  # edgelist:   two-column matrix with edges
  # bands:      numeric matrix with rows=nodes and columns=numbers
  # bars:       list of numeric tables with propotions for bar-charts
  # sorted:     logical to indicate if nodes should be sorted
  # decreasing: logical to indicate type of sorting (used only when sorted=TRUE)
  # lwd:        widths for the arcs (default 1)
  # col:        color for the arcs (default "gray50")
  # cex:        magnification of the nodes labels (default 1)
  # col.nodes:  color of the nodes labels (default "gray50")
  # cex.terms:  magnification of the terms in bar charts
  # col.terms:  color of the terms in bar charts
  # lend:       the line end style for the arcs (see par)
  # ljoin:      the line join style for the arcs (see par)
  # lmitre:     the line mitre limit fort the arcs (see par)
  # bg:         background color (default "white")
  # mar:        numeric vector for margins (see par)
  
  # make sure edgelist is a two-col matrix
  if (!is.matrix(edgelist) || ncol(edgelist)!=2)
    stop("argument 'edgelist' must be a two column matrix")
  edges = edgelist
  # how many edges
  ne = nrow(edges)
  # get nodes
  nodes = unique(as.vector(edges))
  nums = seq_along(nodes)
  # how many nodes
  nn = length(nodes)  
  # ennumerate
  if (sorted) {
    nodes = sort(nodes, decreasing=decreasing)
    nums = order(nodes, decreasing=decreasing)
  }
  # make sure bands is correct
  if (!is.matrix(bands) && !is.data.frame(bands))
    stop("argument 'bands' must be a numeric matrix or data frame")
  if (is.data.frame(bands))
    bands = as.matrix(bands)
  if (nrow(bands) != nn)
    stop("number of rows in 'bands' is different from number of nodes")
  
  # check default argument values
  if (is.null(lwd)) lwd = rep(1, ne)
  if (length(lwd) != ne) lwd = rep(lwd, length=ne)
  if (is.null(col)) col = rep("gray50", ne)
  if (length(col) != ne) col = rep(col, length=ne)
  if (is.null(col.nodes)) col.nodes = rep("gray50", nn)
  if (length(col.nodes) != nn) col.nodes = rep(col.nodes, length=nn)
  if (!is.null(cex) && length(cex) != nn) cex = rep(cex, length=nn)
  if (is.null(bg)) bg = "white"
  
  # nodes frequency from bands
  nf = rowSums(bands) / sum(bands)
  # words center coordinates
  fin = cumsum(nf)
  ini = c(0, cumsum(nf)[-nn])
  centers = (ini + fin) / 2
  names(centers) = nodes
  # node radiums
  nrads = nf / 2
  
  alt <- 0.9/length(centers)
  centers <- cumsum(rep(alt,length(centers)))
  centers <- centers*(3/max(centers))

  # arcs coordinates
  # matrix with numeric indices
  e_num = matrix(0, nrow(edges), ncol(edges))
  for (i in 1:nrow(edges))
  {
    e_num[i,1] = centers[which(nodes == edges[i,1])]
    e_num[i,2] = centers[which(nodes == edges[i,2])]
  }
  # max arc radius
  radios = abs(e_num[,1] - e_num[,2]) / 2
  max_radios = which(radios == max(radios))
  max_rad = unique(radios[max_radios] / 2)
  # arc locations
  locs = rowSums(e_num) / 2
  
  # function to get pie segments
  t2xy <- function(x1, y1, u, rad)
  {
    t2p <- pi * u + 0 * pi/180
    list(x2 = x1 + rad * cos(t2p), y2 = y1 + rad * sin(t2p))
  }
  
  # plot
  shrink=0.5
  par(mar = mar, bg=bg)
  plot.new()
  xlims=c(-0.025, max(centers)+.7)
  plot.window(xlim=xlims, ylim=c(-0.7*max_rad*shrink, 1*max_rad*2*shrink+0.1))
  # plot connecting arcs
  z = seq(0, pi, l=100)
  for (i in 1:ne)
  {
    radio = radios[i]
    x = locs[i] + radio * cos(z)
    y = 0.04+(radio * sin(z)) *shrink
    lines(x, y, col=col[i], lwd=lwd[i], 
          lend=lend, ljoin=ljoin, lmitre=lmitre)
  }
  # plot node bands
  constant=mean(nrads)*2
  print (constant)
  for (i in 1:nn)
  {
    radius = nrads[i]
    p = c(0, cumsum(bands[i,] / sum(bands[i,])))
    dp = diff(p)
    np = length(dp)
    angle <- rep(45, length.out = np)
    for (k in 1:np)
    {
      n <- max(2, floor(200 * dp[k]))
      r2=constant+(radius/2)
      P <- t2xy(centers[i], 0, seq.int(p[k], p[k+1], length.out=n), rad=r2)
      polygon(c(P$x2, centers[i]), c(P$y2, 0.005), angle=angle[i], 
              border=NA, col=col.bands[k], lty=0)
    }
    #draw white circles
    theta = seq(0, pi, length=100)
    #x3 = centers[i] + 0.7*nrads[i] * cos(theta)
    #y3 = 0 + 0.7*nrads[i] * sin(theta)
    x3 = centers[i] + (0.95-(nrads[i]*3))*constant * cos(theta)
    y3 = 0.00 + (0.95-(nrads[i]*3))*constant * sin(theta)
    polygon(x3, y3, col=bg, border=bg, lty=1, lwd=2)    
  }
  # add node names
  if (is.null(cex)) {
    cex = nf
    cex[nf < 0.01] = 0.01
    cex = cex * 5
  }
  # add node names
  text(centers, 0, nodes, cex=0.6, adj=c(0.5,0), col="black")
  return(xlims)
}

