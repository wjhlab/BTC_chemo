###############################################
#  do_Cellular_Influence_(CI)_quantification  #
###############################################


do_CI_quantification<- function(expr = dfgk, # should contain sample_id, X, Y position, celltype in the df
                                sample_id = sample_id,
                                kernels ="gaussian",
                                clusterdelevels = clusterdelevels,
                                sigma = 10){
  require(spatstat);
  
  spatwt<- c()
  for (i in sample_id){  
    expr_k<- expr[expr$sample_id==i, ]
    coord <- expr_k[, c("X_position", "Y_position")]
    
    rownames(coord)<- paste0(rownames(expr_k), ":", expr_k$celltype)
    
    gaussian<- data.frame(matrix(0, nrow = nrow(coord), ncol = length(clusterdelevels)))
    colnames(gaussian) <- clusterdelevels
    
    coord <- cbind(coord, gaussian)
    
    sapply(1:nrow(coord), function(r)
      coord[r, as.character(lapply(strsplit(as.character(rownames(coord)[r]), split = ":"),"[[",2))] <<- 1)
    
    # calculation of spatial weight 
    mywin <- owin(xrange = range(coord$X_position), yrange = range(coord$Y_position))
    spObj <- ppp(x=coord$X_position,y=coord$Y_position,window = mywin, marks = coord[,3:length(coord)])
    
    # adjust sigma 
    # leaveoneout = T -> removing the self
    K <- Smooth(spObj,kernel = kernels,at = "points",leaveoneout = T, sigma=sigma) 
    rownames(K)<-rownames(coord)
    K<- as.data.frame(K)
    K$sample_id <- i 
    
    spatwt<- rbind(spatwt, K)
  } 
  
  return (spatwt) 
}
  
