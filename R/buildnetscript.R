rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc() #free up memory and report the memory usage.

args <- commandArgs(trailingOnly = TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}
filename <- paste(paste("SubNx_",args[1],sep=""),".paj",sep="")
print(paste("Processing filename",filename))


library(network)
library(GERGM)

loadNetwork <- function(fileName){
  net <- read.paj(fileName)
  print(class(net))
  net %v% "rft" <- net %v% "5"
  net %v% "defect" <- net %v% "7"
  net %v% "cv" <- net %v% "9"
  net %v% "con" <- net %v% "11"
  net %v% "vm" <- net %v% "13"
  net %v% "pm" <- net %v% "15"
  edgeAttrName = list.edge.attributes(net)[2]
  net %e% "crossCorr" <- net %e% edgeAttrName
  
  delete.vertex.attribute(net,"4")
  delete.vertex.attribute(net,"5")
  delete.vertex.attribute(net,"6")
  delete.vertex.attribute(net,"7")
  delete.vertex.attribute(net,"8")
  delete.vertex.attribute(net,"9")
  delete.vertex.attribute(net,"10")
  delete.vertex.attribute(net,"11")
  delete.vertex.attribute(net,"12")
  delete.vertex.attribute(net,"13")
  delete.vertex.attribute(net,"14")
  delete.vertex.attribute(net,"15")
  delete.vertex.attribute(net,"shape")
  delete.edge.attribute(net,edgeAttrName)
  return (net)
}

plotNetwork <- function(net){
  net.ecol <- (net %e% "crossCorr"-min(net %e% "crossCorr"))*100
  net.vcol <- (net %v% "vm" - min(net %v% "vm"))*100
  plot(net,edge.col=net.ecol,vertex.col = net.vcol,coord=cbind(net %v% "x"*10,net %v% "y"*10))
}
targetDir <- c("../python/output/");
filename <- paste(targetDir,filename,sep="")

mynet <- loadNetwork(filename);

net<-as.matrix(mynet,matrix.type="adjacency");
network_covariate<-as.matrix(mynet,matrix.type="adjacency",attrname = "crossCorr");

#Covariance data
node_level_covariates <- data.frame(rft = mynet %v% "rft", cv = mynet %v% "cv", con = mynet %v% "con", pm = mynet %v% "pm");
#GERGM requires the network to be a matrix unlike ergm which works with a network object
# net is the matrix version of this
formula <- net ~ edges + netcov("network_covariate") + sender("con") + receiver("con") + out2stars("con") + in2stars("con") + ctriads  
calculate_network_statistics(net)

mynet.nmcc <- gergm(formula, 
          covariate_data = node_level_covariates,
          network_is_directed = TRUE,
          hyperparameter_optimization = TRUE,
          number_of_networks_to_simulate = 10000,
          thin = 1/100,
          proposal_variance = 0.2,
          MCMC_burnin = 5000,
          seed = 456,
          output_directory  = targetDir,
          output_name = basename(tools::file_path_sans_ext(filename)),     
          parallel = TRUE,
          cores = 2,       
          convergence_tolerance = 0.01,
          force_x_theta_update = 4);
fn <- paste(tools::file_path_sans_ext(filename),".RDS",sep="");
saveRDS(mynet.nmcc,fn);
