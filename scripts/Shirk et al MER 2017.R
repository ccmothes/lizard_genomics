# load packages
library(MASS)
library(dplyr)
library(ResistanceGA)

# define file structures
project_dir = "C:/GD_methods"
CDPOP_dir = file.path(project_dir,"CDPOP")
LDmat_dir = file.path(project_dir,"LD_matrices")
GDmat_dir = file.path(project_dir,"GD_matrices")
results_dir = file.path(project_dir,"results")

# set variables
reps= c(1:10)
Disp_pmax = c(20,100)
metrics = c("Kc_Louiselle","Kc_Ritland","Rc_Queller","Rc_Lynch","Rc_Wang","Rc_Li","Fc_Lynch","Fc_Wang","Ra","Euc","DPS",
            "PCA_1ax","PCA_4ax","PCA_16ax","PCA_64ax")  



nhyp = 11 # number of hypotheses per batch

output_fn = "GD_methods_eval.csv"

params = expand.grid(pmax=Disp_pmax,metric=metrics,reps=reps,i=c(1:nhyp),n=c(200,1085),stringsAsFactors = F)

# create results file to append
header = data.frame("metric","Dispersal","n","r","transformation","true_model","true_model_index","topmodel","topmodel_index","pca_axes")

if(file.exists(file.path(results_dir,output_fn))==TRUE)
{
  results = read.csv(file.path(results_dir,output_fn),stringsAsFactors = F)
  params$combo = paste0(params[[1]],params[[2]],params[[3]],params[[4]],params[[5]])
  params_completed = paste0(results$Dispersal,results$metric,results$r,results$true_model_index,results$n)
  params_remaining = params[which(!(params$combo %in% params_completed)),]
} else {
  write.table(header,file.path(results_dir,output_fn),col.names = F,row.names = F,sep=",")
  params_remaining = params
}

IBD_LD_mat = file.path(LDmat_dir,"IBD_distmat.csv")


for(line in c(1:nrow(params_remaining)))
{
  pmax = params_remaining[line,"pmax"]
  metric = params_remaining[line,"metric"]
  r = params_remaining[line,"reps"]
  i = params_remaining[line,"i"]
  n = params_remaining[line,"n"]
  
  if(n=1085){
    sample_locs=c(1:1085)
  } else {
    set.seed(5468)
    sample_locs = sample(c(1:1085),200)
  }

  print(paste("dispersal:",pmax,"%max","n:",n))
  
  print(paste("   GD_metric:",metric))
  print(paste("      rep:",r))
  
  IBR_LD_mats = file.path(LDmat_dir,paste0("IBR_r",reps,"_distmat.csv"))
  LD_mats = c(IBR_LD_mats,IBD_LD_mat)
  
  IBR_genepop_files = file.path(CDPOP_dir,paste0("IBR_",pmax,"max"),paste0("batchrun",reps-1,"mcrun",r-1),"genepopgrid101.gen")
  IBD_genepop_file = file.path(CDPOP_dir,paste0("IBD_",pmax,"max"),paste0("batchrun0mcrun",r-1),"genepopgrid101.gen")
  genepop_files = c(IBR_genepop_files,IBD_genepop_file)
  
  IBR_grid_files = file.path(CDPOP_dir,paste0("IBR_",pmax,"max"),paste0("batchrun",reps-1,"mcrun",r-1),"grid101.csv")
  IBD_grid_file = file.path(CDPOP_dir,paste0("IBD_",pmax,"max"),paste0("batchrun0mcrun",r-1),"grid101.csv")
  grid_files = c(IBR_grid_files,IBD_grid_file)
  
  grid.file = read.csv(grid_files[i])
  grid.file$ID2 = ifelse(grid.file$ID == "OPEN",0,1)
  grid.file$ID2 = cumsum(grid.file$ID2)
  
  open_locs = which(grid.file$ID == "OPEN")
  
  valid_ind = sample_locs[!(sample_locs %in% open_locs)]
  GD_mat_ind = grid.file$ID2[valid_ind]
  
  xy = read.csv(file.path(points_dir,"CDPOP_xy.csv"))[valid_ind,c(4,2,3)]
  
  print(paste("         true model:",basename(LD_mats[i])))
  
  if(i==11)
  {
    GD_fn = file.path(GDmat_dir,paste0(metric,"_r",r,"_pmax",pmax,"_IBD.csv"))
  } else GD_fn = file.path(GDmat_dir,paste0(metric,"_r",r,"_pmax",pmax,"_IBR",i,".csv"))
  
  GD = as.matrix(read.csv(GD_fn,header=F))[GD_mat_ind,GD_mat_ind]
  
  # create an empty vector for AIC results
  pAIC = c()

  for(j in c(1:nhyp))
  {
    print(paste("            candidate model:",basename(LD_mats[j])))
    
    LD = log(as.matrix(read.csv(LD_mats[j],header=F))[valid_ind,valid_ind])
    diag(LD)=0
    LD = lower(LD)
    GD = lower(GD)
    pAIC = c(pAIC,AIC(MLPE.lmm(LD,GD,REML = FALSE)))
    rm(LD)
    
  } # close j loop
  
  top.model.index = which.min(pAIC)
  top.model = basename(LD_mats[top.model.index])
  
  output = data.frame(metric,pmax,n,r,"log",basename(LD_mats[i]),i,top.model,top.model.index,pca.axes)
  write.table(output,file.path(results_dir,output_fn),append=T,col.names = F,row.names = F,sep=",")
  
  rm(GD)
  
}  

results = read.csv(file.path(results_dir,output_fn))
results$correct = ifelse(results$true_model_index == results$topmodel_index,1,0)

method_perf = aggregate(correct~metric+n+pmax,results,mean)

