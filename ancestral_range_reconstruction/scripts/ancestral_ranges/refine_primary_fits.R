# Tighten search and cross-check two optimizers; retain the highest observed likelihood.
args<-commandArgs(trailingOnly=TRUE);stopifnot(length(args)==1)
base<-normalizePath(args[1],winslash="/",mustWork=TRUE)
library(BioGeoBEARS);library(ape)
pilot<-file.path(base,"output/ancestral_ranges/refit_initial")
out<-file.path(base,"output/ancestral_ranges/refit_refined")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
tab<-read.csv(file.path(pilot,"model_summary.csv"))
areas<-c("EA","SA","NAM","EU")
states<-cladoRcpp::rcpp_areas_list_to_states_list(areas=areas,maxareas=4,include_null_range=TRUE)
labels<-vapply(states,function(s)if(anyNA(s))"empty" else paste(areas[s+1],collapse="+"),character(1))
summary<-list();checks<-list();nodeprobs<-list()
desc<-function(tr,n){if(n<=Ntip(tr))return(tr$tip.label[n]);unlist(lapply(tr$edge[tr$edge[,1]==n,2],function(k)desc(tr,k)),use.names=FALSE)}
for(i in seq_len(nrow(tab))) {
 r<-tab[i,];key<-paste(r$input,gsub("+","j",r$model,fixed=TRUE),sep="_")
 isbaseline<-grepl("_extant_only$",r$input)
 source<-if(isbaseline)file.path(base,"outputs/ancestral_area/baseline_multistart_001_v2",paste0("Walbeck_minimum_",gsub("+","j",r$model,fixed=TRUE),"_best_fit.rds")) else file.path(pilot,paste0(key,"_best_fit.rds"))
 best<-readRDS(source);initial<-best;oldlnL<-get_LnL_from_BioGeoBEARS_results_object(initial);bestlnL<-oldlnL
 for(method in c("bobyqa_full","L_BFGS_B")) {
  cache<-file.path(out,paste0(key,"_",method,".rds"))
  if(file.exists(cache))fit<-readRDS(cache) else {
   run<-initial$inputs;p<-run$BioGeoBEARS_model_object@params_table;opt<-initial$optim_result
   free<-rownames(p)[p$type=="free"]
   vals<-if(is.data.frame(opt))as.numeric(opt[1,seq_along(free)]) else opt$par
   stopifnot(length(vals)==length(free),all(is.finite(vals)))
   for(j in seq_along(free))p[free[j],c("init","est")]<-vals[j]
   run$BioGeoBEARS_model_object@params_table<-p;run$speedup<-FALSE;run$print_optim<-FALSE
   run$use_optimx<-method=="bobyqa_full"
   fit<-bears_optim_run(run);saveRDS(fit,cache)
  }
  lnL<-get_LnL_from_BioGeoBEARS_results_object(fit);stopifnot(is.finite(lnL))
  opt<-fit$optim_result
  code<-if(is.data.frame(opt))opt$convcode[1] else opt$convergence
  checks[[paste(key,method)]]<-data.frame(input=r$input,model=r$model,method=method,lnL=lnL,improvement_over_pilot=lnL-oldlnL,convergence_code=code,kkt1=if(is.data.frame(opt))opt$kkt1[1] else NA)
  if(lnL>bestlnL){best<-fit;bestlnL<-lnL}
  write.csv(do.call(rbind,checks),file.path(out,"optimizer_comparison.csv"),row.names=FALSE)
  writeLines(c("# 祖先区分析进度",paste("更新时间：",format(Sys.time(),tz="Asia/Shanghai")),"化石第一轮 84 次拟合已完成。",paste("优化交叉核查：",length(checks),"/",nrow(tab)*2,"完成；同时导出各节点概率供绘树。"),"随后检查年代、位置及后验树敏感性。"),file.path(base,"results/ancestral_area_progress.md"))
 }
 saveRDS(best,file.path(out,paste0(key,"_best_fit.rds")))
 tr<-read.tree(best$inputs$trfn)
 p<-best$ML_marginal_prob_each_state_at_branch_top_AT_node
 stopifnot(all(is.finite(p)),min(p)>-1e-10,max(abs(rowSums(p)-1))<1e-7)
 # Match homologous living-taxon clades, never numerical node labels across trees.
 for(n in (Ntip(tr)+1):(Ntip(tr)+tr$Nnode)) {
  living<-sort(desc(tr,n));living<-living[!grepl("^OBS_",living)]
  if(length(living)<2)next
  nodeprobs[[paste(key,n)]]<-data.frame(input=r$input,model=r$model,node=n,clade_key=paste(living,collapse="|"),area_state=labels,probability=p[n,])
 }
 summary[[key]]<-data.frame(input=r$input,model=r$model,lnL=bestlnL,k=r$k,AIC=2*r$k-2*bestlnL)
 t<-do.call(rbind,summary);t$deltaAIC<-ave(t$AIC,t$input,FUN=function(x)x-min(x));t$AIC_weight<-ave(exp(-.5*t$deltaAIC),t$input,FUN=function(x)x/sum(x))
 write.csv(t,file.path(out,"model_summary.csv"),row.names=FALSE)
 write.csv(do.call(rbind,nodeprobs),file.path(out,"node_probabilities.csv"),row.names=FALSE)
 cat("REFINED",key,bestlnL,"\n");flush.console()
}
writeLines("48 model/input combinations refined by two optimizers; highest observed likelihood retained.",file.path(out,"completed.txt"))
