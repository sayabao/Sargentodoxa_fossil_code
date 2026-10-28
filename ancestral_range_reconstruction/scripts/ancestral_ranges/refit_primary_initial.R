# Conditional genus-lineage observations; presence-only fossil geography.
args<-commandArgs(trailingOnly=TRUE);stopifnot(length(args)==1)
base<-normalizePath(args[1],winslash="/",mustWork=TRUE)
library(BioGeoBEARS);library(ape)
input<-file.path(base,"data/ancestral_ranges")
out<-file.path(base,"output/ancestral_ranges/refit_initial")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
trees<-file.path(input,"Walbeck_minimum_all_dated.nwk")
areas<-c("EA","SA","NAM","EU")
states<-cladoRcpp::rcpp_areas_list_to_states_list(areas=areas,maxareas=4,include_null_range=TRUE)
labels<-vapply(states,function(s)if(anyNA(s))"empty" else paste(areas[s+1],collapse="+"),character(1))
records<-list();summary<-list();roots<-list();errors<-list()
for(trfn in trees) {
 name<-sub("\\.nwk$","",basename(trfn))
 tree<-read.tree(trfn)
 for(model in c("DEC","DEC+J","DIVALIKE","DIVALIKE+J","BAYAREALIKE","BAYAREALIKE+J")) {
  key<-paste(name,gsub("+","j",model,fixed=TRUE),sep="_")
  result<-tryCatch({
   best<-NULL;bestlnL<--Inf
   for(start in 1:2) {
    cache<-file.path(out,paste0(key,"_start",start,".rds"))
    if(file.exists(cache)){fit<-readRDS(cache)} else {
     run<-define_BioGeoBEARS_run(trfn=trfn,geogfn=sub("\\.nwk$",".txt",trfn),max_range_size=4,include_null_range=TRUE,useAmbiguities=TRUE,min_branchlength=1e-6,num_cores_to_use=1,force_sparse=FALSE,return_condlikes_table=TRUE,calc_ancprobs=TRUE)
     run<-readfiles_BioGeoBEARS_run(run)
     # Verify '?' allows all supersets of observed areas, rather than creating absences.
     tipranges<-getranges_from_LagrangePHYLIP(run$geogfn)
     likes<-tipranges_to_tip_condlikes_of_data_on_each_state(tipranges,tree,states_list=states,maxareas=4,include_null_range=TRUE,useAmbiguities=TRUE)
     obs<-which(grepl("^OBS_",tree$tip.label))
     stopifnot(all(rowSums(likes[obs,,drop=FALSE])>=2))
     p<-run$BioGeoBEARS_model_object@params_table
     fixed<-function(p,key,value){p[key,"type"]<-"fixed";p[key,c("init","est")]<-value;p}
     if(grepl("DIVALIKE",model,fixed=TRUE)){p<-fixed(p,"s",0);p["ysv","type"]<-"2-j";p[c("ys","y","v"),"type"]<-"ysv*1/2";p<-fixed(p,"mx01v",0.5)}
     if(grepl("BAYAREALIKE",model,fixed=TRUE)){p<-fixed(p,"s",0);p<-fixed(p,"v",0);p["ysv","type"]<-"1-j";p["ys","type"]<-"ysv*1/1";p["y","type"]<-"1-j";p<-fixed(p,"mx01y",0.9999)}
     withj<-grepl("+J",model,fixed=TRUE)
     p["d",c("init","est")]<-c(.0001,.01)[start];p["e",c("init","est")]<-c(.0001,.001)[start]
     if(withj){p["j","type"]<-"free";p["j",c("init","est")]<-c(.001,.1)[start];p["j","max"]<-if(grepl("BAYAREALIKE",model,fixed=TRUE)).99999 else if(grepl("DIVALIKE",model,fixed=TRUE))1.99999 else 2.99999}
     run$BioGeoBEARS_model_object@params_table<-p
     check_BioGeoBEARS_run(run)
     fit<-bears_optim_run(run);saveRDS(fit,cache)
    }
    lnL<-get_LnL_from_BioGeoBEARS_results_object(fit);stopifnot(is.finite(lnL))
    records[[paste(key,start)]]<-data.frame(input=name,model=model,start=start,lnL=lnL,convcode=fit$optim_result$convcode[1],kkt1=fit$optim_result$kkt1[1],kkt2=fit$optim_result$kkt2[1])
    if(lnL>bestlnL){best<-fit;bestlnL<-lnL}
    write.csv(do.call(rbind,records),file.path(out,"start_comparison.csv"),row.names=FALSE)
   }
   saveRDS(best,file.path(out,paste0(key,"_best_fit.rds")))
   p<-best$ML_marginal_prob_each_state_at_branch_top_AT_node
   stopifnot(all(is.finite(p)),min(p)>-1e-10,max(abs(rowSums(p)-1))<1e-7)
   roots[[key]]<-data.frame(input=name,model=model,area_state=labels,probability=p[Ntip(tree)+1,])
   k<-if(grepl("+J",model,fixed=TRUE))3 else 2
   summary[[key]]<-data.frame(input=name,model=model,lnL=bestlnL,k=k,AIC=2*k-2*bestlnL,kkt1=best$optim_result$kkt1[1],kkt2=best$optim_result$kkt2[1])
   write.csv(do.call(rbind,roots),file.path(out,"root_probabilities.csv"),row.names=FALSE)
   tab<-do.call(rbind,summary);tab$deltaAIC<-ave(tab$AIC,tab$input,FUN=function(x)x-min(x));tab$AIC_weight<-ave(exp(-.5*tab$deltaAIC),tab$input,FUN=function(x)x/sum(x));write.csv(tab,file.path(out,"model_summary.csv"),row.names=FALSE)
   TRUE
  },error=function(e){errors[[key]]<<-conditionMessage(e);FALSE})
  writeLines(capture.output(str(errors)),file.path(out,"errors.txt"))
  cat("PROGRESS",key,result,"\n");flush.console()
 }
}
stopifnot(length(errors)==0,length(summary)==length(trees)*6)
writeLines("All planned fossil-observation pilot fits completed; interpretation remains conditional on genus-lineage scaffold and age choices.",file.path(out,"completed.txt"))
