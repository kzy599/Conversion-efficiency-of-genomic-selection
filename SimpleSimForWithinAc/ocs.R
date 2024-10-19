G_A_ibd <- function(pop,pedmore){

    AnimalID <- pop@id
    
    ped_need = data.table(id = pop@id, sir = pop@father, dam = pop@mother)

    colnames(ped_need) = c("id","sir","dam")

    pedfinal <-rbind(pedmore,ped_need)
    pedfinal <-visPedigree:: tidyped(pedfinal,AnimalID)
    pedfinal<- pedfinal[,.(Ind,Sire,Dam)]

    Amatrix <- makeA(pedfinal)
    
    Amatrix <- as.matrix(Amatrix)

    tempBool_row = match(AnimalID,rownames(Amatrix))
    tempBool_col = match(AnimalID,colnames(Amatrix))
    
    AA<-Amatrix[tempBool_row,tempBool_col]
    
    AA <- cbind(pop@id,AA)
    
    colnames(AA) <- NULL
    
    rownames(AA) <- NULL
    
    return(AA)
}
get_os = function(){
  sysinf = Sys.info()
  if(!is.null(sysinf)){
    os = sysinf[['sysname']]
  }else{ # Sys.info not set up
    os = .Platform$OS.type
    if(grepl("^darwin", R.version$os)){
      os = "darwin"
    }else if(grepl("linux-gnu", R.version$os)){
      os = "linux"
    }
  }
  os = tolower(os)
  stopifnot(os=="windows" | os=="darwin" | os=="linux")
  return(os)
}

runProgram = function(progName, specFileName, runDir = "runDir") {
  #set up files
  dir.create(runDir,recursive = FALSE,showWarnings = FALSE)
  spec_file_path = file.path(runDir, specFileName)
  file.copy(specFileName, spec_file_path, overwrite = TRUE, recursive=FALSE)
  os = get_os()
  package_dir = find.package("AlphaLearn")
  if (os == "windows") {
    dll_dir = file.path(package_dir,os)
    flist <- list.files(dll_dir, "^(?i).+[.]dll$", full.names = TRUE)
    #  copy the DLLs
    file.copy(flist,paste(runDir,"\\.",sep=""),overwrite = TRUE, recursive=FALSE)
    prog_path = paste(package_dir,"\\",os,"\\",progName,".exe",sep="")
    destination = paste(file.path(runDir,progName),".exe",sep="")
    runCommand = paste(progName,".exe ",specFileName, sep="")
  } else {
    prog_path = paste(package_dir,"/",os,"/",progName,sep="")
    destination = paste(runDir,"/",progName,sep="")
    runCommand = paste(progName," ",specFileName,sep="")
  }
  file.copy(prog_path,destination,overwrite = TRUE, recursive=FALSE)
  old_dir = getwd()
  setwd(runDir)
  system(runCommand)
  setwd(old_dir)
}

runAlphaImpute = function(pedigree, genotypes) {
  runDir = "runDir"
  dir.create(runDir,recursive = FALSE,showWarnings = FALSE)
  ped_path = file.path(runDir,"pedigree.txt")
  write.table(pedigree,file=ped_path, row.names=F, col.names=F, quote=F)
  write.table(cbind(pedigree[,1], genotypes),file=paste(runDir,"/","genotypes.txt",sep=""), row.names=F, col.names=F, quote=F)
  
  spec <-ImputeParam$new()
  spec$write_out_spec(runDir)
  runProgram("AlphaImpute", "AlphaImputeSpec.txt",runDir)
  return(readInAI(path=runDir))
}

fixGenotypes = function(genotypes,pop) {
  ids = genotypes[,1]
  newLoci = sapply(pop@id, function(id) {
    which(ids == id)
  })
  genotypes = genotypes[newLoci,-1]
  return(genotypes)
}

readInAI = function(path=".") {
  genotypes <- as.matrix(fread(paste(path,"/Results/ImputeGenotypes.txt",sep="")))
  phase <- as.matrix(fread(paste(path,"/Results/ImputePhase.txt",sep="")))
  gdosages <- as.matrix(fread(paste(path,"/Results/ImputeGenotypeProbabilities.txt",sep="")))
  return(list(genotypes = genotypes, phase = phase, gdosages=gdosages))
}

#' @title Optimal Contribution Selection
#'
#' @description Perform Optimal Contribution Selection via AlphaMate
#' .
#' @param pop population
#' @param nCrosses number of matings/crosses
#' @param nFemalesMax maximum number of females
#' @param nMalesMax maximum number of   males
#' @param minFemaleContribution minimum number of matings/crosses per female
#' @param maxFemaleContribution maximum number of matings/crosses per female
#' @param minMaleContribution minimum number of matings/crosses per male
#' @param maxMaleContribution maximum number of matings/crosses per male
#' @param targetDegree targeted trigonometric degrees between genetic gain and group coancestry
#' @param targetCoancestryRate targeted rate of group coancestry
#' @param nProgenyPerCross number of progeny per mating/cross
#' @param use character specifiying, which type of criterion to use, either "pheno", "ebv", or "gv"
#'
#' @export
ocs = function(pop, nCrosses,usePed = FALSE, ped_parents = NULL,
               nFemalesMax=NULL, minFemaleContribution=NULL, maxFemaleContribution=NULL, equalizeFemaleContributions=NULL,
               nMalesMax=NULL, minMaleContribution=NULL, maxMaleContribution=NULL, equalizeMaleContributions=NULL,
               targetDegree=NULL, targetCoancestryRate=NULL,
               nProgenyPerCross, use) {
  # ---- Prepare data ----
  
  runDir = "runDir"
  dir.create(path=runDir, showWarnings=FALSE)
  # Criterion
  if (use == "gv") {
    tmp = data.frame(id=pop@id, crit=pop@gv)
  }
  if (use == "ebv") {
    tmp = data.frame(id=pop@id, crit=pop@ebv)
  }
  if (use == "pheno") {
    tmp = data.frame(id=pop@id, crit=pop@pheno)
  }
  write.table(x=tmp, file=file.path(runDir, "SelCriterion.txt"), col.names=FALSE, row.names=FALSE, quote=FALSE)
  
  
if(usePed){
  if(is.null(ped_parents)){stop("The ped of parents can't be NULL when usePed to perform OCS")}

  Amatrix <- G_A_ibd(pop = pop, pedmore = ped_parents)
  
  MASS::write.matrix(x=Amatrix, file=file.path(runDir, "GenotypeNrm.txt"))
  
  }else{
  G = AlphaMME::calcGIbs(X=pullSnpGeno(pop))
  
  MASS::write.matrix(x=cbind(pop@id, G), file=file.path(runDir, "GenotypeNrm.txt"))
  }
  
  
  # Gender
  tmp = data.frame(id=pop@id, genderRole=pop@sex)
  tmp$genderRole = as.numeric(factor(pop@sex, levels=c("M", "F")))
  write.table(x=tmp, file=file.path(runDir, "genderRole.txt"), col.names=FALSE, row.names=FALSE, quote=FALSE)
  
  # ---- Prepare spec file ----
  
  sink(file="AlphaMateSpec.txt")
  cat("NrmMatrixFile , GenotypeNrm.txt\n")
  if (use %in% c("gv", "ebv", "pheno")) {
    cat("SelCriterionFile , SelCriterion.txt\n")
  } else {
    stop("use muste be gv, ebv, or pheno")
  }
  cat("GenderFile , genderRole.txt\n")
  cat("NumberOfMatings , ", nCrosses, "\n")
  if (!is.null(nFemalesMax)) {
    cat("NumberOfFemaleParents , ", nFemalesMax, "\n")
  }
  if (!is.null(maxFemaleContribution) | !is.null(minFemaleContribution)) {
    cat("LimitFemaleContributions , Yes\n")
    if (!is.null(minFemaleContribution)) {
      cat("LimitFemaleContributionsMin , ", minFemaleContribution, "\n")
    }
    if (!is.null(maxFemaleContribution)) {
      cat("LimitFemaleContributionsMax , ", maxFemaleContribution, "\n")
    }
  }
  if (!is.null(equalizeFemaleContributions)) {
    cat("EqualizeFemaleContributions , Yes\n")
  }
  if (!is.null(nMalesMax)) {
    cat("NumberOfMaleParents , ", nMalesMax, "\n")
  }
  if (!is.null(maxMaleContribution) | !is.null(minMaleContribution)) {
    cat("LimitMaleContributions , Yes\n")
    if (!is.null(minMaleContribution)) {
      cat("LimitMaleContributionsMin , ", minMaleContribution, "\n")
    }
    if (!is.null(maxMaleContribution)) {
      cat("LimitMaleContributionsMax , ", maxMaleContribution, "\n")
    }
  }
  if (!is.null(equalizeMaleContributions)) {
    cat("EqualizeMaleContributions , Yes\n")
  }
  if (!is.null(targetCoancestryRate)) {
    cat("TargetCoancestryRate , ", targetCoancestryRate, "\n")
  }
  if (!is.null(targetDegree)) {
    cat("TargetDegree , ", targetDegree, "\n")
  }
  cat("Stop\n")
  sink()
  runProgram("AlphaMate", "AlphaMateSpec.txt", runDir=runDir)
  
  # ---- Get crossing plan ----
  
  crossPlan = read.table(file=file.path(runDir, "MatingPlanModeOptTarget1.txt"),
                         header=TRUE, colClasses="character")
  # crossPlan[[2]] = match(x=crossPlan[[2]], table=pop@id)
  # crossPlan[[3]] = match(x=crossPlan[[3]], table=pop@id)
  
  pedigree = matrix(data="", ncol=2, nrow=nrow(crossPlan) * nProgenyPerCross)
  k = 0
  for (i in 1:nrow(crossPlan)) {
    for (j in 1:nProgenyPerCross) {
      k = k + 1
      pedigree[k, 1] = crossPlan[i, 3]
      pedigree[k, 2] = crossPlan[i, 2]
    }
  }
  return(makeCross(pop=pop, crossPlan=pedigree))
}

