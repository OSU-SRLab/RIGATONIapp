#' Filter graph data frames using user provided parameters
#'
#' @param Master A master dataframe from a RIGATONI output
#' @param MasterRNA RNA counts from a RIGATONI output
#' @param ControlRNA RNA counts from control samples within RIGATONI output
#' @param ControlCan Cancer type metadata for ControlRNA
#' @param alt.data Alteration dataframe from RIGATONI object
#' @param Cancers list of cancer types to analyze from the user
#' @param Mutations list of mutations to analyze from the user
#'
#' @return A list of dataframes for graph use
#' @export
filterResults <- function(Master, MasterRNA, ControlRNA, ControlCan, alt.data, Cancers = 'all', Mutations = 'all'){
  if (Cancers != 'all'){
    Master = Master[Master$Cancer %in% Cancers, ]
    MasterRNA = MasterRNA[, colnames(MasterRNA) %in% Master$CaseIDs]
    ControlCan = ControlCan[ControlCan[, 2] %in% Cancers, ]
    ControlRNA = ControlRNA[, colnames(ControlRNA) %in% ControlCan[, 1]]
  }
  if (Mutations != 'all'){
    Master = Master[Master$Alteration %in% Mutations, ]
    MasterRNA = MasterRNA[, colnames(MasterRNA) %in% Master$CaseIDs]
  }
  if (ncol(ControlRNA) < 500){
    rna = cbind(ControlRNA, MasterRNA)
    funs = Master[!(duplicated(Master$CaseIDs)), 1]
    funs = funs[match(colnames(MasterRNA), funs)]
    for (i in 1:nrow(Master)) {
      funs = replace(funs, which(funs == Master[i, 1]), Master[i, 6])
    }
    imm = Master[!(duplicated(Master$CaseIDs)), 1]
    imm = imm[match(colnames(MasterRNA), imm)]
    for (i in 1:nrow(Master)) {
      imm = replace(imm, which(imm == Master[i, 1]), Master[i, 7])
    }
    metadata1 = c(rep('Control', ncol(ControlRNA)), funs)
    metadata2 = c(rep('Control', ncol(ControlRNA)), imm)
  } else {
    rna = cbind(ControlRNA[, sample(1:ncol(ControlRNA), 500)], MasterRNA)
    funs = Master[!(duplicated(Master$CaseIDs)), 1]
    funs = funs[match(colnames(MasterRNA), funs)]
    for (i in 1:nrow(Master)) {
      funs = replace(funs, which(funs == Master[i, 1]), Master[i, 6])
    }
    imm = Master[!(duplicated(Master$CaseIDs)), 1]
    imm = imm[match(colnames(MasterRNA), imm)]
    for (i in 1:nrow(Master)) {
      imm = replace(imm, which(imm == Master[i, 1]), Master[i, 7])
    }
    metadata1 = c(rep('Control', 500), funs)
    metadata2 = c(rep('Control', 500), imm)
  }
  metadata = cbind(metadata1, metadata2, colnames(rna))
  metadata = as.data.frame(metadata)
  colnames(metadata) = c('Function', 'Immune', 'IDs')
  l = list(rna, metadata, Master, ControlRNA, MasterRNA)
  names(l) = c('rna', 'metadata', 'Master', 'ControlRNA', 'MasterRNA')
  return(l)
}

#' Recalcute alteration level function and immune status
#'
#' @param Master A master dataframe from a RIGATONI output
#'
#' @return new alt.data data frame
#' @export
getAltFunctionImmune <- function(Master){
  if (length(unique(Master$Alteration)) == 1){
    alts = unique(Master$Alteration)
  } else {
    master.split <- split(Master, Master$Alteration)
    master.split <- lapply(master.split, function(o){
      o = o[!(duplicated(o[, 1])), ]
      if (nrow(o) < 5){
        return(NA)
      } else {
        return(o)
      }
    })
    master.split = master.split[!(is.na(master.split))]
    master.split = master.split[names(master.split) != ""]
    alts = names(master.split)
  }
  nGOFs = c()
  nLOFs = c()
  nHots = c()
  nColds = c()
  totals = c()
  sigsF = c()
  sigsI = c()
  fun = c()
  im = c()
  for (x in alts){
    if (length(alts) == 1){
      tmp = Master
    } else {
      tmp = master.split[[x]]
    }
    total = nrow(tmp)
    nGOF = nrow(tmp[tmp$Function == 'GOF',])
    nLOF = nrow(tmp[tmp$Function == 'LOF',])
    nHot = nrow(tmp[tmp$Immune == 'Hot', ])
    nCold = nrow(tmp[tmp$Immune == 'Cold', ])
    z = stats::prop.test(nGOF, total, p = NULL, alternative = "two.sided", correct = TRUE)
    if (z$p.value < .05){
      if (nGOF > nLOF){
        nGOFs = c(nGOFs, nGOF)
        nLOFs = c(nLOFs, nLOF)
        totals = c(totals, total)
        sigsF = c(sigsF, z$p.value)
        fun = c(fun, 'GOF')
      } else {
        nGOFs = c(nGOFs, nGOF)
        nLOFs = c(nLOFs, nLOF)
        totals = c(totals, total)
        sigsF = c(sigsF, z$p.value)
        fun = c(fun, 'LOF')
      }
    } else {
      nGOFs = c(nGOFs, nGOF)
      nLOFs = c(nLOFs, nLOF)
      totals = c(totals, total)
      sigsF = c(sigsF, z$p.value)
      fun = c(fun, 'Unknown')
    }
    z = stats::prop.test(nHot, total, p = .05, alternative = "greater", correct = TRUE)
    if (z$p.value < .05){
      nHots = c(nHots, nHot)
      nColds = c(nColds, nCold)
      sigsI = c(sigsI, z$p.value)
      im = c(im, 'Hot')
    } else {
      z = stats::prop.test(nCold, total, p = .60, alternative = "greater", correct = TRUE)
      if (z$p.value < .05){
        nHots = c(nHots, nHot)
        nColds = c(nColds, nCold)
        sigsI = c(sigsI, z$p.value)
        im = c(im, 'Cold')
      } else {
        nHots = c(nHots, nHot)
        nColds = c(nColds, nCold)
        sigsI = c(sigsI, z$p.value)
        im = c(im, 'Unknown')
      }
    }
  }
  alt <- cbind(alts, totals, nGOFs, nLOFs, sigsF, fun, nHots, nColds, sigsI, im)
  alt = as.data.frame(alt)
  colnames(alt) = c('Alt.ID', 'Total.Samples', 'nGOFs', 'nLOFs', 'Func.Sig', 'Function', 'nHots', 'nColds', 'Immu.Sig', 'Immune')
  return(alt)
}

#' Get gene level expression boxplots
#'
#' @param rna rna for graphs from RIGATONI output
#' @param metadata metadata for graphs from RIGATONI output
#' @param gene gene provided by user to analyze
#' @param group X axis parameter
#'
#' @import ggplot2
#' @return ggplot2 object
#' @export
getExpressionBoxplot<- function(rna, metadata, gene, group){
  cd = rna[which(rownames(rna) == gene), ]
  cd = t(cd)
  cd = as.data.frame(cd)
  if (group != 'Mutational Status'){
    cd = cbind(metadata[, which(colnames(metadata) == group)], cd)
  } else {
    metadata[metadata != 'Control'] = 'Mutated'
    cd = cbind(metadata[, 1], cd)
  }
  colnames(cd) = c('Mut.Stat', 'Exp')
  gene_e_g = ggplot2::ggplot(cd, aes(x = Mut.Stat, y = log(round(Exp, 0), 2), fill = Mut.Stat))
  return(gene_e_g)
}

#' Get gene level expression boxplots
#'
#' @param quant quant data to use for analysis
#' @param metadata metadata for graphs from RIGATONI output
#' @param includeOther should the algorithm include "other" cell types in the result
#' @param group X axis parameter
#'
#' @import ggplot2
#' @import reshape2
#' @return ggplot2 object
#' @export
getQuantFigures <- function(quant, metadata, group, includeOther = 'Yes'){
  rnaquant = quant[quant[, 1] %in% metadata[, 3], ]
  rnaquant = rnaquant[match(metadata[, 3], rnaquant[, 1]), ]
  if (group != 'Mutational Status'){
    rnaquant = cbind(rnaquant, metadata[, which(colnames(metadata) == group)])
  } else {
    metadata[metadata != 'Control'] = 'Mutated'
    rnaquant = cbind(rnaquant, metadata[, 1])
  }
  colnames(rnaquant)[ncol(rnaquant)] = 'metadata'
  rnaquant = rnaquant[, 2:ncol(rnaquant)]
  quantm = reshape2::melt(rnaquant, id.vars = c('metadata'))
  if (includeOther != 'Yes'){
    quantm = quantm[quantm$variable != 'Other', ]
  }
  quantm$value = as.numeric(as.character(quantm$value))
  agquant = stats::aggregate(. ~ variable + metadata, quantm, mean)
  maov = stats::manova(as.matrix(rnaquant[, 2:(ncol(rnaquant)-1)]) ~ rnaquant$metadata, data = rnaquant)
  pval = summary(maov, tol=0)
  pval = pval$stats[1, 6]
  if (includeOther != 'Yes'){
    test = aggregate(value ~ metadata, agquant, sum)
    test = max(test[, 2])
    quant_g = ggplot2::ggplot(agquant, aes(fill = variable, x = metadata, y = value)) +
      ggplot2::geom_bar(position="stack", stat="identity") +
      ggplot2::annotate(geom="text", x=1, y=test+.05, label=paste0('Manova: ', round(pval, 3))) +
      ggplot2::theme_bw()
  } else {
    quant_g = ggplot2::ggplot(agquant, aes(fill = variable, x = metadata, y = value)) +
      ggplot2::geom_bar(position="stack", stat="identity") +
      ggplot2::annotate(geom="text", x=1, y=1.05, label=paste0('Manova: ', round(pval, 3))) +
      ggplot2::theme_bw()
  }
  l = list(quant_g)
  names(l) = 'ag_quant'
  for (c in unique(quantm$variable)){
    indiv = quantm[quantm$variable == c, ]
    gene_e_g = ggplot2::ggplot(indiv, aes(x = metadata, y = log(value, 2), fill = metadata))
    l[[c]] = gene_e_g
  }
  return(l)
}

#' Get cancer type distribution
#'
#' @param Master master dataframe from RIGATONI output
#' @param cancersFreq cancer type frequency data
#'
#' @import ggplot2
#' @return ggplot2 object
#' @export
getCancerTypeGraph <- function(Master, cancersFreq){
  Master = Master[!(duplicated(Master$CaseIDs)), ]
  cancers = table(Master$Cancer)
  cancers = as.data.frame(cancers)
  cancersFreq = cancersFreq[cancersFreq[, 1] %in% cancers[, 1], ]
  cancers = cancers[match(cancersFreq[, 1], cancers[, 1]),]
  cancers[, 2] = as.numeric(as.character(cancers[, 2]))/as.numeric(as.character(cancersFreq[, 2]))
  colnames(cancers) = c('Cancer', 'Prevalence')
  cancers = cancers[order(-cancers$Prevalence), ]
  cancers$Cancer = factor(cancers$Cancer, levels = cancers$Cancer)
  cancer_g = ggplot2::ggplot(cancers, aes(x=Cancer, y=Prevalence)) +
    ggplot2::geom_bar(stat = 'identity', aes(fill = Cancer)) +
    ggplot2::theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          legend.position = 'none')
  return(cancer_g)
}



