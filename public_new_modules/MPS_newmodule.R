# --
#-- conda activate MPS_html
#--- 


#--- Load all the library ---#

library('survminer')
library('survival')
library('shiny')
library('ggplot2')
library('gridExtra')
library('enrichplot')
library('grid')
library('clusterProfiler')
library('htmltools')
library('htmlwidgets')
library('rmarkdown')
library('foreach')
library('doMC')

organism = "org.Hs.eg.db"
library(org.Hs.eg.db)

#-- Library for rsconnet and google drive --#

library('rsconnect')
library('googledrive')
library('readr')
library('httr')
library('bslib')

#--- Google drive authentication ---#
# -- Run only one time to create .secrets file --#

#options(gargle_oauth_cache = '.secrets')

#check the value of the option, if you like
#gargle::gargle_oauth_cache()
#googledrive::drive_auth()
#list.files('.secrets/')

options(
  gargle_oauth_cache = '.secrets',
  gargle_oauth_email = 'tavazoie.ru.lab@gmail.com'
)

#drive_auth()


tcga_key = list('AML'='aml', 'bladder'='blca', 'breast'='brca', 'cervical'='cesc', 
                'colon'='coad', 'esophageal'='esca', 'GBM'='gbm', 'glioma'='glm', 
                'head_neck'='hnsc', 'kidney_cc'='kicc', 'kidney_ch' = 'kich', 
                'kidney_pa' = 'kipa', 'liver'='lihc', 'lung_ad'='luad',
                'lung_sq'='lusq', 'melanoma'='skcm','ovarian'='ovsc', 'pancreatic'='paad', 
                'paraganglioma'='pcpg', 'prostate'='prad', 'rectal'='read', 'sarcoma'='sarc', 
                'stomach'='stad', 'testicular'='tgct', 'thymoma'='thym', 'thyroid'='thca', 
                'uterine_endometrial'='ucec')

col_surv = rev(c(rgb(0.86, 0.2, 0.3, 0.75), rgb(0, 0.5, 1, 0.75)))
quant_th = 0.1 ; MI_bins = 10

#Define a function to read files from Google Drive
read_drive_file <- function(file_name) {
  
  file_type <- tools::file_ext(file_name)
  temp_file <- tempfile(fileext = paste0(".", file_type))
  
  drive_download(file_name, path = temp_file, overwrite = TRUE)
  
  data <- switch(file_type,
                 "csv" = read_csv(temp_file),
                 "txt" = read.table(temp_file, header = TRUE, sep = "\t", quote = "\"", comment.char = ''),
                 "rds" = readRDS(temp_file),
                 stop("Unsupported file type"))
  unlink(temp_file)
  return(data)
}


# Function to get gene info data from Google Drive
get_gene_info_from_drive <- function(file_name) {
  gene_info <- read_drive_file(file_name)
  gene_info$Entrez.Gene.ID <- as.character(gene_info$Entrez.ID)
  gene_info_filt <- gene_info[which(gene_info$Entrez.Gene.ID != '' | gene_info$Approved.Symbol != ''), grep('Entrez|Approved.Symbol', colnames(gene_info))]
  return(gene_info_filt)
}

gene_info_file_name <- "EntrezIDs_To_ApprovedSymbol_20221107.txt"
common_genes_set_file_name <- "common_TCGA_genes_rds.rds"
cmd_mkdir = "mkdir -p DIR"


cmd_mkdir = "mkdir -p DIR"



#------ 
# Functions
#-------

#--- Function 1: Calculate Mutual Information Values

calculateMI_v2 <- function(x, y) {
  rs_cMatrix = x + y ; num_genes = sum(rs_cMatrix)
  log_t1 = log((num_genes*x)/(sum(x) * rs_cMatrix)) ; log_t1[is.infinite(log_t1)] = 0
  log_t2 = log((num_genes*y)/(sum(y) * rs_cMatrix)) ; log_t2[is.infinite(log_t2)] = 0
  val_ret = sum((1/num_genes)* ((x)*log_t1 + (y)*log_t2))
}

#--- FUNCTION 2 : get survival Shiny

getSurv_shiny <- function(all_clin_df, surv_type = 'OS', rand_iter = 1000, samp_name = 'MPS+', ctrl_name = 'MPS-') {
  library('survival')
  all_clin = all_clin_df
  s_name = samp_name #; if(s_name == 'positive') { s_name = 'MPS+' }
  c_name = ctrl_name #; if(c_name == 'negative') { c_name = 'MPS–' }
  if (surv_type == 'OS') { all_clin$mod_censor_time = as.numeric(all_clin$OS.time) ; all_clin$death_event_binary = as.numeric(all_clin$OS) }
  if (surv_type =='PFI') { all_clin$mod_censor_time =as.numeric(all_clin$PFI.time) ; all_clin$death_event_binary = as.numeric(all_clin$PFI) }
  len_s = sum(all_clin$group == samp_name) ; len_c = sum(all_clin$group == ctrl_name)
  s_all = survfit(Surv(as.numeric(all_clin$mod_censor_time), all_clin$death_event_binary) ~ 1)
  ss = Surv(as.numeric(as.character(all_clin$mod_censor_time)),all_clin$death_event_binary)
  val_ret = NA
  pv = s_fit = pv_cox = hz_ratio = s_coxph = z = c_ix = NA
  if (! (sum(is.na(ss)) == length(is.na(ss)))) {
    s_fit = survfit(Surv(mod_censor_time, death_event_binary) ~ group, data = all_clin)
    s_diff = tryCatch({ 
      survdiff(Surv(mod_censor_time, death_event_binary) ~ group, data = all_clin)
    }, error = function(err) { list(chisq=NA, n=dim(all_clin)[1]) } )
    pv <- ifelse(is.na(s_diff),1,(round(1 - pchisq(s_diff$chisq, length(s_diff$n) - 1),50)))[[1]]
    
    s_coxph = coxph(Surv(mod_censor_time, death_event_binary) ~ group, data = all_clin)
    if (! is.na(s_coxph[['coefficients']])) { 
      s_c = summary(s_coxph)
      hz_ratio = as.numeric(data.frame(s_c$coefficients)['exp.coef.'])
      pv_cox = as.numeric(data.frame(s_c$coefficients)['Pr...z..'])
      z = as.numeric(data.frame(s_c$coefficients)['z'])
      c_ix = as.numeric(summary(s_coxph)$concordance['C'])
    }
    
    pv_r = c()
    
    for(r_i in seq(rand_iter)) {
      all_clin_r = all_clin ; prob_v = as.numeric(prop.table(table(all_clin$group)))
      all_clin_r$group = sample(c(s_name, c_name), dim(all_clin)[1], replace=T, prob=prob_v)
      s1_r = tryCatch({ 
        survdiff(Surv(mod_censor_time, death_event_binary) ~ group, data = all_clin_r)
      }, error = function(err) { list(chisq=NA, n=dim(all_clin_r)[1]) } )
      pv_r = c(pv_r, round(1 - pchisq(s1_r$chisq, length(s1_r$n) - 1),50))
    }
    fp = (sum(pv_r < pv)/rand_iter)
    median_surv = summary(s_fit)$table[,'median']
    med_ab = round(as.numeric(median_surv[paste('group=',s_name,sep='')]),1)
    med_re = round(as.numeric(median_surv[paste('group=',c_name,sep='')]),1)
    m_sub = paste(s_name,' ', len_s, '(', med_ab,' mo)',
                  '| ', c_name,' ', len_c, '(', med_re, ' mo)', sep='')
    
    out_df = data.frame(pval_surv = pv, pval_cox = pv_cox, C_index = c_ix, z_cox = z,HZR_cox = hz_ratio, info_survival = m_sub)
    
    val_ret = list(KM_pv = pv, KM_fit = s_fit, clin_data = all_clin, cox_pv = pv_cox, hzr = hz_ratio, 
                   cox_fit = s_coxph, avg_fit = s_all, rand_pv = pv_r, z=z,
                   concordance = c_ix, out_df = out_df, fdr = fp)
  }
}


#----- PART 1 : Survival Analysis using new module ----#

#--- FUNCTION 1 : MPS+ and MPS- with clinical data \# name : get_MPS_newModule for SHINY NOW

get_MPS_newModule <- function(input_genes, module_name = '', collection = 'tcga', select_cohorts = 'breast', number_bins = 10, pval_thresh = 0.01, number_rand = 10000, MPS_thresh = 0) {
  if (select_cohorts == 'full') { coh_list = names(tcga_key) }
  if (select_cohorts != 'full') { coh_list = intersect(select_cohorts, names(tcga_key)) }
  
  #cont_genes = disc_genes = data.frame(gene = as.character(gene2module$Approved.Symbol))
  T_1 = Sys.time()
  common_genes_set <- read_drive_file(common_genes_set_file_name)
  g_set = common_genes_set#readRDS(paste(parent_dir, 'shiny_test/data/common_TCGA_genes_rds', sep=''))
  inp_g = toupper(as.character(input_genes[,1]))
  
  #inp_g = as.character(read.csv(input_genes,sep='\t',header=T)[,1])
  genes_in_mod = intersect(inp_g, g_set) ; div_f = as.numeric(length(genes_in_mod))
  if (length(genes_in_mod) < 1) {
    stop("Module empty. Please select another module")
  }
  module_id = paste('id_', sample(1000,1), '_', gsub('-','',Sys.Date()),sep='') ; if (module_name == '') { module_name = module_id }
  set.seed(108)
  num_rand = number_rand
  total_genes = length(g_set)
  total_genes_bin = rmultinom(1, total_genes, rep((1/number_bins), number_bins))
  num_in_path = div_f
  path_genes_bin = rmultinom(num_rand, num_in_path, rep((1/number_bins), number_bins))
  c_vec_r = path_genes_bin ; nc_vec_r = total_genes_bin[,1] - c_vec_r
  v_rand = sapply(seq(num_rand), function(i) calculateMI_v2(c_vec_r[,i], nc_vec_r[,i]))
  mu_r = mean(v_rand) ; sd_r = sd(v_rand)
  v_rand_q = as.numeric(quantile(v_rand, (1-pval_thresh)))
  len_rand = length(v_rand) 
  #system(gsub('DIR', paste(parent_dir, 'tmp_new_module', sep=''), cmd_mkdir))
  #saveRDS(genes_in_mod, paste(parent_dir, 'tmp_new_module/input_genes', sep=''))
  all_new_clin = c()
  for(dis_ in coh_list) {
    d_t = tcga_key[[dis_]]
    
    con_g_name = paste0(d_t, '_quant_', quant_th, '_primary_zscore.rds')
    dis_g_name = paste0(d_t, '_quant_', quant_th, '_primary_zscore_bins', MI_bins,'.rds')
    #con_g = readRDS(paste(parent_dir, 'data/', 'disease_gene_expression/', collection, '/', d_t, '_quant_', quant_th, '_primary_zscore.rds', sep=''))
    #dis_g = readRDS(paste(parent_dir, 'data/', 'disease_gene_expression/', collection, '/', d_t, '_quant_', quant_th, '_primary_zscore_bins.rds', MI_bins, sep=''))
    
    con_g = read_drive_file(con_g_name)
    dis_g = read_drive_file(dis_g_name)
    rownames(con_g) = as.character(con_g$gene) ; con_g = con_g[,-grep('gene', colnames(con_g))]
    rownames(dis_g) = as.character(dis_g$gene) ; dis_g = dis_g[,-grep('gene', colnames(dis_g))]
    
    g_row = intersect(g_set, intersect(rownames(con_g), rownames(dis_g)))
    s_set = intersect(colnames(con_g), colnames(dis_g))
    con_g = con_g[g_row, s_set] ; dis_g = dis_g[g_row, s_set]
    
    Nrow = dim(con_g)[1] ; per_bin = round(Nrow/number_bins)
    t1 = Sys.time()
    #registerDoMC(cores=4)
    out_list <- foreach(sam_ = s_set) %dopar% {
      #for(sam_ in s_set) {
      vv = con_g[sam_] ; vv$bin = 0 ; vv[genes_in_mod,]$bin = 1 ; sign_fact = sign(cor(vv)[1,2])
      vv = dis_g[sam_]
      l_ = lapply(seq(number_bins), function(i) rownames(vv)[vv[sam_][,1] == i])
      c_vec = unlist(lapply(lapply(l_, intersect, genes_in_mod), length))
      nc_vec= as.numeric(sapply(l_, length) - c_vec )
      mi_ = calculateMI_v2(c_vec, nc_vec) ; mi_sign = sign_fact*mi_
      #mi_across_samples_params = c(mi_across_samples_params, mi_sign)
    }
    #tmp_ = data.frame(sample_analyzed = s_set, MI_signed = mi_across_samples_params)
    mi_vec = unlist(out_list) ; abs_mi = abs(mi_vec) ; ix_0 = which(abs_mi < v_rand_q)
    v_z = (abs_mi - mu_r)/sd_r ; v_z[v_z < 0] = 0 ; v_z[is.na(v_z)] = 0 #; v_z[ix_0] = 0
    v_z = v_z*sign(mi_vec)
    tmp_ = data.frame(SAMPLE_ID = s_set, MPS = v_z) ; rownames(tmp_) = s_set
    #saveRDS(tmp_, paste(parent_dir, 'tmp_new_module/', dis_,'_MPS', sep=''))
    t2 = Sys.time() ; print(paste(dis_, (t2-t1)))
    d_mps = tmp_
    samps = as.character(d_mps$SAMPLE_ID)
    #d_clin = readRDS(paste(parent_dir, 'data/', 'disease_clinical/', collection, '/', d_t, '_clinical_primary_forMPS.rds', sep=''))
    
    d_clin_name = paste0(d_t, '_clinical_primary_forMPS.rds')
    d_clin = read_drive_file(d_clin_name)
    
    s_p = samps[which(d_mps$MPS > ((1)*MPS_thresh))] ; s_n = samps[which(d_mps$MPS < ((-1)*MPS_thresh))]
    tmp_df_p = data.frame(SAMPLE_ID = s_p) ; tmp_df_p$group = 'MPS+' ; tmp_df_p$MPS_groups = 'MPS+'
    tmp_df_n = data.frame(SAMPLE_ID = s_n) ; tmp_df_n$group = 'MPS-' ; tmp_df_n$MPS_groups = 'MPS-'
    all_clin = merge(d_clin, rbind(tmp_df_p, tmp_df_n))
    all_clin = all_clin[! is.na(all_clin$group),]
    all_clin$OS.time = all_clin$OS.time/30 ; all_clin$PFI.time = all_clin$PFI.time/30
    #all_clin$path_id = path_id
    all_new_clin = rbind(all_new_clin, all_clin)
  }
  all_new_clin$genes_in_module = paste(genes_in_mod, collapse='|')
  all_new_clin$module_name = module_name ; all_new_clin$module_id = module_id
  #saveRDS(all_new_clin, paste(parent_dir, 'tmp_new_module/all_clin_MPSth_',MPS_thresh, sep=''))
  T_2 = Sys.time() ; print(T_2 - T_1)
  val_ret = all_new_clin
}

split_genes <- function(genes, chunk_size = 15) {
  n <- length(genes)
  split(genes, rep(1:ceiling(n / chunk_size), each = chunk_size, length.out = n))
}

#-- FUNCTION 2: Survival Plot for new module

plot_Surv_get_clinicalParam_newModule <- function(clin_data) {
  pl_list_ovs = pl_list_pfs = list()
  genes_in_module = unlist(strsplit(as.character(unique(clin_data$genes_in_mod)), '\\|'))
  total_patients = dim(clin_data)[1]
  #readRDS(paste(parent_dir, 'tmp_new_module/input_genes', sep=''))
  num_genes_in_module = length(genes_in_module)
  name_module_trim = as.character(unique(clin_data$module_name))
  for(dis_ in as.character(unique(clin_data$cohort))) {
    all_clin = clin_data[which(clin_data$cohort == dis_),]
    pl_o = pl_p = nullGrob()
    if (sum(table(all_clin$group) > 20) == length(unique(all_clin$group))) {
      surv_ovs = getSurv_shiny(all_clin, 'OS',1) ; surv_pfs = getSurv_shiny(all_clin, 'PFI',1)
      fdr_ovs = as.numeric(surv_ovs[[12]]) ; fdr_pfs = as.numeric(surv_pfs[[12]])
      if (fdr_ovs == 0) { f_ovs = '1e-3' } ; if (fdr_ovs != 0) { f_ovs = formatC(fdr_ovs,format='e',digits=1) }
      if (fdr_pfs == 0) { f_pfs = '1e-3' } ; if (fdr_pfs != 0) { f_pfs = formatC(fdr_pfs,format='e',digits=1) }
      
      m_ovs= paste(dis_, 'p=', formatC(as.numeric(surv_ovs[[11]]$pval_surv),format='e',digits=1), 
                   ' | fdr<', f_ovs,' | HR=', signif(as.numeric(surv_ovs[[11]]$HZR_cox), 2))
      m_ovs= paste(dis_, 'p=', formatC(as.numeric(surv_ovs[[11]]$pval_surv),format='e',digits=1), 
                   ' | HR=', signif(as.numeric(surv_ovs[[11]]$HZR_cox), 2))
      m_sub_ovs = paste('OVS | MPS+ ',sum(all_clin$group == "MPS+"),'| MPS- ',sum(all_clin$group == "MPS-")," | n= ",total_patients, " | ", name_module_trim, sep='')
      m_pfs= paste(dis_, 'p=', formatC(as.numeric(surv_pfs[[11]]$pval_surv),format='e',digits=1), 
                   ' | fdr<', f_pfs,' | HR=', signif(as.numeric(surv_pfs[[11]]$HZR_cox), 2))
      m_pfs= paste(dis_, 'p=', formatC(as.numeric(surv_pfs[[11]]$pval_surv),format='e',digits=1), 
                   ' | HR=', signif(as.numeric(surv_pfs[[11]]$HZR_cox), 2))
      m_sub_pfs = paste('PFS | MPS+ ',sum(all_clin$group == "MPS+"),'| MPS- ',sum(all_clin$group == "MPS-"), "| n= ",total_patients," | ", name_module_trim, sep='')
      
      pl_ovs = ggsurvplot(surv_ovs[[2]], surv_ovs[[3]], pval=F, title=m_ovs, font.title=12, censor.shape=124,censor.size=2,
                          subtitle=m_sub_ovs, font.subtitle=10,surv.median.line='hv', palette =col_surv, risk.table=F)
      
      pl_pfs = ggsurvplot(surv_pfs[[2]], surv_pfs[[3]], pval=F, title=m_pfs, font.title=12, censor.shape=124,censor.size=2,
                          subtitle=m_sub_pfs, font.subtitle=10,surv.median.line='hv', palette =col_surv, risk.table=F)
      pl_o = pl_ovs$plot ; pl_p = pl_pfs$plot
    }
    #pl_list_ovs[[dis_]] = pl_o ; pl_list_pfs[[dis_]] = pl_p
    pl_list_ovs = pl_o ; pl_list_pfs = pl_p
  }
  val_ret = list(pl_list_ovs, pl_list_pfs)
  
}


#---- PART 2 : Survival Analysis for Existing Module ----#

#---Under construction for optimization --#



#--- PART 3 : Common plots for Genes --------#


plot_moduleGenes <- function(clin_data) {
  all_clin = clin_data
  dis_ = as.character(unique(all_clin$cohort)) ; module_id = as.character(unique(all_clin$module_id))
  name_module = as.character(unique(all_clin$module_name)) ; name_module_trim = strtrim(name_module, 14)
  genes_in_module = (unlist(strsplit(as.character(unique(all_clin$genes_in_mod)), '\\|')))
  num_genes_in_module = length(genes_in_module)
  gene_info_filt <- get_gene_info_from_drive(gene_info_file_name)
  ez_ids = unique(as.character(merge(gene_info_filt, data.frame(Approved.Symbol = genes_in_module))$Entrez.Gene.ID))
  ez_ids <- unique(as.character(merge(gene_info_filt, data.frame(Approved.Symbol = genes_in_module))$Entrez.Gene.ID))
  
  # Check if ez_ids is empty, return NULL plot and message
  if (length(ez_ids) == 0) {
    return(list(NULL, num_genes_in_module, name_module, NULL, genes_in_module))
  }
  
  # Perform enrichment analysis only if ez_ids are available
  
  #-- Problem with cluster profiler package --#
  enr_go <- enrichGO(gene = ez_ids, organism, pAdjustMethod = "fdr", ont = 'BP',
                     pvalueCutoff = 1e-3, minGSSize = 30, maxGSSize = 5000)
  
  # Check if enrichment analysis returned valid results
  if (is.null(enr_go) || nrow(enr_go@result) == 0) {
    pl_enrich <- NULL
  } else {
    pl_enrich <- enrichplot::dotplot(enr_go, x = 'count', showCategory = 15)
  }
  
  # Return a list with plot, metadata, and enrichment results
  val_ret <- list(pl_enrich, num_genes_in_module, name_module, enr_go, genes_in_module)
  return(val_ret)
}


get_geneExp_info <- function(clin_data) {
  all_clin = clin_data
  dis_ = as.character(unique(all_clin$cohort)) ; path_id = as.character(unique(all_clin$path_id))
  s_p = as.character(all_clin$SAMPLE_ID[which(all_clin$group == 'MPS+')])
  s_n = as.character(all_clin$SAMPLE_ID[which(all_clin$group == 'MPS-')])
  # path_id = 'FIRE1KBup_397' ; dis_ = 'breast' ; zs_thresh = 0
  d_module = readRDS(paste(parent_dir, 'shiny_test/data/module_info/', path_id, '_rds', sep=''))
  
  genes_in_module = intersect(as.character(gene2module$Approved.Symbol), unlist(strsplit(as.character(d_module$genes), '\\|')))
  num_genes_in_module = length(genes_in_module)
  
  d_genev = vroom(paste(parent_dir, 'shiny_test/data/disease_gene_expression/', dis_, '_tumdat_zscore', sep=''), delim = "\t",col_names=T)
  #d_gene_binsv = vroom(paste(parent_dir, 'shiny_test/data/disease_gene_expression/', dis_, '_tumdata_zscore_bins10', sep=''), delim = "\t",col_names=T)
  
  #genes_full = intersect(d_genev$gene, d_gene_binsv$gene)
  d_gene_module = as.data.frame(d_genev[which(d_genev$gene %in% genes_in_module),])
  rownames(d_gene_module) = as.character(d_gene_module$gene)
  d_gene_module = d_gene_module[,-grep('gene',colnames(d_gene_module))]
  genes_in_module = intersect(rownames(d_gene_module), genes_in_module)
  dg_p = d_gene_module[,s_p] ; dg_n = d_gene_module[,s_n]
  pv_vec = sapply(genes_in_module, function(g_) { 
    v_ab = as.numeric(dg_p[g_,]) ; 
    v_be = as.numeric(dg_n[g_,]) ; 
    w_ = wilcox.test(v_ab, v_be) ; w_$p.value 
  })
  pv_vec[is.na(pv_vec)] = 1 ; fdr_vec = p.adjust(pv_vec, 'fdr')
  df_g = data.frame(gene = names(pv_vec), pv = pv_vec, fdr = fdr_vec)
  
  d_genev = readRDS(paste(parent_dir, 'shiny_test/data/disease_gene_expression/', dis_, '_log2dat_rds', sep=''))
  #intersect(as.character(d_genev, genes_in_module))
  d_gene_module = d_genev[intersect(as.character(d_genev$gene), genes_in_module), grep('TCGA',colnames(d_genev))]
  #median_gene_module = as.numeric(apply(d_gene_module[,c(s_p, s_n)], 2, mean))
  
  dg_p = d_gene_module[,s_p] ; dg_n = d_gene_module[,s_n]
  lR_vec = (apply(dg_p, 1, median) - apply(dg_n, 1, median))
  tmp_ = data.frame(gene = names(lR_vec), lR = as.numeric(lR_vec)) ; df_g = merge(df_g, tmp_)
  df_g$path_id = path_id ; df_g$cohort = dis_ ; rm(tmp_)
  eps_th = min(c(1e-16, min(df_g$fdr[df_g$fdr !=0]))) ; df_g$fdr[which(df_g$fdr == 0)] = eps_th
  val_ret = df_g
}

plot_geneExp <- function(gene_DE_data) {
  df_g = gene_DE_data
  pl_volc = nullGrob()
  if (dim(df_g)[1] > 25) {
    pl_volc = ggplot(df_g, aes(lR, -log10(fdr))) + geom_point(col=rgb(0.5, 0.5, 0.5, 0.25)) + theme_classic() +
      geom_hline(yintercept = -log10(0.01),col=rgb(1,0,0,0.2)) + geom_vline(xintercept=0, col=rgb(1,0,0,0.2)) +
      xlab('log2-ratio MPS+ vs. MPS-') + ylab('-log10(fdr)')
    
  }
  val_ret = list(pl_volc)
}


#--- PART 4 : Common plots for Histopathy --------#

get_commonHistopath <- function(disease_name) {
  dis_ = disease_name
  hist_v = hist_v_ori = c('histological_type', 'coarse_pathologic_stage', 'age_category')
  if (dis_ == 'breast') { hist_v = union(c('HR_category', 'HER2_category', 'TN_category'), hist_v_ori) }
  
  if (dis_ == 'prostate') { hist_v = union(c('gleason_score_category','PSA_value'), hist_v_ori) }
  
  if (dis_ == 'AML') { hist_v = union(c('mol_test_status'), hist_v_ori) }
  
  if (dis_ == 'head_neck') { hist_v = union(c('HPV_status'), hist_v_ori) }
  
  if (dis_ == 'cervical') {  hist_v = union(c('HPV_status'), hist_v_ori) }
  
  if (dis_ == 'colon') { hist_v = union(c('MSI_status','tumor_side'), hist_v_ori) }
  
  val_ret = hist_v
}


plot_Histopath <- function(clin_data) {
  all_clin = clin_data
  dis_ = as.character(unique(all_clin$cohort)) ; path_id = as.character(unique(all_clin$module_id))
  
  
  hist_v = get_commonHistopath(dis_) #c('histological_type', 'coarse_pathologic_stage', 'age_groups')
  dis_hist = all_clin[,c(hist_v, 'group')]
  pl_hist = nullGrob()
  if (sum(table(dis_hist$group) > 20) == length(unique(dis_hist$group))) {
    pl_hist = list()
    for(h_v in sort(hist_v)) {
      if (h_v != 'PSA_value') {
        tmp_h = dis_hist[,union('group', h_v)] ; colnames(tmp_h)[which(colnames(tmp_h) == h_v)] = 'var'
        t_ = table(as.character(tmp_h$var)) ; tmp_h$var = as.character(tmp_h$var)
        if (length(t_) > 1) {
          for(i in names(t_)) { ix_ = which(tmp_h$var == i) ; tmp_h$var[ix_] = paste(i, ' (', as.numeric(t_[i]), ')',sep='') }
          pl_ = ggplot(tmp_h, aes(var, fill=group)) + geom_bar(stat = "count",position = "fill", width=0.25) + 
            scale_y_continuous(labels=scales::percent) + theme_bw() + xlab('') + #coord_flip() +
            scale_fill_manual(values=col_surv) + ggtitle(gsub('_', ' ', h_v)) + ylab('percent') + coord_flip()#+ facet_wrap(~ hist_var, scales = "free_x")
          pl_hist[[h_v]] = pl_
        }
      }
      if (h_v == 'PSA_value') {
        tmp_h = dis_hist[,union('group', h_v)] ; colnames(tmp_h)[which(colnames(tmp_h) == h_v)] = 'var'
        pl_ = ggplot(tmp_h, aes(group, log2(var), fill=group)) + theme_bw() + xlab('') + ylab(paste('log2(', gsub('_', ' ', h_v), ') a.u')) + 
          geom_violin(alpha=0.5, draw_quantiles=c(0.25, 0.5, 0.75)) + scale_fill_manual(values=col_surv) + 
          ggtitle(gsub('_', ' ', h_v)) + coord_flip()
        pl_hist[[h_v]] = pl_
      }
    }
  }
  val_ret = pl_hist
  
}

plot_Histopath_multivariate <- function(clin_data) {#1 = 'group', var2 = 'histology') {
  all_clin = clin_data ; rownames(all_clin) = all_clin$sample_analyzed
  #dis_ = as.character(unique(all_clin$disease_type)) ; path_id = as.character(unique(all_clin$path_id))
  dis_ = as.character(unique(all_clin$cohort)) ; path_id = as.character(unique(all_clin$module_id))
  var_vec = c('group', get_commonHistopath(dis_)) #c('histological_type', 'coarse_pathologic_stage', 'age_groups')
  ix_rm = c() ; var_vec_ex = c()
  for(v_v in var_vec) {
    all_clin[v_v][,1] = as.character(all_clin[v_v][,1])
    ix_na = which(is.na(all_clin[v_v][,1]))
    if (dim(all_clin)[1] - length(ix_na) < 30) { var_vec_ex = c(var_vec_ex, v_v) }
    if (dim(all_clin)[1] - length(ix_na) >= 30) {
      ix_rm = union(ix_rm, which(is.na(all_clin[v_v][,1]))) #; print(ix_rm) ; print(v_v)
    }
    t_ = table(as.character(all_clin[v_v][,1])) ; cat_n = names(which(t_ <= 10))
    if ((length(t_) - length(cat_n)) > 1 & length(cat_n) >= 1) {
      for(n_ in cat_n) { ix_rm = union(ix_rm, which(all_clin[v_v][,1] == n_)) }
    }
    if ((length(t_) - length(cat_n)) <= 1) { var_vec_ex = c(var_vec_ex, v_v) }
  }
  var_vec = setdiff(var_vec, var_vec_ex) ; all_clin = all_clin[setdiff(seq(dim(all_clin)[1]), ix_rm), ]
  #all_clin$mod_censor_time = all_clin$mod_censor_time/30
  all_clin_num = all_clin ; for(v_v in var_vec) { all_clin_num[v_v][,1] = as.numeric(factor(all_clin_num[v_v][,1])) }
  ss_ovs = Surv(as.numeric(as.character(all_clin$OS.time)),all_clin$OS)
  ss_pfs = Surv(as.numeric(as.character(all_clin$PFI.time)),all_clin$PFI)
  val_ret = list(NA, NA, NA) ; p_for_ovs = p_for_pfs = nullGrob()
  if (! (sum(is.na(ss_ovs)) == length(is.na(ss_ovs)))) {
    su_ovs = Surv((all_clin$OS.time), (all_clin$OS))
    s_cox_ovs = coxph(as.formula(paste('su_ovs ~', paste(var_vec, collapse='+'))), data=all_clin_num)
    p_for_ovs = ggforest(s_cox_ovs,all_clin_num, main='OS',refLabel = "reference")
    
  }
  if (! (sum(is.na(ss_pfs)) == length(is.na(ss_pfs)))) {
    su_pfs = Surv((all_clin$PFI.time), (all_clin$PFI))
    s_cox_pfs = coxph(as.formula(paste('su_pfs ~', paste(var_vec, collapse='+'))), data=all_clin_num)
    p_for_pfs = ggforest(s_cox_pfs,all_clin_num, main='PFI',refLabel = "reference")
    if (dis_ == 'AML') { p_for_pfs = nullGrob() }
  }
  val_ret = list(p_for_ovs, p_for_pfs)
  
}


plot_Histopath_subset <- function(clin_data) {#1 = 'group', var2 = 'histology') {
  all_clin = clin_data ; rownames(all_clin) = all_clin$SAMPLE_ID
  dis_ = as.character(unique(all_clin$cohort)) ; path_id = as.character(unique(all_clin$module_id))
  var_vec = c('group', get_commonHistopath(dis_))
  total_patients = dim(clin_data)[1]
  MPSp = sum(clin_data$group == "MPS+")
  MPSn = sum(clin_data$group == "MPS-")#c('histological_type', 'coarse_pathologic_stage', 'age_groups')
  
  ix_rm = c() ; var_vec_ex = c()
  for(v_v in var_vec) {
    all_clin[v_v][,1] = as.character(all_clin[v_v][,1])
    ix_na = which(is.na(all_clin[v_v][,1]))
    if (dim(all_clin)[1] - length(ix_na) < 30) { var_vec_ex = c(var_vec_ex, v_v) }
    if (dim(all_clin)[1] - length(ix_na) >= 30) {
      ix_rm = union(ix_rm, which(is.na(all_clin[v_v][,1]))) #; print(ix_rm) ; print(v_v)
    }
    t_ = table(as.character(all_clin[v_v][,1])) ; cat_n = names(which(t_ <= 10))
    if ((length(t_) - length(cat_n)) > 1 & length(cat_n) >= 1) {
      for(n_ in cat_n) { ix_rm = union(ix_rm, which(all_clin[v_v][,1] == n_)) }
    }
    if ((length(t_) - length(cat_n)) <= 1) { var_vec_ex = c(var_vec_ex, v_v) }
  }
  var_vec = setdiff(var_vec, var_vec_ex) ; all_clin = all_clin[setdiff(seq(dim(all_clin)[1]), ix_rm), ]
  c_ = 1 ; pl_list_ovs = pl_list_pfs = list()
  for(v_v in setdiff(var_vec, 'group')) {
    var_cat = as.character(unique(all_clin[v_v][,1]))
    for(v_c in var_cat) {
      all_clin_var = all_clin[which(all_clin[v_v][,1] == v_c),]
      pl_o = pl_p = nullGrob()
      if (sum(table(all_clin_var$group) > 5) == length(unique(all_clin_var$group))) {
        surv_ovs = getSurv_shiny(all_clin_var, 'OS',1) ; surv_pfs = getSurv_shiny(all_clin_var, 'PFI',1)
        m_ovs= paste(v_v,  ' | p=', formatC(as.numeric(surv_ovs[[11]]$pval_surv),format='e',digits=1), 
                     ' | HR=', signif(as.numeric(surv_ovs[[11]]$HZR_cox), 2))
        m_sub_ovs = paste('OVS (', v_c,')', "| MPS+ ", MPSp, "| MPS- ", MPSn, "| n= ", total_patients, sep='')
        m_pfs= paste(v_v, ' | p=', formatC(as.numeric(surv_pfs[[11]]$pval_surv),format='e',digits=1), 
                     ' | HR=', signif(as.numeric(surv_pfs[[11]]$HZR_cox), 2))
        m_sub_pfs = paste('PFS (', v_c,')', "| MPS+ ", MPSp, "| MPS- ", MPSn, "| n= ", total_patients, sep='')
        pl_ovs = ggsurvplot(surv_ovs[[2]], surv_ovs[[3]], pval=F, title=m_ovs, font.title=12, censor.shape=124,censor.size=2,
                            subtitle=m_sub_ovs, font.subtitle=10,surv.median.line='hv', palette =col_surv, risk.table=F)
        pl_pfs = ggsurvplot(surv_pfs[[2]], surv_pfs[[3]], pval=F, title=m_pfs, font.title=12, censor.shape=124,censor.size=2,
                            subtitle=m_sub_pfs, font.subtitle=10,surv.median.line='hv', palette =col_surv, risk.table=F)
        pl_o = pl_ovs$plot ; pl_p = pl_pfs$plot ; pl_list_ovs[[c_]] = pl_o ; pl_list_pfs[[c_]] = pl_p
        if (dis_ == 'AML') { pl_list_pfs[[c_]] = nullGrob() } ; c_ = c_ + 1
      }
    }
  }
  val_ret = list(pl_list_ovs, pl_list_pfs) 
  
}



#-- Generate the html --#

generatehtml <- function(file, data) {
  tmp_geneEnrichPlot <- tempfile(fileext = ".png")
  tmp_survPlot <- tempfile(fileext = ".png")
  tmp_histPlot <- tempfile(fileext = ".png")
  tmp_histPlot2 <- tempfile(fileext = ".png")
  tmp_histPlot3 <- tempfile(fileext = ".png")
  tmp_histPlot4 <- tempfile(fileext = ".png")
  
  pl <- data
  
  if (!is.null(pl$geneEnrichPlot)) {
    ggsave(tmp_geneEnrichPlot, plot = pl$geneEnrichPlot, width = 7, height = 7)
  } else {
    tmp_geneEnrichPlot <- NULL
    warning("Gene enrichment plot object is NULL")
  }
  
  if (!is.null(pl$survPlot)) {
    ggsave(tmp_survPlot, plot = grid.arrange(grobs = pl$survPlot, ncol = 2), width = 7, height = 7)
  } else {
    warning("Survival plot object is NULL")
  }
  
  if (!is.null(pl$histPlot)) {
    ggsave(tmp_histPlot, plot = grid.arrange(grobs = pl$histPlot, ncol = ceiling(sqrt(length(pl$histPlot)))), width = 7, height = 7)
  } else {
    warning("Histopathology plot object is NULL")
  }
  
  if (!is.null(pl$histPlot2) && length(pl$histPlot2) > 0) {
    n_col <- ceiling(sqrt(length(pl$histPlot2)))
    ggsave(tmp_histPlot2, plot = grid.arrange(grobs = pl$histPlot2, ncol = n_col), height = 8, width = 9)
  } else {
    warning("Histopathology plot 2 object is NULL or empty")
  }
  
  if (!is.null(pl$histPlot3) && length(pl$histPlot3) > 0) {
    n_col <- ceiling(sqrt(length(pl$histPlot3)))
    ggsave(tmp_histPlot3, plot = grid.arrange(grobs = pl$histPlot3, ncol = n_col), height = 8, width = 9)
  } else {
    warning("Histopathology plot 3 object is NULL or empty")
  }
  
  if (!is.null(pl$histPlot4) && length(pl$histPlot4) > 0) {
    n_col <- ceiling(sqrt(length(pl$histPlot4)))
    ggsave(tmp_histPlot4, plot = grid.arrange(grobs = pl$histPlot4, ncol = n_col), height = 8, width = 9)
  } else {
    warning("Histopathology plot 4 object is NULL or empty")
  }
  
  rmarkdown::render("/lustre/fs4/home/schhabria/pipelines/R-Shiny-SC/report.Rmd", output_file = file,
                    params = list(
                      #data = data(),
                      gene_enrich_plot = tmp_geneEnrichPlot,
                      surv_plot = tmp_survPlot,
                      hist_plot = tmp_histPlot,
                      hist_plot2 = tmp_histPlot2,
                      hist_plot3 = tmp_histPlot3,
                      hist_plot4 = tmp_histPlot4
                      #moduleInfo = pl[[3]],
                      #geneCount = pl[[2]],
                      #genes = pl[[5]]
                    ))
}

#---- Calling the app ----#

ui <- fluidPage(
  title = 'icamp (user)',
  titlePanel("Inferring Clinical Associations of Module Perturbations (user)"),
  sidebarLayout(
    sidebarPanel(
      selectInput("disease_name", "Select Patient Cohort",
                  choices = c("pancreatic", "AML", "bladder", "breast", "cervical", "colon", "esophageal", "GBM",
                              "glioma", "head_neck", "kidney", "liver", "lung_ad", "lung_sq",
                              "melanoma", "ovarian", "pancreatic", "paraganglioma", "prostate",
                              "rectal", "sarcoma", "stomach", "testicular", "thymoma", "thyroid","uterine_endometrial")
      ),
      fileInput("module_inp", "Choose CSV File",
                accept = c( "text/csv",
                            "text/comma-separated-values,text/plain",
                            ".csv")
      ),
      actionButton("submit", "Submit")
      #downloadButton("downloadhtml", "Download html Report")
    ),
    mainPanel(
      uiOutput("message"),
      tabsetPanel(type = "tabs",
                  tabPanel("Dashboard", htmlOutput("readme")),
                  tabPanel("Module Info (Genes)", htmlOutput('module_info'),
                           #plotOutput(outputId = 'geneEnrichPlot', click = 'plot_click'),
                           verbatimTextOutput('volc_info')
                           #downloadButton("download_geneEnrichPlot", "Download Gene Enrichment Plot")
                  ),
                  tabPanel("Clinical Data", htmlOutput('table_info'),                           
                           htmlOutput('clinical_info'),
                           downloadButton("downloadData", "Download clinical information of the patients")
                  ), 
                  
                  tabPanel(
                    "Survival",
                    div(
                      class = "plot-container",
                      downloadButton("download_survPlot", "Download Survival Plot"),
                      plotOutput("survPlot", height = "400px")  
                    ),
                    tags$style(
                      type = "text/css",
                      ".plot-container { position: relative; width = 100%; height: 400px; text-align: center; } 
                                  #download_survPlot { margin-button: 10px; }"
                    )
                  ),
                  tabPanel(
                    "histopathology",
                    div(
                      class = "plot-container",
                      downloadButton("download_histPlot", "Download Histopathology plot"),
                      plotOutput("histPlot", height = "400px")
                    ),
                    tags$style(
                      type = "text/css",
                      ".plot-container { position: relative; width: 100%; text-align: center; } 
                                    #download_histPlot { margin-bottom: 10px; }"
                    )
                  ),
                  
                  
                  tabPanel(
                    "histopathology (multivariate)",
                    div(
                      class = "plot-container",
                      downloadButton("download_histPlot2", "Download histopathology (multivariate) plot"),
                      plotOutput("histPlot2", height = "400px")
                    ),
                    tags$style(
                      type = "text/css",
                      ".plot-container { position: relative; width: 100%; text-align: center; } 
                                    #download_histPlot2 { margin-bottom: 10px; }"
                    )
                  ),
                  tabPanel(
                    "histopathology (subset: OVS)",
                    div(
                      class = "plot-container",
                      downloadButton("download_histPlot3", "Download histopathology (subset: OVS) plot"),
                      plotOutput("histPlot3", height = "400px")
                    ),
                    tags$style(
                      type = "text/css",
                      ".plot-container { position: relative; width: 100%; text-align: center; } 
                                    #download_histPlot3 { margin-bottom: 10px; }"
                    )
                  ),
                  
                  tabPanel(
                    "histopathology (subset: PFS)",
                    div(
                      class = "plot-container",
                      downloadButton("download_histPlot4", "Download histopathology (subset: OVS) plot"),
                      plotOutput("histPlot4", height = "400px")
                    ),
                    tags$style(
                      type = "text/css",
                      ".plot-container { position: relative; width: 100%; text-align: center; } 
                                    #download_histPlot4 { margin-bottom: 10px; }"
                    )
                  )
                  #tabPanel("histopathology", column(12, plotOutput(outputId='histPlot')),
                  #downloadButton("download_histPlot", "Download Histopathology Plot")),
                  
                  #tabPanel("histopathology (multivariate)", column(12, plotOutput(outputId='histPlot2')),
                  #downloadButton("download_histPlot2", "Download Histopathology Plot 2"),
                  #tags$style(type="text/css", "#download_histPlot2 { position: absolute; bottom: 10px; left: 50%; transform: translateX(-50%); }")),
                  
                  # tabPanel("histopathology (subset: OVS)", 
                  #     column(12, plotOutput(outputId='histPlot3')),
                  #     tags$style(type="text/css", "#download_histPlot3 { position: absolute; bottom: 10px; left: 50%; transform: translateX(-50%); }"),
                  #     downloadButton("download_histPlot3", "Download Histopathology OVS Plot")
                  #     ),
                  
                  # tabPanel("histopathology (subset: PFS)", column(12, plotOutput(outputId='histPlot4')),
                  #     downloadButton("download_histPlot4", "Download Histopathology PFS Plot"))
      )
      
    )
    
  )
)



server <- function(input, output,session) {
  
  # Reactive expression to store plot data
  plot_data <- reactiveVal(NULL)
  data <- reactiveVal(NULL)
  error_message <- reactiveVal(NULL)
  
  observeEvent(input$submit, {
    # Show modal dialog when submit button is clicked
    showModal(modalDialog("Processing data. Please wait...", footer = NULL))
    
    req(input$module_inp)
    inFile <- input$module_inp
    input_genes <- read.csv(inFile$datapath, sep = ',')
    disease_name <- isolate(input$disease_name)
    clin_data_val <- isolate(get_MPS_newModule(input_genes, select_cohorts = disease_name))
    
    if (!is.null(clin_data_val)) {
      # Generate plots based on clinical data
      pl <- plot_moduleGenes(clin_data_val)
      survPlot <- plot_Surv_get_clinicalParam_newModule(clin_data_val)
      histPlot <- plot_Histopath(clin_data_val)
      histPlot2 <- plot_Histopath_multivariate(clin_data_val)
      histPlot3 <- plot_Histopath_subset(clin_data_val)
      
      # Store plot data in reactive value
      plot_data_list <- list(
        geneEnrichPlot = pl[[1]],
        moduleInfo = pl[[3]],
        geneCount = pl[[2]],
        genes = pl[[5]],
        survPlot = survPlot,
        histPlot = histPlot,
        histPlot2 = histPlot2,
        histPlot3 = histPlot3[[1]],
        histPlot4 = histPlot3[[2]]
      )
      
      plot_data(plot_data_list)
      data(clin_data_val)
      removeModal()
      
    } else {
      error_message("An error occurred while processing the new module.")
      showModal(modalDialog(
        title = "Error",
        div(style = "color: red; font-weight: bold;", error_message()),
        easyClose = FALSE,
        footer = modalButton("Close")
      ))
    }
})
  
  output$error_message <- renderText({
    error_message()
  })
  
  output$readme <- renderUI({
    readme_text <- "
    <h2><strong>Welcome to Your Shiny Dashboard!</strong></h2>
    
    <p>This tool is crafted to facilitate your exploration of the expression patterns of biologically coherent genes defined as modules within the TCGA cohort. 
    It analyzes the coordinated shifts in mRNA expression among the genes listed in the module and conducts survival analysis based on these findings. 
    Here’s a quick guide to help you get started::</p>
    
    <h3 style=\"color: blue;\"> ** Terms used in this dashboard ** </h3>
    <ul>
        <li><strong>Modules:</strong> List of biologically coherant genes like pathway.</li>
        <li><strong>Module Pertubation Score (MPS):</strong> Coordinated shift in the mRNA expression of a set of genes in a patient’s tumor transcriptome</li>
        <ul>
                    <li> MPS+ : up-regulation of genes in the module in the patient</li>
                    <li> MPS- : down-regulation of genes in the module in the patient</li>
        </ul>
    </ul>

    <h3 style=\"color: blue;\"> ** Module Definitions ** </h3>
    <ul>
        <li><h4><strong>New module:</strong></h4>
           Upload the CSV file containing a list of your genes of interest.
        </li>
    </ul>
    
    <p><strong>Click on SUBMIT and switch to the other tabs. We hope you find this dashboard useful and intuitive. Happy exploring!</strong></p>
    
    "
    HTML(readme_text)
  })
  
  #--- Module Data ---#
  output$module_info <- renderText({ 
    pl = plot_data()
    if (is.null(pl)) return()
    
    chunked_genes <- split_genes(pl$genes)
    genes_str <- paste(sapply(chunked_genes, paste, collapse = ", "), collapse = "<br>")
    
    paste('<b> Name of module: </b>', pl$moduleInfo, "<br>",
          '<b> Number of genes in the module: </b>', pl$geneCount, "<br>",         
          '<b> Genes in the module: </b><br>', genes_str, "<br>")
  })
  

  
  
  #--- geneEnrichPlot ---#
  output$geneEnrichPlot <- renderPlot({ 
    pl <- plot_data()
    if (is.null(pl) || is.null(pl$geneEnrichPlot)) {
      plot.new()
      text(0.5, 0.5, "No significance was found.")
      return()
    }  
    plot(pl$geneEnrichPlot)
  }, height = 600, width = 900)
  
  
  output$download_geneEnrichPlot <- downloadHandler(
    filename = function() {
      "gene_enrichment_plot.pdf"
    },
    content = function(file) {
      pl <- plot_data()
      if (!is.null(pl$geneEnrichPlot)) {
        ggsave(file, pl$geneEnrichPlot)
      }
    }
  )
  #--- Clinical data ---#
  
  output$table <- renderTable({
    data()
  })
  
  output$downloadData <- downloadHandler(
    filename = function() {
      "clinical_info.csv"
    },
    content = function(file) {
      x <- data()
      write.csv(x, file, row.names = FALSE)
    })
  
  output$clinical_info <- renderText({
    x = data()
    x = paste('<b> number of patients: ', dim(x)[1], "<br>", 'MPS+ ', sum(x$group == 'MPS+'), ' | MPS- ', sum(x$group == 'MPS-'),"<br>","<br>", "</b>")
    
  })
  
  #-- Survival Plot ---#
  
  output$survPlot <- renderPlot({
    pl <- plot_data()
    if (is.null(pl)) return()
    
    grid.arrange(pl$survPlot[[1]], pl$survPlot[[2]], ncol=2)
  }, height = 600, width = 600)
  
  output$download_survPlot <- downloadHandler(
    filename = function() {
      "survival_plot.pdf"
    },
    content = function(file) {
      pl <- plot_data()
      if (!is.null(pl$survPlot)) {
        ggsave(file, plot = grid.arrange(pl$survPlot[[1]], pl$survPlot[[2]], ncol = 2), width = 7, height = 7)
      }
    }
  )
  
  #--- Histopathy Plot for general category ---#
  output$histPlot <- renderPlot({ 
    pl <- plot_data()
    if (is.null(pl)) return()
    
    hist_plots <- pl$histPlot
    n_col <- ceiling(sqrt(length(hist_plots)))
    grid.arrange(grobs = hist_plots, ncol = n_col)
  }, height = 800, width = 900)
  
  output$download_histPlot <- downloadHandler(
    filename = function() {
      "histopathology_plot.pdf"
    },
    content = function(file) {
      pl <- plot_data()
      if (!is.null(pl$histPlot)) {
        ggsave(file, plot = grid.arrange(grobs = pl$histPlot, ncol = ceiling(sqrt(length(pl$histPlot)))),  width = 9, height = 8)
      }
    }
  )
  #--- Histopathy Plot for Mutli category data ---#
  output$histPlot2 <- renderPlot({ 
    pl <- plot_data()
    p_ <- pl$histPlot2
    n_col <- ceiling(sqrt(length(p_)))
    if (n_col == 0) { 
      grid.arrange(nullGrob()) 
    } else {
      grid.arrange(grobs = p_, ncol = n_col) 
    }
  }, height = 800, width = 900)
  
  # Download histPlot2
  output$download_histPlot2 <- downloadHandler(
    filename = function() {
      "histopathology_plot2.pdf"
    },
    content = function(file) {
      pl <- plot_data()
      if (!is.null(pl$histPlot2)) {
        ggsave(file, plot = grid.arrange(grobs = pl$histPlot2, ncol = 2),  width = 8, height = 9)
      }
    }
  )
  
  #--- Histopathy Plot for OVS ---#
  
  output$histPlot3 <- renderPlot({ 
    pl <- plot_data()
    if (is.null(pl)) return()
    
    p_ <- pl$histPlot3
    n_col <- ceiling(sqrt(length(p_)))
    if (n_col == 0) { 
      grid.arrange(nullGrob()) 
    } else {
      grid.arrange(grobs = p_, ncol = n_col) 
    }
  }, height = 800, width = 900)
  
  output$download_histPlot3 <- downloadHandler(
    filename = function() {
      "histopathology_plot_OVS.pdf"
    },
    content = function(file) {
      pl <- plot_data()
      if (!is.null(pl$histPlot3)) {
        ggsave(file, plot = grid.arrange(grobs = pl$histPlot3, ncol = ceiling(sqrt(length(pl$histPlot3)))), width = 9, height = 8)
      }
    }
  )
  
  #--- Histopathy Plot for PFS ---#
  
  output$histPlot4 <- renderPlot({ 
    pl <- plot_data()
    if (is.null(pl)) return()
    
    p_ <- pl$histPlot4
    n_col <- ceiling(sqrt(length(p_)))
    if (n_col == 0) { 
      grid.arrange(nullGrob()) 
    } else {
      grid.arrange(grobs = p_, ncol = n_col) 
    }
  }, height = 800, width = 900)
  
  
  output$download_histPlot4 <- downloadHandler(
    filename = function() {
      "histopathology_plot_PFS.pdf"
    },
    content = function(file) {
      pl <- plot_data()
      if (!is.null(pl$histPlot4)) {
        ggsave(file, plot = grid.arrange(grobs = pl$histPlot4, ncol = ceiling(sqrt(length(pl$histPlot4)))), width = 9, height = 8)
      }
    }
  )
  

}


shinyApp(ui = ui, server = server)
