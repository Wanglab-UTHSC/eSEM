#' @export getFa
#' @export subKinase
#' @export getSeq
#' @title Build MEME input
#' @description script to extract MEME input .fa file
#' @param ks Vector; Kinase species; ("mouse", "rat", "human")
#' @param filePath The path is used to save .fa files. Recommend to create a new folder to save all .fa.files.
#' @examples getFa(ks = "human", filePath = getwd())
#'
#'
#'
#'
#'
#'



getFa <- function(ks = NULL,
                  filePath){
  if (is.null(ks)) {
    stop("Please set your kinase organism!")
  }else{
    species_kinase <- getSeq(kins = ks)
    kinase_list <- levels(factor(species_kinase$GENE))
    lapply(kinase_list, subKinase, psp_table = species_kinase, savePath = filePath)
  }
}




getSeq <- function(kins){
  # Load psp
  psp <- psp_data()
  # Toupper gene name
  psp$GENE <- toupper(psp$GENE)
  # Replace "-" with "." in SUB_ACC_ID
  psp$SUB_ACC_ID <- as.matrix(gsub("-",".",as.matrix(psp$SUB_ACC_ID)))
  # Extract ORGANISM psp
  psp_species <- psp[which(psp$organism %in% kins), c(1, 7, 10, 12)]

  # Replace "_" with "" in peptide (remove empty amino acid sites)
  psp_species$SITE_...7_AA <- as.matrix(gsub("_", "",as.matrix(psp_species$SITE_...7_AA)))

  # Modify PSP_human
  psp_species$sites <- paste(psp_species$SUB_ACC_ID, psp_species$SUB_MOD_RSD, sep = "_")
  psp_species$nchar <- nchar(psp_species$SITE_...7_AA)

  psp_species$site_name <- paste(">", psp_species$sites, sep = "")
  psp_species$seq <- paste(psp_species$site_name, "\n", psp_species$SITE_...7_AA)
  return(psp_species)
}


subKinase <- function(kinase, psp_table, savePath){
  kinase_table <- psp_table[which(psp_table$GENE == kinase), ]
  kinase_file <- paste(kinase_table$seq)
  kinase_file <- gsub(" ","", kinase_file)
  file_path <- file.path(savePath, paste(kinase, ".fa", sep = ""))
  #file_path <- paste(savePath, "\\", kinase, ".fa", sep = "")
  write.table(kinase_file, file = file_path,
              row.names = F, col.names = F, quote = F)
}


##
