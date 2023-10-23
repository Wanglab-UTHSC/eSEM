#' @export motifEx
#' @import stringr
#' @description This script is used to extract kinase-substrate relationships from motif scanning and build motif_input file for KSEM().
#' @param mastFile The output file "extract_result.txt" from Step: 03step_ExtractMASTresult.sh.
#' @param savePath Path to save motif_input.txt file. Default is current directory.
#' @importFrom utils write.csv
#' @title motifEx



motifEx <- function(mastFile,
                    savePath=getwd()){
  mastFile_fix <- str_split_fixed(mastFile$V1, "_MASTout.", 2) %>% data.frame()
  colnames(mastFile_fix) <- c("GENE", "sites")

  file_path <- file.path(savePath, "motif_input.csv")

  write.csv(mastFile_fix, file = file_path, row.names = F)
  return(mastFile_fix)
}


# mastFile <- read.table("C:\\Project_Dehui_2022\\KSEM\\02_KSEM_code\\KSEM_v2\\KSEM\\data_raw\\extract_result.txt")
# motif_example <- motifEx(mastFile = mastFile, savePath = ".\\data_raw")


