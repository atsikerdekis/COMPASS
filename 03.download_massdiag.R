args <- commandArgs(trailingOnly=TRUE)

if (length(args) != 4) {
  stop("Usage: Rscript 03.download_massdiag.R EXPNAME DATESTART DATEEND PATH_DATA")
}

expname   <- args[1]
startDate <- as.Date(args[2],format="%Y%m%d")
endDate   <- as.Date(args[3],format="%Y%m%d")
path_data <- args[4]

outdir <- paste0(path_data,expname,"/massdiag/")
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)

dates <- format(seq.Date(startDate,endDate,by="day"),"%Y%m%d")
wanted <- paste0("massdia_chem__",expname,"_",dates,"00.txt")

existing <- file.exists(paste0(outdir,wanted))

if (all(existing)) {
  message("All requested MASSDIA files already exist.")
  quit(save="no",status=0)
}

tempdir <- tempfile(paste0("massdiag_",expname,"_"))
dir.create(tempdir)

on.exit(unlink(tempdir,recursive=TRUE,force=TRUE),add=TRUE)

command <- paste0("ecp 'ec:/xnl/public/",expname,"/log/m*.tar' ",tempdir,"/")
message(command)

status <- system(command)

if (status != 0) stop("ecp failed for experiment ",expname)

tarfiles <- list.files(tempdir,pattern="\\.tar$",full.names=TRUE)

if (length(tarfiles) == 0) stop("No MASSDIA tar files downloaded for ",expname)

for (tarfile in tarfiles) {
  message("Reading: ",tarfile)

  contents <- untar(tarfile,list=TRUE)

  for (filename in wanted) {
    matches <- contents[basename(contents) == filename]

    if (length(matches) == 0) next

    untar(tarfile,files=matches[1],exdir=tempdir)

    extracted <- file.path(tempdir,matches[1])

    if (file.exists(extracted)) {
      file.copy(extracted,paste0(outdir,filename),overwrite=TRUE)
      message("Extracted: ",paste0(outdir,filename))
    }
  }
}

missing <- wanted[!file.exists(paste0(outdir,wanted))]

if (length(missing) > 0) {
  warning("Missing MASSDIA files: ",paste(missing,collapse=", "))
}
