#' Import raw data for HTI 291 Acoustic Tag Reviever
#'
#' @param file The path to .RAT file holding raw data from reciever.
#' @param fs The sample rate of the file.
#'
#' @return A dataframe containing all detected pulses from all hydrophones.
#' @export
#'
#' @examples
#' file <- 'HTI_Example_Data/AT16S012840914 (G8 4,0,48).RAT'
#' data.raw <- read.HTI.RAT(file)
#'
read.HTI.RAT <- function(file,fs){
  require(R.utils)
  ## Open file for reading
  conn <- file(file,open = 'r')
  total.lines <- countLines(file)

  # Initialize data.frames
  buffer <- 10000
  data.buffer <- data.frame(Sample = integer(length = buffer),
                            Hydrophone = integer(length = buffer),
                            Channel = integer(length = buffer),
                            PulseWidth_3dB = integer(length = buffer),
                            PulseWidth_6dB = integer(length = buffer),
                            PulseWidth_12dB = integer(length = buffer),
                            Peak_Amplitude = integer(length = buffer),
                            Noise_Level = integer(length = buffer),
                            Auto_Threshold = integer(length = buffer))
  data.raw <- data.buffer
  data.raw.size <- nrow(data.raw)

  i <- 1 # Initialize file line number
  j <- 1 # Initialize row index number
  readtable <- F

  cat('Reading HTI file...\n')
  pb <- txtProgressBar(min = 0, max = total.lines,initial = 1,width = 20, style = 3)

  while(i<= total.lines){

    # Read line
    line <- readLines(con = conn,1)

    ## Process data
    if(readtable == T){
      line.vector <- strsplit(line,split = ',')[[1]]
      ## If line looks good, write to data table and increment row counter
      if (length(line.vector) == 9){
        data.raw[j,] <- as.integer(line.vector)
        j <- j + 1
      }
    }else if(length(grep(pattern = 'Start',x = line,fixed = T)) > 0){
      # When start prcessing is read, start reading table values
      readtable <- T
      # Get date and time on which the logging started
      StartDate <- paste(substring(line, regexpr("at", line) + 11, regexpr("at", line) + 12),
                         substring(line, regexpr("at", line) + 7, regexpr("at", line) + 9),
                         substring(line, regexpr("at", line) + 23, regexpr("at", line) + 26), sep = "/")

      lct <- Sys.getlocale("LC_TIME"); Sys.setlocale("LC_TIME", "C")
      StartDate <- as.Date(StartDate, "%d/%b/%Y")
      Sys.setlocale("LC_TIME", lct)

      StartTime <- substring(line, regexpr("at", line) + 14, regexpr("at", line) + 21)
      StartDateTime <- as.POSIXct(paste(StartDate, StartTime, sep = " "))
    }

    # Increment line counter
    i <- i + 1

    # Expand data table if table end is reached
    if(j >= data.raw.size){
      data.raw <- rbind(data.raw,data.buffer)
      data.raw.size <- nrow(data.raw)
    }
    ## Update status bar every 1000 records
    if(i %% 1000 == 0){
      setTxtProgressBar(pb,value = i)
    }
  }
  # Close file connection
  close(pb)
  close(conn)

  ## Trim final datasize (remove empty buffer space at end of table)
  data.raw <- data.raw[1:j,]
  data.raw$seconds <- data.raw$Sample/fs
  data.raw <- data.raw[!(data.raw$Sample==0 & data.raw$Peak_Amplitude==0),]
  # Ascribing a date and time to each detection
  StartSeconds <- min(data.raw$seconds)
  attr(x = data.raw, which = 'StartTime') <- as.character(StartDateTime)

  return(data.raw)
}


#' Scan HTI metadata in folder
#'
#' This function will scan all PAT files in a given folder and return the metadata accociated with them.
#' This can be used for filtering your files by recording time before analysing.
#'
#' @param folder
#'
#' @return data.frame of metadata
#' @export
#'
#' @examples
HTI.meta <- function(folder){
  require(R.utils)

  #folder <- '/mnt/Backup/DaysAgo.0/Jacobahaven_2017/Experiment_2017_FloatingPen/HTI'
  files <- dir(path = folder, pattern = '.PAT')

  ## Match Locale with data set (Englsih, US)
  Sys.setlocale(category = 'LC_TIME','en_US.UTF-8')
  #Sys.getlocale()

  #i <- files[1]
  meta <- data.frame()
  for(i in files){
    full.path <- paste0(folder,'/',i)
    conn <- file(full.path, open = 'r')

    ## Read all lines in file
    total.lines <- countLines(full.path)
    lines <- readLines(con = conn, n = total.lines, ok = T, encoding = 'UTF-8')
    close(conn)

    ## Scan for metadata
    #### Grab start time
    idx.start <- grep(x = lines,pattern = 'Start')
    if(length(idx.start > 0)){
        if(gregexpr('Processing', lines[idx.start])[[1]] > -1){
          start.type <- 'Processing'}else{start.type <- 'Sequence'}
        idx.cut <- c(13+nchar(start.type),
                     gregexpr(pattern = ' Peak Location',text = lines[idx.start])[[1]][1] - 1)
        start.string <- substr(x = lines[idx.start], start = idx.cut[1], stop = idx.cut[2])
        start.posix <- as.POSIXct(strptime(x = start.string,
                                           format = "%a %b %d %H:%M:%S %Y",
                                           tz = ""))
    }else{
      start.posix <- NA
      start.type <- 'NA'
    }

    #### Grab end time
    idx.stop <- grep(x = lines,pattern = 'Stop')
    if(length(idx.stop > 0)){
      if(gregexpr('Processing', lines[idx.stop])[[1]] > -1){
        stop.type <- 'Processing'
      }else{
        stop.type <- 'Sequence'
      }
      idx.cut <- c(12+nchar(stop.type),
                   gregexpr(pattern = ' Peak Location',text = lines[idx.stop])[[1]][1] - 1)
      stop.string <- substr(x = lines[idx.stop], start = idx.cut[1], stop = idx.cut[2])
      stop.posix <- as.POSIXct(strptime(x = stop.string,
                                         format = "%a %b %d %H:%M:%S %Y",
                                         tz = ""))
    }else{
      stop.posix <- NA
      stop.type <- 'NA'
    }

    meta <- rbind(meta,
                  data.frame(file = substr(i,start = 1, stop = nchar(i)-4),
                             start = start.posix,
                             start.type = start.type,
                             stop = stop.posix,
                             stop.type = stop.type))
  }

  return(meta)
}
