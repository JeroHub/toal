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
    }else if(length(grep(pattern = 'Start Processing',x = line,fixed = T)) > 0){
      # When start prcessing is read, start reading table values
      readtable <- T
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
  return(data.raw)
}
