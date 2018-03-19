#' Detect pulse onset frequencies of acoustic tags
#'
#' @param detections A vector of detection times (seconds) from a single hydrophone.
#' @param frequencies A vector of frequencies to test.
#' @param n.tags The number of tags with unique frequencies.
#' @param miNGap The minimum gap between unique tag frequencies (Hz).
#' @param plot Set to T to plot results with ggplot.
#'
#' @details This function uses convolution to find the unique tag frequencies in setups where each tagged individual has a unique pulse interval on their acoustic tag.
#' You can set \code{plot = T} for error checking.
#'
#' @return
#' @export
#'
TagFreq.detect <- function(detections, frequencies, n.tags, minGap, plot = F){

  ## Function for finding r.bars.circular per frequency
  frequency.rbar <- function(detections, frequency){
    dt <- 1/frequency
    ## Convert detections into radians relative to sample frequency
    detections.mod.rad <- (detections %% dt)*2*pi/dt
    #browser()
    histo <- hist(detections.mod.rad,breaks = seq(0,2*pi,2*pi/100),plot = F)
    return(max(histo$counts))
  }

  Freqs <- data.frame(Frequency = frequencies)
  Freqs$r.bar <- sapply(X = Freqs$Frequency,
                        FUN = function(x){frequency.rbar(detections,x)})

  # Detect 4 highest peaks
  #browser()
  # List all peaks
  peaks.idx <- which(diff(sign(diff(Freqs$r.bar))) < 0) + 1
  # Order peaks by histogram values
  peaks.order <- order(Freqs$r.bar[peaks.idx], decreasing = T)

  ## Remove peaks which are too close to each other
  i <- 1
  while((i < n.tags) & (length(peaks.order) > n.tags)){
    idx <- i + which(abs(frequencies[peaks.idx[peaks.order[i]]] - frequencies[peaks.idx[peaks.order[(i+1):length(peaks.order)]]]) <= minGap)
    if(length(idx) > 0){
      ## remove clustered peaks
      peaks.order <- peaks.order[-idx]
    }
    i <- i + 1
    head(frequencies[peaks.idx[peaks.order]])
  }

  peaks.idx.filt <- peaks.idx[peaks.order[1:n.tags]]
  frequencies.peak <- Freqs[peaks.idx.filt,1]

  if(plot == T){
    require(ggplot2)
    print(
      ggplot(Freqs) + geom_line(aes(x = Frequency, y = r.bar)) +
        geom_point(data = Freqs[peaks.idx.filt,],
                   aes(x = Frequency, y = r.bar),
                   col = 'red') +
        theme_bw() + ylab('Peak Histogram Density') +
        ggtitle('Frequency Detection'))
  }

  return(frequencies.peak)
}

#' Filter pulse onsets by frequency
#'
#'This functions returns the indicies of those pulse detections which match with the provided pulse period frequencies.
#'Set \code{plot = T} to show filtering results in ggplot.
#'
#' @param detections A vector of detection times (seconds).
#' @param frequencies A vector of pulse onset frequencies (Hz) to filter.
#' @param sensitivity Bandwidth in radians around target frequency.
#' @param plot Set to T to plot results with ggplot
#' Increase this value, if you are missing detections due to the tag changing frequency over time.
#'
#' @export
#'
TagFreq.filter <- function(detections, frequencies,
                           sensitivity = 0.1,
                           plot = F){
  ## Add column for rads
  detections <- data.frame(seconds = detections,
                           rads = NA)

  ## Initalize results list
  results <- list()
  results.plot <- data.frame()

  ## Loop through each frequency
  for(i in frequencies){

    # Calc tag pulse period
    dt <- 1/i

    maxError.rad <- sensitivity

    ## Convert seconds to radians
    detections$rads <- (((detections$seconds %% dt)/dt)*2*pi)

    ## Find radian value for expected pulse
    # Generate histogram
    histo <- hist(detections$rads, breaks = 100, plot = F)
    # Check for time shift
    ## Detect 2 highest peaks
    peaks.order <- order(histo$density, decreasing = T)

    rad.peak <- histo$mids[order(histo$density,decreasing = T)[1]]

    ## Filter frequency of interest (get indicies for matches)
    window <- ((rad.peak+2*pi) + c(-1,1)*maxError.rad/2) %% (2*pi)
    if(window[1] > window[2]){ # Around 1,2
      idx <-  which(detections$rads >= window[1] |
                      detections$rads <= window[2])
    }else{ # Between 1,2
      idx <- which((detections$rads >= window[1]) &
                     (detections$rads <= window[2]))
    }

    ## Append results list
    results[[length(results) + 1]] <- idx

    if(plot == T){
      temp <- cbind(Freq = i, Detected = F, detections)
      temp$Detected[idx] <- T
      results.plot <- rbind(results.plot,temp)
    }
  }

  if (plot == T){
    ggtitle.TagFreq <- paste("Frequency Filtering | hydrophone",
                            mean(dataset.sub$Hydrophone), sep = " ")
    print(
      ggplot(results.plot) +
        geom_point(aes(x = seconds, y = rads, color = factor(Detected)),
                   alpha = 0.5) +
        theme_bw() + coord_polar(theta = 'y') + ylim(0,2*pi) +
        facet_wrap(~factor(Freq), ncol = 2) +
        scale_color_manual(values = c('black','red')) +
        guides(color = F) + ylab('Radians') + xlab('Seconds') +
        ggtitle(ggtitle.TagFreq))
  }

  return(results)

}

#' Label detections according to pulse period
#'
#' @param data data.frame with columns `c("Tag", "Hydrophone", "Seconds")`.
#' @param plot Show ggplot graphics for error checking segments.
#' @param sensitivity An optional value in radians for plotting the sensitivity bandwidth in the resulting ggplot.
#' This is just a approximate visual aid for selecting filter sensitivities for the `Freq.filter` function and should not be interpreted as exact.
#'
#' @return
#' @export
TagFreq.label <- function(data, plot = F, sensitivity = NA){
  ## Initialize results dataframe
  dataset <- data.frame()

  ## Get tag frequencies
  tag.freqs <- sort(unique(data$Tag))

  # Loop different tag frequencies
  for(j in tag.freqs){

    ## Make table for tag/hydrophone combination
    temp <- subset(data, subset = Tag == j)

    ########################################
    # Convert to normalized circular domain
    ########################################
    ## Calc pulse period
    dt <- 1/j
    ## Calculate raw periods
    period <- (temp$seconds / dt)
    ## Convert seconds to radians (relative to pulse period)
    rads.raw.unwrapped <- temp$Seconds*pi*2/dt
    rads.raw <- rads.raw.unwrapped %% (2*pi)
    ## Calc correction for moving rads to period center
    rads.center.const <- (-((circular.mean(rads.raw,
                                           weight = NULL,
                                           verbose = F)$theta + 2*pi) %%
                              (2*pi)) + 3*pi) %% (2*pi)
    ## Center rads around pi
    rads.center <- (rads.raw + rads.center.const) %% (2*pi)
    rads.center.unwrapped <- rads.raw.unwrapped + rads.center.const
    ## Get periods (starting at 1 for relative periods)
    period.abs <- floor(rads.center.unwrapped/(2*pi))
    period.rel <- period.abs - min(period.abs) + 1

    # Save values
    temp$period.start <- min(period.abs)*(dt) - ((rads.center.const/(2*pi))*dt)
    temp$dt <- dt
    temp$period <- as.integer(period.rel)
    temp$rads <- rads.center

    ## Append to results table
    dataset <- rbind(dataset,temp)
  }

  if(plot == T){
    require(ggplot2)
    print(
      ggplot(dataset) +
        geom_path(aes(x = Seconds, y = rads, group = factor(period))) +
        geom_point(aes(x = Seconds, y = rads, color = factor(Hydrophone))) +
        theme_bw() + scale_color_discrete(name = 'Hydrophone') +
        ylab('Relative Detection Time (Radians per period)') +
        xlab('Absolute Detection time (Seconds)') +
        geom_hline(yintercept = c(pi + sensitivity/2, pi - sensitivity/2),
                   linetype = 2)
    )
  }

  return(dataset)
}

#' Acoustic Tag Localization: Smooth detections
#'
#' This will remove reflections and apply a loess smoothing to the dataset, interpolating missed points.
#' For detecting reflections, reflections are assumed to arrive shortly after the direct signal.
#' Any detection which has an arrival time within \code{sensitivity} an earlier detection is considered a reflection.
#' You can set \code{plot = T} to get a visaul overview of the reflection filtering.
#'
#' @param data.tags A object returned from \code{atl.label}
#' @param sensitivity Sensitivity of the reflection detection (radians).
#' @param span Parameter for LOESS smoothing.
#' @param plot Set to \code{T} to plot results with \code{ggplot}.
#' @param freq.range Vector holding minimum and maximum pulse periods to scan.
#' @param freq.resolution Vector holding the dt of the pulse period search.
#'
#' @return
#'
#' @examples
atl.clean <- function(data.tags, sensitivity = 0.02, span = 0.5, plot = T){

  ## Intialize results data.frame
  data.tags.filtReflect <- data.frame()
  data.smooth <- data.frame()

  for(f in unique(data.tags$Freq.id)){
    for(h in unique(data.tags$Hydrophone)){
      for(p in unique(data.tags$period)){
        ## Subset dataset
        temp <- subset(data.tags,
                       subset = Hydrophone == h &
                         Freq.id == f &
                         period == p)
        if (nrow(temp) > 0){
          temp$Reflections <- F
          ## Search for refelctions if multiple pings in period
          if(nrow(temp) > 1){
            ## find differences in radiants between detections in this period
            mat <- as.data.frame(t(combn(x = 1:nrow(temp),2)))
            mat$diff <- abs(temp$rads[mat[,1]] - temp$rads[mat[,2]])
            ## remove later detections from pairs of detections
            idx <- which(mat$diff < sensitivity)
            idx.remove <- c()
            for(i in idx){
              if(temp$rads[mat[i,1]] > temp$rads[mat[i,2]])
                idx.remove <- append(idx.remove,mat[i,1])
              else{
                idx.remove <- append(idx.remove,mat[i,2])
              }
            }
            ## Mark reflected pulses
            temp$Reflections[unique(idx.remove)] <- T
          }
          # Append to results array
          data.tags.filtReflect <- rbind(data.tags.filtReflect, temp)
        }
      }
    }
  }

  if(plot == T){
    require(ggplot2)
    print(

      ggplot(data.tags.filtReflect) +
        geom_smooth(data = subset(data.tags.filtReflect, subset = Reflections == F),
                    aes(y = rads, x = period), method = 'loess',span = span,
                    method.args = list(family = "symmetric")) +
        geom_vline(aes(xintercept = period), alpha = 0.2) +
        geom_point(aes(y = rads, x = period,
                       fill = 20*log10(Peak_Amplitude),
                       color = Reflections,alpha = Reflections
        ),
        pch = 21) +
        theme_bw() + scale_fill_gradientn(colors = rainbow(4)) +
        scale_color_manual(values = c('black','red')) +
        guides(fill = guide_colorbar(title = 'Peak Amplitude (dB)')) +
        facet_grid(factor(Freq.id) ~ factor(Hydrophone),
                   scales = 'free', space = 'free') +
        scale_alpha_manual(values = c(1,0.5)))
  }

  return(data.tags.filtReflect)
}


#' Detection and Identification of Pulses with Variable Periods
#'
#' For acoustic tag localization experiments where the unique tag identity is defined by the pulse period.
#' This is intended for processing datasets from experiments which used acoustically tagged fish surrounded by a hydrophone array.
#'
#' This function will:
#' 1. Automatically detect the periods of each unique tag.
#' 2. Assign each pulse detection to the matching individual.
#' 3. Filter out detected pulse reflections.
#' 4. Apply loess smoothing to the resulting dataset.
#'
#' To verify the parameters you chose, set \code{plot = T} to see a visual overview of the dataset as it's processed.
#'
#' @param data Dataset returned from \code{read.HTI.RAT} function.
#' @param n.tags NUmber of tags with a unique pulse period.
#' @param sensitivity Sensitivity in radians for pulse period filtering.
#' @param window Window size (seconds)
#' @param span Amount of smoothing to apply to loess curve.
#' @param plot Show ggplot of smoothed data as the analysis runs.
#' @param proof Show verbose plots and pause analysis for each plot.
#' This is intended for helping choose initial analysis parameters.
#' @param pulsePeriodDistance The minimum distance (index units) between pulse period frequencies.
#' @param time.shift.sensitivity A parameter used to detect shifts in pulse detections.
#' Higher values will make the analysis more sensitive to correcting time shifts.
#'
#' @return
#' @export
#'
#' @examples
tagDetectionPulse <- function(data = data.raw, n.tags = 4, sensitivity = 0.2,
                              window = 60*3, span = 0.3, plot = T, proof = F,
                              pulsePeriodDistance.idx = 10,
                              time.shift.sensitivity = 3){

  ## Get start and ent times of dataset
  seconds.range <- range(data.raw$seconds)

  ## Initalize data.frame for smoothed data
  smoothed <- data.frame()
  col.names <- c('hydrophone', 'pulsePeriod.id', 'pulsePeriod.hz',
                 'seconds', 'reflection', 'smoothed')

  ## Get start and end seconds for smooting wthin window segments
  ## We'll only save the middle 50th quantile of the smoothed data
  smooth.range <- window*c(0.25,0.75)

  # Loop through overlapping windows of dataset
  window.starts <- 1:(floor(seconds.range[2]/window)*2 - 1)
  cat('Analysing windows...\n')
  pb <- txtProgressBar(min = seconds.range[1],
                       max = seconds.range[2],
                       initial = seconds.range[1],style = 3, width = 20)
  ptm <- proc.time()
  window.start.seconds <- 0
  while((window.start.seconds + window) < seconds.range[2]){

    # Show progress
    setTxtProgressBar(pb = pb,value = window.start.seconds)
    cat(paste0(', Analysis at: ', round(window.start.seconds/60,digits = 1),
               ' minutes' ))
    time <- (proc.time() - ptm)[3]/60
    cat(paste0(', run-time: ', round(time,2), ' minutes, Speed: x',
               round(((window.start.seconds- seconds.range[1])/60)/as.double(time), digits = 2)))

    ## Find each tag frequency and filter
    data.tags <- atl.label(data = data,start = window.start.seconds,
                           duration = window, n.tags = n.tags,
                           sensitivity = sensitivity, plot = proof,
                           pulsePeriodDistance.idx = pulsePeriodDistance.idx,
                           time.shift.sensitivity = time.shift.sensitivity
    )

    ## Check for empty window segment
    if(!is.null(nrow(data.tags))){

      ## Mark reflections
      data.tags.clean <- atl.clean(data.tags = data.tags,plot = plot,
                                   sensitivity = sensitivity, span = span)

      ## Apply loess smooth to middle 50th percentile of dataset
      ## apply to each freq/hydro combination
      for(f in unique(data.tags.clean$Freq.id)){
        for(h in unique(data.tags.clean$Hydrophone)){
          # Find matching samples
          idx <- which(data.tags.clean$Hydrophone == h &
                         data.tags.clean$Freq.id == f)
          ## Get frequency
          freq <- data.tags.clean$Freq[idx[1]]
          dt <- 1/freq
          ## Get period range of data
          period.range <- range(data.tags.clean$period[idx])

          # Build loess model with all data (Tukey smoothing)
          ## One point for each period (expected pulse)
          model.loess <- loess(family = 'sym', model = T, span = span,
                               data = data.tags.clean[idx,],
                               formula = rads ~ period)
          predict.loess <- data.frame(period = period.range[1]:period.range[2])
          predict.loess$rad <- predict(model.loess,newdata = predict.loess)

          # Convert period/rads into seconds
          offset <- data.tags.clean$period.start[idx[1]]
          ## Conversion factor for converting radians to seconds
          predict.loess$seconds <- (predict.loess$period)*dt + offset +
            (predict.loess$rad*dt/(2*pi))

          # Append results to returning data.frame
          ## Filter out values to export
          idx2 <- which(
            (predict.loess$seconds >= (window.start.seconds + (window*0.25)) &
               (predict.loess$seconds < (window.start.seconds) + (window*0.75)))
          )
          ## Check that there are values within window center
          if(length(idx2) > 0){
            temp <- data.tags.clean[idx, c('Hydrophone','Freq.id','Freq',
                                           'seconds', 'Reflections')]
            temp$Smoothed <- F
            names(temp) <- col.names
            ## Append raw data to results
            smoothed <- rbind(smoothed,
                              temp)
            names(smoothed) <- col.names
            ## Append smoothed data to results
            smoothed <- rbind(smoothed,
                              data.frame(hydrophone = as.integer(h),
                                         pulsePeriod.id = as.integer(f),
                                         pulsePeriod.hz = as.double(
                                           data.tags.clean$Freq[idx][1]),
                                         seconds = as.double(
                                           predict.loess$seconds[idx2]),
                                         reflection = F,
                                         smoothed = T
                              ))
          }
        }
      }
    }
    window.start.seconds <- window.start.seconds + floor(window/2)
  }
  class(smoothed) <- append('tagDetectionsPP',class(smoothed))
  return(smoothed)
}


#' Plot filtered and smoothed acoustic tag detections
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
plot.tagDetectionsPP <- function(x, ...){
  require(ggplot2)

  dataset <- x
  data.plot <- data.frame()

  ## Loop through hydro/Freq combinations
  for(f in unique(dataset$pulsePeriod.id)){
    for(h in unique(dataset$hydrophone)){
      ## Subset dataset
      idx <- which((dataset$pulsePeriod.id == f) &
                     (dataset$hydrophone == h))

      ## Get average tag periods (hz)
      hz <- mean(dataset$pulsePeriod.hz[idx])
      dt <- 1/hz

      ## Convert to radians
      temp <- dataset[idx,]
      temp$rad <- ((dataset$seconds[idx] %% dt)*2*pi) / dt

      # ## Get peak radian values and center values
      # histo <- hist(temp$rad, plot = F,breaks = 100)
      # peak.idx <- order(histo$density,decreasing = T)[1]
      # peak.rad <- histo$mids[peak.idx]
      # temp$rad <- (temp$rad - peak.rad + pi) %% (2*pi)
      data.plot <- rbind(data.plot,
                         temp)
    }
  }

  print(
    ggplot(data = data.plot) +
      geom_line(data = subset(data.plot, subset = smoothed == F &
                                reflection == F),
                aes(x = seconds/60, y = rad,
                    color = factor(pulsePeriod.id, ordered = F), alpha = 'real',
                    size = 'real', linetype = 'real')) +
      geom_line(data = subset(data.plot, subset = smoothed == T),
                aes(x = seconds/60, y = rad, col = factor(pulsePeriod.id, ordered = F),
                    alpha = 'smooth', size = 'smooth', linetype = 'smooth')) +
      facet_grid(hydrophone~.) + theme_bw() +
      scale_color_discrete(name = 'Tag') +
      scale_alpha_manual(values = c(real = 0.5, smooth =  1),
                         labels = c('Real', 'Smoothed'), name = 'Data') +
      scale_size_manual(values = c(real = 1, smooth =  0.5),
                        labels = c('Real', 'Smoothed'), name = 'Data') +
      scale_linetype_manual(values = c(real = '51', smooth =  'solid'),
                            labels = c('Real', 'Smoothed'), name = 'Data') +
      xlab('Minutes') + ylab('Pulse Period (radians)')
  )
}
