require(toal)

# Load data
## samplerate of HTI dataset = 12000Hz
# Find 4 tag frequencies
#data.raw <- read.HTI.RAT('~/Repos/HydrophoneLocalization/HTI_Example_Data/AT16S012840914 (G8 4,0,48).RAT',
#                         fs = 12000)

data.raw <- read.HTI.RAT('~/Dropbox/PCAD4Cod/HTI data for James/AT16S012840914 (G8 4,0,48).RAT',
                         fs = 12000)

data.smooth <- tagDetectionPulse(data = data.raw, n.tags = 4, sensitivity = 0.2,
                                 window = 3*60, span = 0.2, pulsePeriodDistance.idx = 10,
                                 plot = T)

# data.smooth <- tagDetectionPulse(data = subset(data.raw,
#                                                subset = seconds > 0 &
#                                                  seconds < 10*60),
#                                  n.tags = 4, sensitivity = 0.2,
#                                  window = 3*60, span = 0.2, plot = T,
#                                  proof = F,
#                                  pulsePeriodDistance.idx = 10,
#                                  time.shift.sensitivity = 3)

plot(data.smooth)
