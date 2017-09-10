### Time-Of-Arrival Localization

`toal` is a package intended for biologists working with tagged animals monitored by stationary sensor arrays (minimum 4 sensors).
Using time-of-arrival differences, `toal` implements the spherical interpolation method of localizing sound sources.
This method is closed form, and is thus quite quick to process.

Additionally, this package provides functions geared towards marine biologisits which implement methods for:
- identifying acoustic tags based on the pulse period
- Automatically removing surface reflections
- Reading raw data from HTI acoustic tag monitoring systems
- Applying LOESS smoothing to acoustic pulse detections and estimating in missing detections

To get a quick look at what this package does, you can install it and load the example vignettes with the code below:

```r
# Install from github
require(devtools)
install_github('RTbecard/toal', build_vignettes = T, force = T)

# Open vignette for sound source localization
vignette(package = 'toal', topic = 'Localization')
```
