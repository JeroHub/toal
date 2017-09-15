require(ggplot2)
require(toal)

## Define a simulated sensor array.
sensorLocations <- data.frame(x = c(0,10, 10, 0),
y = c(10, 0, 10, 0),
z = c(5, 5, 0, 0),
row.names = paste0('H',1:4))


## Create 5 random source locations within array
n <- 7 # point sources to simulate
sourceLocation <- data.frame(x = runif(n = n,min = -10, max = 20),
y = runif(n = n,min = 0, max = 10),
z = runif(n = n,min = 0, max = 5))

## Calculate simulated arrival delays
N <- nrow(sensorLocations) # Number of sensors
c <- 1500 # Speed of sound in water

## Function for calculating distance between points
dist <- function(a,b){
## Calc distance between two locations
dif <- abs(a - b)
return(sqrt(sum(dif^2)))
}

# Initalize empty array of TOAs
sensorLatency = matrix(ncol = N,nrow = nrow(sourceLocation))
# Calculate distances from each
for(i in 1:N){
for(j in 1:nrow(sourceLocation)){
sensorLatency[j,i] <- dist(sourceLocation[j,],sensorLocations[i,1:3])
}
}
rownames(sensorLatency) <- paste0('Source ',1:nrow(sourceLocation))
colnames(sensorLatency) <- paste0('Hydrophone ',1:N)

## Convert distances to arrival latencies (from meters to seconds)
sensorLatency <- sensorLatency/c

## Plot resulting source locations and location estimates
nudge <- 0.5


sourceEstimates <- TOA.localization(toa = sensorLatency,
hydrohpone.positions = sensorLocations,
c = 1500)

ggplot(sensorLocations) +
geom_point(aes(x = x, y = y, size = z)) +
geom_point(data = sourceLocation,
aes(x = x, y = y, size = z),
color = 'red', pch = 4) +
geom_text(aes(x = x, y = y, label = rownames(sensorLocations)),
nudge_x = nudge, nudge_y = nudge) +
theme_bw() + coord_equal() +
scale_size_continuous(name = 'Depth',trans = 'reverse') +
ggtitle('Simulated Sensor array and sound sources') +
geom_point(data = sourceEstimates,
           aes(x = x, y = y, size = z, color = factor(eq)),
           pch = 3, stroke = 1) +
scale_color_discrete(name = 'Equation')

