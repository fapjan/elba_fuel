## Simulation of woodland consumption on ancient Elba
## Authors: Fabian Becker, Raphael Eser
## Contact: fabian.becker@fu-berlin.de

## Freie Universitaet Berlin
## TOPOI A-5-4 "Iron mining and smelting on Elba Island"

## v1.1: 2017-10-07 Initial model
## v1.2: 2017-11-15 Phase included
## v1.3: 2017-11-27 some minor changes
## v1.4: 2017-11-30 including site specific production
## v1.5: 2017-12-22 minor changes in parameters, nicer plots, codeprep for publishing

## P_R_E_F_A_C_E #########################################################

## PACKAGES
library(rgdal)

## clean workspace
rm(list = ls())
start = Sys.time()

## local WORKING DIRECTORY
setwd("Y:/ELBA/paper/paper3_fuel/elba_fuel_code")

## READ DATA (current annual increment forests in Tuscany)
cai_raw <- read.table("cai.csv", header = TRUE, sep = ";")
cai <- as.vector(unlist(cai_raw))

## SETTINGS
set.seed(1814) # set 1814 for reproducability of results
n_iter = 1000 # number of Monte-Carlo simulations
ns_iter = rep(n_iter, 1)
bin_width = 100 # i.e. temporal resolution of the model

## DUMMIES
sim_data <- list()
area_km2 <- NULL
sensitivity <- NULL

## F_U_N_C_T_I_O_N_S #####################################################

## Function to simulate area necessary for wood consumption and wood production

simulate.area  <- function(x, steps = FALSE){
  
  ## Define vectors
  v <- vector(mode = "numeric", length = length(x$slg))
  area_km2_thinning <- v * NA
  area_km2_thinning_phase <- v * NA
  area_km2_thinning_site_phase <- v * NA
  area_km2_clearing <- v * NA
  area_km2_clearing_phase <- v * NA
  area_km2_clearing_site_phase <- v * NA
  area_km2_clearing2 <- v * NA
  area_km2_clearing_phase2 <- v * NA
  area_km2_clearing_site_phase2 <- v * NA
    consumption <- v * NA
  production_thinning <- v * NA
  production_clearfelling1 <- v * NA
  production_clearfelling <- v * NA
  consumption_phase <- v * NA
  consumption_site_phase <- v * NA
  
  ## Estimate consumption/production for each simulation run
  for (i in 1:length(x$slg)) {
    
    ## calculate charcoal consumption from estimated slag amount and furnace efficiency (slag:charcoal)
    charcoal = x$slg[i]  * x$furnace[i]
    
    ## Estimate wood necessary for the production of the charcoal required
    ## using a kiln efficiency factor (charcoal:wood ratio)
    wood = charcoal / x$kiln[i]
    
    ## the same for a phase of high production
    ## multiplying charcoal need with proportion of slag deposited in phase of high production
    wood_phase = (charcoal * x$si[i]) / x$kiln[i]
    
    ## same for differnt approach to obtain proportion in high production phase
    wood_site_phase = (charcoal * x$sis[i]) / x$kiln[i]
    
    ## annual wood consumption
    wood_yr = wood / x$duration[i]
    wood_yr_phase = wood_phase / x$i[i]
    wood_yr_site_phase = wood_site_phase / x$i[i]
    
    ## calculate standing weight per roation period from increment and length of rotation period
    stand_weight = x$rotation[i] * x$increment[i]
    stand_weight2 = x$rotation[i] * x$increment2[i]
    
    ## estimate area under management from annual wood demand and annual sustainable regrowth
    ## for smelting period and high production phases
    ## (coppicing assumption)
    area = wood_yr / x$productivity[i]
    area_phase = wood_yr_phase / x$productivity[i]
    area_site_phase = wood_yr_site_phase / x$productivity[i] 
    
    ## calulate area from standing weight (thinning assumption)
    area_yr = wood_yr / stand_weight
    area_yr2 = wood_yr / stand_weight2
    area_yr_phase = wood_yr_phase / stand_weight
    area_yr_phase2 = wood_yr_phase / stand_weight2
    area_yr_site_phase = wood_yr_site_phase / stand_weight
    area_yr_site_phase2 = wood_yr_site_phase / stand_weight2
    
    ## include non-use of area under recovery
    area_cycle = area_yr * x$rotation[i]
    area_cycle2 = area_yr2 * x$rotation[i]
    area_cycle_phase = area_yr_phase * x$rotation[i]
    area_cycle_phase2 = area_yr_phase2 * x$rotation[i]
    area_cycle_site_phase = area_yr_site_phase * x$rotation[i]
    area_cycle_site_phase2 = area_yr_site_phase2 * x$rotation[i]
    
    ## write estimated data from one simulation run in a vector with results of all simulation runs
    consumption[i] <- wood_yr
    consumption_phase[i] <- wood_yr_phase
    consumption_site_phase[i] <- wood_yr_site_phase
    production_thinning[i] <- (x$area[i] * 100) * x$productivity[i]
    production_clearfelling1[i] <- (x$area[i] * 100) * x$increment[i] 
    production_clearfelling[i] <- (x$area[i] * 100) * x$increment2[i] 
    
    ## change units
    area_km2_thinning[i] = area * 0.01
    area_km2_thinning_phase[i] = area_phase * 0.01
    area_km2_thinning_site_phase[i] = area_site_phase * 0.01
    
    area_km2_clearing[i] = area_cycle * 0.01
    area_km2_clearing_phase[i] = area_cycle_phase * 0.01
    area_km2_clearing_site_phase[i] = area_cycle_site_phase * 0.01
    
    area_km2_clearing2[i] = area_cycle2 * 0.01
    area_km2_clearing_phase2[i] = area_cycle_phase2 * 0.01
    area_km2_clearing_site_phase2[i] = area_cycle_site_phase2 * 0.01
    
  }
  
  ## add all estimated data to al list
  
  ## area
  if (steps == FALSE) {
    list(thinning = area_km2_thinning,
         thinning_phase = area_km2_thinning_phase,
         thinning_site_phase = area_km2_thinning_site_phase,
         clearing = area_km2_clearing,
         clearing_phase = area_km2_clearing_phase,
         clearing_site_phase = area_km2_clearing_site_phase,
         clearing2 = area_km2_clearing2,
         clearing_phase2 = area_km2_clearing_phase2,
         clearing_site_phase2 = area_km2_clearing_site_phase2
    )
  
    ## annual consumption/production  
  } else {
    list(consumption = consumption,
         consumption_phase = consumption_phase,
         consumption_site_phase = consumption_site_phase,
         production_thinning = production_thinning,
         production_clearfelling  = production_clearfelling1,
         production_clearfelling2 = production_clearfelling
    )
  }
  
}


## P_A_R_A_M_E_T_E_R #####################################################

## SLAG ##################################################################

## values for all ancient sites
## no data available for the Seccheto site
slg_wt <- read.table('slg_wt.csv', header = TRUE, sep = ";")

## randomly select the amount of slag for each site from the given possible range
## sum data from all sites
## repeat n-times

for (k in 1:n_iter) {
  for (i in 1:nrow(slg_wt)) {
    slg_wt$runif[i] <- runif(n = 1,
                             min = min(as.numeric(slg_wt[i,2:3])),
                             max = max(as.numeric(slg_wt[i,2:3]))
    )
  }
  sim_data$slg[k] <- sum(slg_wt$runif)
}



## SLAG:CHARCOAL #########################################################

## Data obtained by literature review
slg_to_chr_lit <- c(
  1/0.75,          # Wertime 1983 cit Pleiner, cf. Williams XXXX
  75/30,           # ca. Benvenuti 2016
  4,               # Nihlen, field notes, 1958
  2.5/0.4,         # Mommerstegg 2011 cit Voss 1988 cit Tylecote (no 19)
  5,               # Joosten (mass bilance)
  6,               # Crew,
  (3 + 1/3) / 0.67 # Thomas and Young (mass bilance)
)

## randomly sample from literature values
## assuming unifrom distribution of possible values between max and min of literature values
sim_data$furnace <- runif(n = n_iter,
                          min = min(slg_to_chr_lit),
                          max = max(slg_to_chr_lit)
)


## CHARCOAL:WOOD #########################################################

## charcoal recovery
## or: wood:charcoal ratio
## or: "charcoal mound productivity"
## or: efficienty factor
## adapted from VEAL 2012
## and VEAL 2017
## and SAREDO PARODI 2013

## read compilation of data from literature
chr_to_wood_lit <- read.table("k_eff.csv", header = TRUE, sep = ";")

## calculate efficency
chr_to_wood_lit_data <- chr_to_wood_lit$charcoal / chr_to_wood_lit$wood

## sample from efficiency factors
## assumption: factor normally distributed between min and max of range from literature
## range covers 3 sigma of true values
sim_data$kiln <- rnorm(n = n_iter,
                       mean = (max(chr_to_wood_lit_data) + min(chr_to_wood_lit_data)) / 2,
                       sd = (max(chr_to_wood_lit_data) - min(chr_to_wood_lit_data)) / 6
)


### SMELTING PERIOD ######################################################

## read data base on availabe dating material related to smelting activities
## based on a literature review
pottery_raw <- read.table("dates_per_sites.csv", 
                          sep = ";",
                          header = TRUE
)

## okay, it is a bit complicated, but ...
## replace empty character columns by NA
## necessary for further processing
pottery_raw$from <- sub("^$", NA, pottery_raw$from)
pottery_raw$to <- sub("^$", NA, pottery_raw$to)

## formatting of dates
## from XXX (B)CE to +/-
for (i in 1:length(pottery_raw$from)) {
  tmp <- unlist(strsplit(as.character(pottery_raw$from[i])," "))
  tmp <- as.data.frame(t(tmp))
  if (is.na(tmp[1]) == TRUE) {
    pottery_raw$start[i] <- NA
  } else {
    sign_factor <- ifelse(tmp[2] == "BCE",-1,1)
    pottery_raw$start[i] <- as.numeric(as.character(tmp[1,1])) * sign_factor
    
  }
}
for (i in 1:length(pottery_raw$to)) {
  tmp <- unlist(strsplit(as.character(pottery_raw$to[i])," "))
  tmp <- as.data.frame(t(tmp))
  if (is.na(tmp[1]) == TRUE) {
    pottery_raw$end[i] <- NA
  } else {
    sign_factor <- ifelse(tmp[2] == "BCE",-1,1)
    pottery_raw$end[i] <- as.numeric(as.character(tmp[1,1])) * sign_factor
    
  }
}
is <- order(pottery_raw$start)

## TIME OF OPERATION
## random sampling of deposition year for each dating material

## dummy list to be filled with simulated deposition dates
sim_dates <- vector("list", n_iter)

## run simulation k times
for (k in 1:n_iter) {
  ## simulate depostion date for run time i
  for (i in 1:nrow(pottery_raw)) {
    tmp <- runif(n = 1,
                 min = pottery_raw$start[i],
                 max = pottery_raw$end[i]
    )
    
    tmp <- round(tmp, 0)
    sim_dates[[k]][i] <- tmp
  }
  #print(paste("...",k,"/",n_iter," ..."))
}

## DURATION

## prepare simulated dates for duration calculation
## exclude oldest 2.5% and youngest 2.5% sherds from calculation
sim_q_95 <- lapply(sim_dates, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
sim_duration <- unname(unlist(lapply(sim_q_95, diff, na.rm = TRUE)))
sim_data$duration <- sim_duration


## DURATION SITES

## list of sites
sites <- unique(pottery_raw$site)

## find MINimum data of each site per simulation run
x <- aggregate(sim_dates[],list(pottery_raw$site), min, na.rm = TRUE)
x <- unname(x) ## proper (= no) naming of list

## lower minimum value per site
x_min <- as.data.frame(apply(x[,-1],1,quantile,probs = 0.025))

## find MAXimum data of each site per simulation run
x <- aggregate(sim_dates[],list(pottery_raw$site), max, na.rm = TRUE)
x <- unname(x) ## proper (= no) naming of list

## upper maximum value per site
x_max <- as.data.frame(apply(x[,-1],1,quantile,probs = 0.975))

## duration of each site 
duration = x_max - x_min

## proper data fromatting
site_data_dates <- cbind(x[[1]], x_min, x_max, duration)
colnames(site_data_dates) <- c("site", "start", "end", "duration")

## merge site duration and site slag
site_data <- merge(site_data_dates,slg_wt, all.y = TRUE)

## formatting  
site_data$end <- floor(site_data$end)
site_data$start <- floor(site_data$start)

## calculate annual production from site duration and site slag
site_data$from_yr <- site_data$from / site_data$duration
site_data$to_yr <- site_data$to / site_data$duration

## dummys 
years <- seq(-600,300,1) ## dummy vector 
slg_years <- matrix(nrow = length(years), ncol = nrow(site_data)) 

## formatting
site_data <- site_data[-which(is.na(site_data$start) == TRUE),]

for (m in 1:n_iter) {
  for (i in 1:nrow(site_data)) {
    slg <- runif(n = 1,
                 min = site_data$from_yr[i],
                 max = site_data$to_yr[i]
    )
    slg_years[which(years == site_data[i,]$start):which(years == site_data[i,]$end),i] <- slg
  }
  site_slg_ct <- aggregate(rowSums(slg_years, na.rm = TRUE), 
                           by = list(round(years/100,0)*100), 
                           sum)
  sim_data$sis[m] <- max(site_slg_ct$x / sum(site_slg_ct$x))
}

## Intensive period (duration)
sim_data$i <- rep(bin_width,n_iter)

## SLAG INTENSIV PERIOD (factor)
## find bin boundaries
## lower left boundary
bin_1 <- floor(min(pottery_raw$start, na.rm = TRUE) / 100) * 100

## upper right boundary
bin_n <- ceiling(max(pottery_raw$start, na.rm = TRUE) / 100) * 100

## get sequence of breaks
bin_breaks <- seq(bin_1, bin_n, bin_width)

## mid value of each bin
bin_mid <- seq(bin_1 + bin_width / 2, bin_n + bin_width / 2, bin_width)

## BIN results of each pottery date simulation run
## dummy list

sim_bin <- vector("list", n_iter)

## bin!
for (i in 1:length(sim_dates)) {
  sim_bin[[i]] <- cut(sim_dates[[i]], 
                      breaks = bin_breaks
  )
}
sim_bin_table <- lapply(sim_bin, table)

## data naming
bin_names <- levels(sim_bin[[1]])

## dummy data frame / naming
df <- data.frame(matrix(unlist(sim_bin_table), nrow = n_iter, byrow = T))
colnames(df) <- bin_names

## calculate statistics of relativ depostion of dating material
df_rel <- df / max(rowSums(df))
df_rel_mean <- apply(df_rel, 2, mean)
df_rel_sd <- apply(df_rel, 2, sd)

phase_intensity <- unname(unlist((df_rel[which(df_rel_mean == max(df_rel_mean))])))
sim_data$si <- phase_intensity


## FORST PRODUCTIVITY ####################################################

## data from (regional) literature
for_prod_lit <- c(1, ## from Veal: From context to economy
                  4#,  ## from Veal: From context to economy
                  #5, ## Craddock 1995 in Saredo Parodi 2013
                  #1.4, ## Sperl 2012 in Saredo Parodi 2013
                  #3.5 ##  Saredo Parodi 2013(cultivated Med. macchia forest)
)

## sample randomly from given literature values
## assuming that true value is between max and min
sim_data$productivity <- runif(n = n_iter,
                               min = min(for_prod_lit),
                               max = max(for_prod_lit)
)

## GROSS ANNUAL INCREMENT
gai_lit <- c(2,10) # Price 2000

sim_data$increment <- runif(n = n_iter,
                            min = min(gai_lit),
                            max = max(gai_lit)
)

## REGIONAL SPECIFIC ANNUAL INCREMENT
cai_raw <- read.table("cai.csv", header = TRUE, sep = ";")
cai <- as.vector(unlist(cai_raw))

## obtain propability function from CAI tables
increment2_pdf <- density(cai,
                          kernel = "rectangular",
                          from = min(gai_lit),
                          to = max(gai_lit)
)
increment2 <-  approx(cumsum(increment2_pdf$y) / sum(increment2_pdf$y),
                      increment2_pdf$x,
                      n = n_iter)$y
sim_data$increment2 <- increment2[sample(seq(1, n_iter,1), size = n_iter)]

## ROTATION PERIOD #######################################################

## Roation period as reported by Pliny tE and Columella 
rotation_lit <- c(5, 10)

sim_data$rotation <- runif(n = n_iter,
                           min = min(rotation_lit),
                           max = max(rotation_lit)
)



## DURATION OF OPERATION (t)

## processing see above, here: plotting

## plot raw chronology data
par(mfrow = c(1,1))
plot(1,
     type = "n",
     xlim =c (min(pottery_raw$start, na.rm = TRUE),
            max(pottery_raw$end,na.rm = TRUE)),
     ylim = c(1,
            nrow(pottery_raw))
)
for (i in 1:nrow(pottery_raw)) {
  segments(y1 = i,y0 = i,
           x0 = pottery_raw$start[is][i],
           x1 = pottery_raw$end[is][i]
  )
}

## AREA OF ISLAND (km2) #################################################


wood_shape <- readOGR(dsn = ".", layer = "060101_bosco_3003")
nonwood_shape <- readOGR(dsn = ".", layer = "060105_ps_inc_3003")

## Getting wood data
wood_data <- wood_shape@data
nonwood_data <- nonwood_shape@data

wood_area <- as.data.frame(aggregate(wood_data$area, list(wood_data$bosco_dty), sum)[,2] / 1000000)
nonwood_area <- as.data.frame(aggregate(nonwood_data$area, list(nonwood_data$ps_inc_dty), sum)[,2] / 1000000)

## some formatting
colnames(wood_area) <- "area_km2"
colnames(nonwood_area) <- "area_km2"
rownames(wood_area) <- aggregate(wood_data$area, list(wood_data$bosco_dty), sum)[,1]
rownames(nonwood_area) <- aggregate(nonwood_data$area, list(nonwood_data$ps_inc_dty), sum)[,1]

## woodland area from present vegetation cover and (conservative) estimations from literature
area_lit <- c(
  sum(wood_area),
  sum(wood_area) + sum(nonwood_area),
  224 * 0.40, ## Saredo Parodi 2013
  100 ## Brambila in in Saredo Parodi 2013:  base of calculation?
)

## sample from above values, assuming normal distrbution
sim_data$area <- rnorm(n = n_iter,
                       mean = (max(area_lit) + min(area_lit)) / 2,
                       sd = (max(area_lit) - min(area_lit)) / 6
)


### PARAMETER PLOT #######################################################

# Parameter naming for plotting
parameters_main <- c("amount of slag",
                     "furnace efficiency",
                     "kiln efficiency",
                     " duration smelting period (A1)",
                     "proportion intensive period (A2)",
                     "duration intensive period",
                     "proportion intensive period (A3)",
                     "productivity",
                     "increment",
                     "adjusted increment",
                     "rotation period",
                     "harvestable area"
)

## Plotting
par(mfrow = c(3,4), mar = c(3,3,2,0.1))
kg <- function(x,data) return(mean(dnorm((x - data)/bw.nrd0(data)))/bw.nrd0(data))
plot_data <- which(!(names(sim_data) %in% c("increment", "i")))
plot_data <- c(1,2,3,4,
               5,6,7,8,
               9,10,11,12
)
def_x_lim <- list()
for(i in plot_data){
  def_x_lim[[i]] <- c(min(sim_data[[i]]),
                      max(sim_data[[i]]))
  h <- hist(sim_data[[i]],
            main = parameters_main[i],
            xlab="", 
            ylab="", yaxt = "n",
            freq = FALSE,
            border = "darkgrey", col = "grey",
            xlim = def_x_lim[[i]]
  )
  lines(density(sim_data[[i]],
                from = min(sim_data[[i]]),
                to = max(sim_data[[i]])
  )
  )
  center <- median(sim_data[[i]])
  segments(x0 = center, x1 = center,
           y0 = 0, y1 = kg(center,sim_data[[i]])
  )
}

### SENSITIVITY ##########################################################

## dummy vector

y = sim_data
l = length(simulate.area(y))
sensitivity = matrix(nrow = length(sim_data), ncol = l)

## calculate mean of each parameter
## assign it to a vector
for (k in 1:length(sim_data)) {
  y[[k]] <- rep(mean(y[[k]]),length(y[[k]]))
}

## run simulation with one parameter changing
## other parameters are set to mean
for (m in 1:length(sim_data)) {
  for (i in 1:l) {
    z <- y
    z[[m]] <- sim_data[[m]]
    a <- simulate.area(z)[i]
    
    sensitivity[m,i] <- sd(a[[1]]) / mean(a[[1]])
  }
  
}
b <- z$area - a[[1]]
sensitivity[m,i] <- sd(b) / mean(b)

## data structuring  
sensitivity <- as.data.frame(sensitivity)

sensitivity <- rowMeans(replace(sensitivity,sensitivity == 0, NA), na.rm = TRUE)
names(sensitivity) <- parameters_main

## plotting
par(mar = c(14,3,0.1,0.1), mfrow=c(1,1))
barplot(sensitivity, las = 2)  


### CALCULATIONS AND PLOTTING ############################################

area_km2 <- simulate.area(sim_data)

hist_breaks <- seq(0,max(unlist(area_km2)) + 10,10)

result_title <- c(
"thinning, period (S1)", "thinning, phase (S2)", "thinning, phase (S3)",
"coppicing, period (S1)", "coppicing, phase (S2)", "coppicing, phase (S3)",
"coppicing, period (S1)", "coppicing, phase (S2)", "coppicing, phase (S3)"
)

par(mfrow = c(3,3), mar = c(3,3,4,0.1))
for (i in 1:length(area_km2)) {
  hist(area_km2[[i]],
       breaks = hist_breaks,
       xlim = c(0,quantile(unlist(area_km2),prob = 0.995)),
       ylim = c(0,0.05),
       border = "darkgrey", col = "grey",
       freq = FALSE,
       main = result_title[i]
  )
  lines(density(area_km2[[i]]), col = "black")
  lines(density(sim_data$area), col = "black")
  alpha = 0.15
  q_s <- quantile(area_km2[[i]], probs = c(alpha, 1 - alpha))
  seg_fac <- 0.0005
  segments(x0 = q_s[1] , x1 = q_s[2],
           y0 = 0 + seg_fac, y1 = 0 + seg_fac,
           lwd = 4, col = "black")
  q_s <- quantile(sim_data$area, probs = c(alpha, 1 - alpha))
  segments(x0 = q_s[1] , x1 = q_s[2],
           y0 = 0 - seg_fac, y1 = 0 - seg_fac,
           lwd = 4, col = "black")
  abline(v = 224, col = "grey")
  
}

## Descriptive statistiscs area consumption
str(area_km2)
area_consumption <- as.data.frame(unlist(lapply(area_km2, mean)))
colnames(area_consumption) <- "mean"
area_consumption$q2.5 <- unlist(lapply(area_km2,quantile,probs = 0.025, na.rm = TRUE))
area_consumption$q16.7 <- unlist(lapply(area_km2,quantile,probs = 1/6, na.rm = TRUE))
area_consumption$q83.3 <- unlist(lapply(area_km2,quantile,probs = 5/6, na.rm = TRUE))
area_consumption$q97.5 <- unlist(lapply(area_km2,quantile,probs = 0.975, na.rm = TRUE))
area_consumption$q95 <- unlist(lapply(area_km2,quantile,probs = 0.95, na.rm = TRUE))
area_consumption$q66.7 <- unlist(lapply(area_km2,quantile,probs = 4/6, na.rm = TRUE))

area_production <- as.data.frame(mean(sim_data$area))
colnames(area_production) <- "mean"
area_production$q2.5 <- quantile(sim_data$area, probs =  0.025)
area_production$q16.7 <- quantile(sim_data$area, probs =  1/6)
area_production$q66.7 <- quantile(sim_data$area, probs =  4/6)
area_production$q83.3 <- quantile(sim_data$area, probs =  5/6)
area_production$q95 <- quantile(sim_data$area, probs =  0.95)
area_production$q97.5 <- quantile(sim_data$area, probs =  0.975)

## Proportion of the island's area which is under management 
def_prop <- list()
for (i in 1:length(area_km2)) {
  def_prop[[i]] <- sum(sim_data$area - area_km2[[i]] < 0) / n_iter
}

## some more data output
overarea_prop <- as.data.frame(matrix(unlist(def_prop),ncol = 3, nrow = 3, byrow = TRUE))
rownames(overarea_prop) <- c("thinning","coppicing_unif","coppicing")
colnames(overarea_prop) <- c("period","phase","phase_sites")
format(overarea_prop, scientific = FALSE)


## calculation
wood_t <- simulate.area(sim_data, steps = TRUE)
str(wood_t)
sum(wood_t$production_thinning - wood_t$consumption < 0) / n_iter

## some addtional results
overproduction_thinning <- sum(wood_t$production_thinning - wood_t$consumption < 0) / n_iter
overproduction_clearfelling <- sum(wood_t$production_clearfelling - wood_t$consumption < 0) / n_iter
overproduction_thinning_phase_site <- sum(wood_t$production_thinning - wood_t$consumption_site_phase < 0) / n_iter
overproduction_clearfelling_phase_site <- sum(wood_t$production_clearfelling - wood_t$consumption_site_phase < 0) / n_iter

annual_wood <- as.data.frame(unlist(lapply(wood_t, mean)))
colnames(annual_wood) <- "mean"
annual_wood$q2.5 <- unlist(lapply(wood_t,quantile,probs = 0.025, na.rm = TRUE))
annual_wood$q16.7 <- unlist(lapply(wood_t,quantile,probs = 1/6, na.rm = TRUE))
annual_wood$q50.0 <- unlist(lapply(wood_t,quantile,probs = 3/6, na.rm = TRUE))
annual_wood$q83.3 <- unlist(lapply(wood_t,quantile,probs = 5/6, na.rm = TRUE))
annual_wood$q97.5 <- unlist(lapply(wood_t,quantile,probs = 0.975, na.rm = TRUE))