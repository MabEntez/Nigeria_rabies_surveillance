rm(list=ls())  #Erase workspace
library(ggplot2)
library(spdep)
library(dplyr)
library(tidyr)
library(stringr)
library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))

##### Preparing Data ######

dt <- read.csv("Nigeria data.csv")
dt <- dt %>%
  filter(Spp == "Dogs")
dt[dt==""] <- NA
dt %>%
      group_by(Year) %>%
         summarise(Total_Cases = sum(Number.of.outbreaks, na.rm = TRUE),
                                     Cases = sum(Cases, na.rm = TRUE),
                                     Death = sum(Death, na.rm = TRUE),
                                     Killed.and.disposed.off = sum(Killed.and.disposed.off, na.rm = TRUE))
dt <- dt %>%
  group_by(Year, Month, State) %>%
  summarise(Total_Cases = sum(Number.of.outbreaks, rm.na = TRUE))
month_mapping <- c(
  "Jan" = "January", "Feb" = "February", "Mar" = "March",
  "Apr" = "April", "May" = "May", "Jun" = "June",
  "Jul" = "July", "Aug" = "August", "Sep" = "September", "Sept" = "September",
  "Oct" = "October", "Nov" = "November", "Dec" = "December"
)
dt <- dt %>%
  mutate(Month = recode(Month, !!!month_mapping))


seasonality <- dt %>%
  group_by(Month) %>%
  summarise(Total_Cases = sum(Total_Cases, na.rm = TRUE))

seasonality <- seasonality %>%
  filter(!Month %in% c("Jan-Jun 2022", "Jan-Jun 2023", "Jul-Dec 2022"))

seasonality$Month <- factor(seasonality$Month, levels = c("January", "February", "March", "April", "May", 
                                            "June", "July", "August", "September", "October", 
                                            "November", "December"))

ggplot(seasonality, aes(x = Month, y = Total_Cases)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  labs(title = "The number of outbreaks of rabies per month from 2014 - 2023", x = "Month", y = "Number of outbreaks") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

dt <- dt %>%
  mutate(Year_Range = case_when(
    Year >= 2014 & Year <= 2016 ~ '2014-2016',
    Year == 2017 ~ '2017',
    Year == 2018 ~ '2018',
    Year == 2019 ~ '2019',
    Year == 2020 ~ '2020',
    Year == 2021 ~ '2021',
    Year == 2022 ~ '2022',
    Year == 2023 ~ '2023',
    TRUE ~ as.character(Year)
  ))

outbreak <- dt %>%
  group_by(State, Year_Range) %>%
  summarize(Total = sum(Total_Cases, na.rm = TRUE)) %>%
  ungroup()

# Reshape the data into wide format
outbreak <- outbreak %>%
  pivot_wider(names_from = Year_Range, values_from = Total, values_fill = NA)

outbreak[outbreak == "FCT"] <- "Federal Capital Territory"
outbreak[outbreak == "Nassarawa"] <- "Nasarawa"
##### Sort out adjacency matrix for States ####
nigeria <- st_read("./Nigeria Shapefile/nga_admbnda_adm1_osgof_20190417.shp") # Loading the Shape File
#summary(nigeria)
#plot (nigeria)

states_list <- data.frame(State = unique(nigeria$ADM1_EN))

# Assuming `outbreak` is your reshaped data
# Perform a full join to ensure all states are included
outbreak <- full_join(states_list, outbreak, by = "State")

# Replace any missing values with NA
outbreak[is.na(outbreak)] <- NA

outbreak[,1] <- nigeria$ADM1_PCODE

map_o <- reshape(outbreak,
                     timevar = "Year",
                     times = c("2014-2016", "2017", "2018", "2019", "2020", "2021"),
                     varying = c("2014-2016", "2017", "2018", "2019", "2020", "2021"),
                     v.names = "outbreak",
                     idvar = "State",
                     direction = "long")
map_o <- merge(nigeria, map_o, by.x = "ADM1_PCODE", by.y = "State")

y <- ggplot(map_o) + geom_sf(aes(fill = outbreak)) +
  facet_wrap(~Year, dir = "h", ncol = 3) +
  ggtitle("Reported Outbreaks") + theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank()
  ) +
  scale_fill_gradient2(
    low = "white", high = "darkgreen"
  )+
  labs(fill = "Number of\nOutbreaks")

png("Reported Outbreaks.png", width = 600, height = 400)
y
dev.off()

nig_adj <- poly2nb(as(nigeria, "Spatial"), row.names = nigeria$ADM1_PCODE)
nb_adj <- nb2WB(nig_adj)
adj <- nb_adj$adj
Ladj <- length(adj)
weightsAdj <- nb_adj$weights
numAdj <- nb_adj$num

nstate <- length(nigeria$ADM1_PCODE)
ntime <- length(1:8)
nat_avg <- mean(as.matrix(outbreak[,2:9]), na.rm = TRUE)
## Set all zero reports to NA
##caseDog[caseDog==0] <- NA

save.image("NimbleDataCom.RData")



##### Nimble Model For Dogs ######
library(rstudioapi)
library(nimble)
setwd(dirname(getActiveDocumentContext()$path)) #set file location as working directory
rm(list=ls())  #Erase workspace

load("NimbleDataCom.RData")



modelCode <- nimbleCode({
  for (state in 1:nstate){
    for (t in 1:ntime){
      outbreak[state,t] ~ dpois(muDog[state,t])
      log(muDog[state,t]) <- intercept + u[state] + v[state] + g[t] + d[state,t] 
      
      d[state,t] ~ dnorm(0,taud)
    }	
    v[state] ~ dnorm(0,tauv)
  }
  
  g[1] ~ dnorm(0,taug)
  for (t in 2:ntime){
    g[t] ~ dnorm(g[t-1],taug)
  }
  
  u[1:nstate]~dcar_normal(adj[1:Ladj],weightsAdj[1:Ladj],numAdj[1:nstate],tauu)
  
  intercept ~ dpois(nat_avg)
  #intercept <- 0
  tauu <- 1/pow(sdu,2)
  sdu ~ dunif(0,6)
  taug <- 1/pow(sdg,2)
  sdg ~ dunif(0,4)
  taud <- 1/pow(sdd,2)
  sdd ~ dunif(0,4)
  tauv <- 1/pow(sdv,2)
  sdv ~ dunif(0,2)
  
})

modelData <- list(outbreak=outbreak[,2:9])


modelConstants <- list(adj=adj, Ladj=Ladj, weightsAdj=weightsAdj, numAdj=numAdj, ## settings for adjacency matrix
                       nstate=nstate, ntime=ntime, nat_avg=nat_avg ## settings for loops
                       ) ## matrix of expectations of cases
modelInits <- list(intercept=0, sdu=.1,sdg=.1, sdd=.1, sdv = .1,
                   u=rep(0,nstate), v=rep(0,nstate), g=rep(0,ntime), d=matrix(0.1,nrow=nstate,ncol=ntime))

## Create model
model <- nimbleModel(code=modelCode,name="model",constants=modelConstants,
                     data=modelData,inits=modelInits)

## Compile model
Cmodel <- compileNimble(model)

## Configure MCMC
modelConfig <- configureMCMC(model,monitors=Cmodel$getNodeNames(stochOnly=T,includeData=F),
                             enableWAIC = T, print=T)

modelConfig$addMonitors(c("outbreak","intercept","g","u","v"))
## Build MCMC
modelMCMC <- buildMCMC(modelConfig)

## Compile MCMC
CmodelMCMC <- compileNimble(modelMCMC,project=model,resetFunctions = T)

## Run MCMC
set.seed(1)
samples <- runMCMC(CmodelMCMC,niter=30000,thin = 5,nburnin=28000,nchains=3,summary=T, WAIC=T,
                   samplesAsCodaMCMC=T, progressBar = T)

save.image("NimbleModelCom.RData")

##### Plotting Results as Test ######
rm(list=ls())  #Erase workspace
library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))
load("NimbleModelCom.RData")
library(dplyr)
library(tmaptools)
library(ggplot2)
library(tmap)
library(BAMMtools)
library(magrittr)
library(tmap)
library(sf)

samples$summary
samples$WAIC
plot(samples$samples[,"intercept"])
plot(samples$samples[,"u[1]"])
plot(samples$samples[,"v[1]"])
plot(samples$samples[,"d[1, 3]"])
plot(samples$samples[,"sdu"])
plot(samples$samples[,"sdg"])

plot(samples$samples[,"sdd"])
plot(samples$samples[,"sdv"])
plot(samples$samples[,"g[1]"])

plot(samples$samples[,"outbreak[1, 7]"])

##### Plotting V & U Mean #####
results_uv <- matrix(nrow = 37, ncol = 3)
colnames(results_uv) <- c("ADM1_PCODE", "Vi Mean", "Ui Mean")
for (i in 1:2){
  if (i == 1) {
    variation <- "v["
    col_num <- 2
  } else {
    variation <- "u["
    col_num <- 3
  }
  state_number <- 0
  for (state in 1:37){
    state_number <- state_number + 1
    input <- paste(variation, as.character(state_number) , "]", sep="")
    exp_data <- sapply(unlist(samples$samples[,input]), exp)
    mean_value <- mean(exp_data)
    results_uv[state_number, col_num] <- mean_value
    if (state_number <= 9){
      results_uv[state_number, 1] <- paste("NG00", as.character(state_number), sep="")
    } else {
      results_uv[state_number, 1] <- paste("NGY0", as.character(state_number), sep="")
    }
  }
}
write.csv(results_uv, file = "Mean V & U Dog.csv", row.names = F)
csv_results_uv = read.csv("Mean V & U Dog.csv")

jenks_v <- getJenksBreaks(csv_results_uv$Vi.Mean, 7, subset = NULL)
map_UV <- cbind(nigeria, csv_results_uv, deparse.level = 1) 
png("Unstructured.png", width = 500,height = 500)
qtm(map_UV, fill = "Vi.Mean",
    fill.palette = "Blues",
    fill.style="fixed",
    fill.breaks=c(jenks_v[1], jenks_v[2], jenks_v[3], jenks_v[4], jenks_v[5], jenks_v[6], jenks_v[7]),
    fill.title= "UnStructured\nVariance Mean") +
  tm_legend(legend.position = c("right", "bottom"))
dev.off()
jenks_u <- getJenksBreaks(csv_results_uv$Ui.Mean, 7, subset = NULL)

range_U <- range(csv_results_uv[,3])
multi_U <- (range_U[2] - range_U[1])/6
png("Structured.png", width = 500,height = 500)
qtm(map_UV, fill = "Ui.Mean",
    fill.palette = "Blues",
    fill.style="fixed",
    fill.breaks=c(jenks_u[1], jenks_u[2], jenks_u[3], jenks_u[4], jenks_u[5], jenks_u[6], jenks_u[7]),
    fill.title="Structured\nVariance Mean") +
  tm_legend(legend.position = c("right", "bottom"))
dev.off()

##### Plotting G Mean #####

time <- 1
results_g <- matrix(nrow = 8, ncol = 3)
year = 2015
for (time in 1:8){
  year <- year + 1
  input <- paste("g[", as.character(time) , "]", sep="")
  exp_data <- sapply(unlist(samples$samples[,input]), exp)
  mean_value <- mean(exp_data)
  sd_value <- sd(exp_data)
  results_g[time, 3] <- sd_value
  results_g[time, 2] <- mean_value
  results_g[time, 1] <- year
}
#results_g[1,1] <- "2014-2016"
write.csv(results_g, file = "Mean G Dog.csv", row.names = F)
csv_results_g = read.csv("Mean G Dog.csv")
par(mfrow = c(1, 1))
plot (csv_results_g[, -3], ylim=c(0.5,1.7), type = "l", col = "red", xlab = "Year", ylab = "Teporal Effect", main = "Temporal Effect Predicted by Model 2014 - 2021")
arrows(y0 = csv_results_g[,2] - csv_results_g[,3], y1 = csv_results_g[,2] + csv_results_g[,3], x0 = csv_results_g[,1], x1 = csv_results_g[,1] , code=3,angle=90, length=0.1, col="black", lwd=0.5)

##### Plotting D Mean #####
results_d <- matrix(nrow = 37, ncol = 9)
colnames(results_d) <- c("State Order", "2014 - 2016", "2017", "2018", "2019", "2020", "2021", "2022", "2023")
for (i in 1:5){
  state_number <- 0
  for (state in 1:37){
    state_number <- state_number + 1
    input <- paste("d[", as.character(state_number), ", ", as.character(i) , "]", sep="")
    summary_dt <- summary(samples$samples[,input])
    mean_value_log <- summary_dt$statistics[1]
    mean_value <- 10 ** mean_value_log
    results_d[state_number, (i + 1)] <- mean_value
    results_d[state_number, 1] <- state_number
  }
}
write.csv(results_d, file = "Mean D Dog.csv", row.names = F)
cvs_results_d = read.csv("Mean D Dog.csv")

##### Plotting Mu #####

results_mu <- matrix(nrow = 37, ncol = 9)
colnames(results_mu) <- c("State Order", "2014 - 2016", "2017", "2018", "2019", "2020", "2021", "2022", "2023")
for (i in 1:8){
  state_number <- 0
  for (state in 1:37){
    state_number <- state_number + 1
    input <- paste("outbreak[", as.character(state_number), ", ", as.character(i) , "]", sep="")
    summary_dt <- summary(samples$samples[,input])
    mean_value <- summary_dt$statistics[1]
    results_mu[state_number, (i + 1)] <- mean_value
    results_mu[state_number, 1] <- state_number
  }
}
write.csv(results_mu, file = "Mean Mu Dog.csv", row.names = F)
csv_results_mu = read.csv("Mean Mu Dog.csv")

csv_results_mu[, 2:9] <- round(csv_results_mu[, 2:9])
map_results_mu <- csv_results_mu
colnames(map_results_mu) <- c("State Order", "2014 - 2016", "2017", "2018", "2019", "2020", "2021", "2022", "2023")
jenks_mu <- getJenksBreaks(c(csv_results_mu[,2], csv_results_mu[,3], csv_results_mu[,4], csv_results_mu[,5], csv_results_mu[,6], csv_results_mu[,7], csv_results_mu[,8], csv_results_mu[,9]), 7, subset = NULL)
map_results_mu$province <- nigeria$ADM1_EN
map_results_mu <- reshape(map_results_mu,
                     timevar = "Year",
                     times = c("2014 - 2016", "2017", "2018", "2019", "2020", "2021", "2022", "2023"),
                     varying = c("2014 - 2016", "2017", "2018", "2019", "2020", "2021", "2022", "2023"),
                     v.names = "mu",
                     idvar = "province",
                     direction = "long")

map_mu <- merge(nigeria, map_results_mu, by.x = "ADM1_EN", by.y = "province")

mu <- ggplot(map_mu) + geom_sf(aes(fill = mu)) +
  facet_wrap(~Year, dir = "h", ncol = 4) +
  ggtitle("Predicted outbreaks") + theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank()
  ) +
  scale_fill_gradient2(
    low = "white", high = "blue", name = "Predicted\nnumber\nof outbreaks"
  )
png("Model predicted outbreaks.png", width = 900, height = 400)
mu
dev.off()

library(classInt)

# Create a color palette based on the number of Jenks breaks
breaks <- jenks_mu
num_colors <- length(breaks) - 1
colors <- colorRampPalette(c("white", "red"))(num_colors)

# Create a new column in your data frame to categorize 'mu' based on Jenks breaks
map_mu$mu_category <- cut(map_mu$mu, breaks = breaks, include.lowest = TRUE, labels = FALSE)

# Map the categories to colors
map_mu$mu_color <- colors[map_mu$mu_category]

labels <- sapply(1:num_colors, function(i) {
  paste0(breaks[i], "-", breaks[i + 1])
})

# Plotting with ggplot2
ggplot(map_mu) + 
  geom_sf(aes(fill = factor(mu_category))) +
  facet_wrap(~Year, dir = "h", ncol = 4) +
  ggtitle("Predicted outbreaks") + 
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank()
  ) +
  scale_fill_manual(
    values = colors,
    name = "Predicted outbreaks",
    breaks = seq_along(breaks[-1]), 
    labels = labels
  )

##### SIR Calculation #####
SIR_cases <- round((csv_results_mu[,-c(1)]))
colnames(SIR_cases) <- c("2014 - 2016", "2017", "2018", "2019", "2020", "2021", "2022", "2023")
# Calculate the mean of all cells in the selected columns
SIR_expected <- mean(as.matrix(modelData$outbreak), na.rm = TRUE)

sir_ratio <- (SIR_cases + 0.0001)/(SIR_expected+ 0.0001)
sir_ratio$province <- nigeria$ADM1_EN
jenks_sir <- getJenksBreaks(c(sir_ratio$`2014 - 2016`, sir_ratio$`2017`, sir_ratio$`2018`, sir_ratio$`2019`, sir_ratio$`2020`, sir_ratio$`2021`, sir_ratio$`2022`, sir_ratio$`2023`), 7, subset = NULL)
sir_ratio <- reshape(sir_ratio,
                  timevar = "Year",
                  times = c("2014 - 2016", "2017", "2018", "2019", "2020", "2021", "2022", "2023"),
                  varying = c("2014 - 2016", "2017", "2018", "2019", "2020", "2021", "2022", "2023"),
                  v.names = "SIR",
                  idvar = "province",
                  direction = "long")
map_sir <- merge(nigeria, sir_ratio, by.x = "ADM1_EN", by.y = "province")

breaks <- c(0, 0.25,  0.5, 0.75,  1,  2.5,  6,  9, 14)
num_colors <- length(breaks) - 1
colors <- c("#007bff", "#7699ff", "#abbaff", "#d7dcff", "#ffdacc", "#ffb59b", "#ff8e6b", "#ff633c", "#ff2600")

# Create a new column in your data frame to categorize 'mu' based on Jenks breaks
map_sir$SIR_category <- cut(map_sir$SIR, breaks = breaks, include.lowest = TRUE, labels = FALSE)

# Map the categories to colors
map_sir$mu_color <- colors[map_sir$SIR_category]

labels <- sapply(1:num_colors, function(i) {
  paste0(breaks[i], "-", breaks[i + 1])
})

# Plotting with ggplot2
SIR <- ggplot(map_sir) + 
  geom_sf(aes(fill = factor(SIR_category))) +
  facet_wrap(~Year, dir = "h", ncol = 4) +
  ggtitle("Outbreak Risk Index of Rabies in Nigerian states 2014 - 2023") + 
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank()
  ) +
  scale_fill_manual(
    values = colors,
    name = "Outbreak Risk Index",
    breaks = seq_along(breaks[-1]), 
    labels = labels
  )


png("SIR.png", width = 900, height = 400)
SIR
dev.off()
