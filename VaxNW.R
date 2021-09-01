library(EpiModel)

set.seed(1)
nw <- network.initialize(n = 1000, directed = FALSE)
age.grp <- sample(1:3, 1000, replace = TRUE)
nw <- set.vertex.attribute(nw, "age.grp", age.grp)

#formation takes into account average total num of edges at each time step, num of edges with same group nodes
formation <- ~edges + nodematch("age.grp", diff = TRUE)

#target stats correspond to formation factors
#Suppose we expect young to have the most matched edges followed by middle, then old. 
target.stats <- c(500, 150, 100, 50)

#Adding departure (death) rate
departure_rate <- 0.008

#Suppose duration of non-matched edges is 10 time steps, matched for young is 50, middle is 40, old is 20
durs <- c(10, 50, 40, 20)
#Dissolution coefficient specifying avg num of time steps for an edge's existence accounting for different age groups
diss <- dissolution_coefs(~offset(edges) +
                            offset(nodematch("age.grp", diff = TRUE)),
                          duration = durs)

#Estimating the network with 3 groups and specified target statistics
est1 <- netest(nw, formation, target.stats, diss, edapprox = TRUE)


#Simulating SEIRS model over network

#Parameter spec
#departure rate mult represents the multiplier for departure rate which indicates probability of death for infected nodes

param <- param.net(inf.prob = 0.6,
                   arrival.rate = 0.01,
                   inf.prob.age.grp2 = 0.4,
                   inf.prob.age.grp3 = 0.3,
                   departure.rate = departure_rate,
                   departure.disease.mult = 2,
                   act.rate = 1,
                   ei.rate = 0.05,
                   ir.rate = 0.05,
                   rs.rate = 0.05,
                   vaccination.rate.initialization = 0.05,
                   protection.rate.initialization = 0.8,
                   vaccination.rate.progression = 0.05,
                   protection.rate.progression = 0.8,
                   vaccination.rate.arrivals = 0.3,
                   protection.rate.arrivals = 0.8,
                   vaccine.efficacy = 0.8)


init <- init.net(i.num = 50)

#Vax model
source("module-fx.R")

#To split output by age group - epi.by = "age.grp"
control <- control.net(type = NULL, nsteps = 200, nsims = 3, epi.by = "age.grp",  infection.FUN = infect,
                       progress.FUN = progress,
                       departures.FUN = dfunc,
                       arrivals.FUN = afunc,
                       resimulate.network = TRUE,
                       verbose = TRUE,
                       module.order = c("resim_nets.FUN", "departures.FUN",
                                       "arrivals.FUN", "infection.FUN",
                                        "progress.FUN", "prevalence.FUN"))

#Final Sim
set.seed(1)
sim <- netsim(est1, param, init, control)

#Dataframe with resulting stats
df <- as.data.frame(sim)

#Network Plot which doesn't produce correct cross-section
par(mfrow = c(1,2), mar = c(0,0,1,0))
plot(sim, type = "network", at = 200, col.status = TRUE, sims = 3,
     main = "Prevalence at t200")

