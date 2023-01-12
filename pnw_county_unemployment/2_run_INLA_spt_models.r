## Last revised: 01/11/2023
## Assume 
# options(repos = 'https://cran.seoul.go.kr')
# Provided that the git repository is set the working directory
if (!require(pacman)) {install.packages('pacman')}
# INLA setup
# install.packages("INLA",repos = c(getOption("repos"), INLA = "https://inla.r-inla-download.org/R/testing"), dep=TRUE)
p_load(INLA, tidyverse, dtplyr, sf, spdep, tmap, rmapshaper)

# num.threads x:y format; x for stage 1, y for stage 2
# x*y <= (total # of threads available in the machine)
inla.setOption(num.threads="8:1",
               inla.mode = "experimental",
               pardiso.license = "/mnt/c/Users/sigma/OneDrive/0_IS_Personal/pardiso_license.lic")
inla.pardiso.check()
# setwd("__your_directory__")

source('./pnw_county_unemployment/1_base_functions.r')

mort = readRDS("./Data/county_boundary.rds")
mort = mort |>
    dplyr::filter(grepl('^(41|53)', GEOID)) 
## Spatial Weight Matrix ####
mort %>% dplyr::select(GEOID) %>% 
    mutate(GEOID = as.character(GEOID)) %>% 
    arrange(GEOID) %>% 
    poly2nb(snap = 0.005) -> mort.w
mort.wn = sapply(mort.w, length)
graph.wm <- nb2mat(mort.w, style = 'B')
graph = graph.wm

## Data Load ####
## load
mortune.covardf <- readRDS('./Data/cleaned_covariates_allsex.rds')
mu.cd <- mortune.covardf %>%
  filter(grepl('^(41|53)', GEOID)) %>%
  arrange(disease, year, GEOID) %>%
  # group_by(disease, GEOID) %>%
  # ungroup %>%
  rename(UELAG1 = UE) %>%
  mutate(MEDINCOME = ifelse(year.covar == 1990, (MEDINCOME/76.90)*100, ifelse(year.covar == 2010, (MEDINCOME/125.96)*100, MEDINCOME))) %>% 
  filter(year >= 2001) %>%
  mutate(GEOID_int1 = as.integer(plyr::mapvalues(GEOID, unique(GEOID), seq(1, length(unique(GEOID))))),
         GEOID_int2 = GEOID_int1,
         GEOID_int3 = GEOID_int1,
         GEOID_year_int = seq(1, nrow(.)),
         year1 = year - 2000,
         year2 = year1)


## VIF and correlation check
mu.cd_ment = mu.cd %>%
  filter(grepl('Mental.*', disease) & year >= 2001) #%>%
  #mutate(Mortrate = scale(Mortrate))

vif_form = Mortrate ~ UELAG1 + MEDINCOME + PPOVTBELOW + PRENTER + PUNIVABOVE + PNONWHITE + PRURAL + PHSINGLE + PELDER
car::vif(lm(formula = vif_form, data = mu.cd_ment %>% mutate_if(is.numeric, list(~as.vector(scale(.))))))
cor(mu.cd_ment %>% 
    model.frame(vif_form, .) %>% 
    dplyr::select(-1) %>% 
    dplyr::mutate_all(.funs = list(~ifelse(is.na(.), median(., na.rm = TRUE), .)))) %>%
    round(., 3)


## Formulae for main models ####
# Default prior
hyp.prior <- list(prec.unstruct = list(prior = 'loggamma', param = c(1, 0.00005)),
                  prec.spatial = list(prior = 'loggamma', param = c(1, 0.00005)))
hyp.prior.b <- list(prec = list(prior = 'loggamma', param = c(1, 0.00005)))
hyp.prior.i <- list(prec = list(prior = 'loggamma', param = c(100, 0.05)))


### Reviewed versions of model fitting
form.m10 <- Mortrate ~ 
  f(GEOID_int1, model = "bym", graph = graph, constr = TRUE) +
  f(year1, model = 'ar1', constr = TRUE) +
  UELAG1 + PELDER + PNONWHITE + MEDINCOME + PUNIVABOVE + PPOVTBELOW + PRURAL + PRENTER + PHSINGLE
# Model 2
form.m20 <- Mortrate ~ 
  f(GEOID_int1, model = "bym", graph = graph, constr = TRUE) +
  f(year1, model = 'ar1', constr = TRUE) +
  f(GEOID_int2, UELAG1, model = "besag", graph = graph, constr = TRUE, hyper = hyp.prior.b) +
  UELAG1 + PELDER + PNONWHITE + MEDINCOME + PUNIVABOVE + PPOVTBELOW + PRURAL + PRENTER + PHSINGLE
# Model 3
form.m30 <- Mortrate ~ 
  f(GEOID_int1, model = "bym", graph = graph, constr = TRUE) +
  f(year1, model = 'ar1', constr = TRUE) +
  f(year2, UELAG1, model = 'rw1', constr = TRUE) +
  UELAG1 + PELDER + PNONWHITE + MEDINCOME + PUNIVABOVE + PPOVTBELOW + PRURAL + PRENTER + PHSINGLE
## Model 4
form.m40 <- Mortrate ~ 
  f(GEOID_int1, model = "bym", graph = graph, constr = TRUE, hyper = hyp.prior, scale.model = TRUE) +
  f(year1, model = 'ar1', constr = TRUE) +
  f(GEOID_int2, UELAG1, model = "besag", graph = graph, constr = TRUE, hyper = hyp.prior.b, scale.model = TRUE) +
  f(year2, UELAG1, model = 'rw1', constr = TRUE) +
  UELAG1 + MEDINCOME + PPOVTBELOW + PRENTER + PUNIVABOVE + PRURAL + PNONWHITE + PHSINGLE + PELDER

## Model 5: Spatiotemporal interaction types I, II, III, and IV
#--- Type IV interaction with besag ICAR ---#
form.m5IV <- update(form.m40, ~ . + f(GEOID_int3, UELAG1, model = "besag", graph = graph, group = year3, control.group = list(model = "rw1"),
                          constr = TRUE, hyper = hyp.prior.i, scale.model = TRUE) )
#--- Type I interaction ---#
form.m5I <- update(form.m40, ~ . + f(GEOID_year_int, UELAG1, model = "iid", constr = TRUE))
#--- Type II interaction ---#
form.m5II <- update(form.m40, ~ . + f(GEOID_int3, UELAG1, model = "iid", graph = graph, group = year3, 
                                        control.group = list(model = "rw1"),
                                        constr = TRUE) )
#--- Type III interaction ---#
## added following generic0 trick
r.def <- 14
A.constr <- kronecker(diag(14),matrix(1,1,75))
A.constr <- A.constr[-1,]
Q.xi = Diagonal(x = mort.wn) - Matrix(mort.wm)
Q.Leroux <- diag(75)-Q.xi


#########################################################
##  Define the temporal structure matrix of a RW1/RW2  ##
#########################################################
D1 <- diff(diag(14),differences=1)
Q.gammaRW1 <- t(D1)%*%D1
R <- kronecker(Diagonal(14),Q.xi)


form.m5III <- update(form.m40, ~ . + f(GEOID_year_int, UELAG1, 
                                      model = "generic0", 
                                      Cmatrix = R, rankdef = r.def,
                                      constr = TRUE, extraconstr = list(A = A.constr, e = rep(0, 13)) ))

# form.m5III <- update(form.m40, ~ . + f(year3, UELAG1, model = "iid", group = GEOID2, 
#                                          control.group = list(graph = graph, model = 'besag', scale.model = TRUE),
#                                       constr = TRUE) )

# Linear combination for type IV
lcIV <- inla.make.lincombs(year2 = kronecker(Diagonal(14), Matrix(1, 75, 1)),
                           GEOID_int2 = kronecker(Matrix(1, 14, 1), Diagonal(75)), 
                           GEOID_int3 = kronecker(Diagonal(14), Diagonal(75)))


# Theta initialization all to zeroes
thetas10 = rep(0, 5)
thetas20 = rep(0, 6)
thetas30 = rep(0, 6)
thetas40 = rep(0, 7)
thetas50 = rep(0, 8)

# Run model fitting
core_config = "8:1"
system.time(inla_mort_mental_10 <- fit_inla_spt(form.m10, mu.cd, dis = 'Mentalandsubstanceusedisorders', 
                                                verbose = TRUE, theta = thetas10, cores = core_config, logarithm = TRUE,
                                                standardize = TRUE))
system.time(inla_mort_mental_20 <- fit_inla_spt(form.m20, mu.cd, dis = 'Mentalandsubstanceusedisorders', 
                                                verbose = TRUE, theta = thetas20, cores = core_config, logarithm = TRUE,
                                                standardize = TRUE))
system.time(inla_mort_mental_30 <- fit_inla_spt(form.m30, mu.cd, dis = 'Mentalandsubstanceusedisorders', 
                                                verbose = TRUE, theta = thetas30, cores = core_config, logarithm = TRUE,
                                                standardize = TRUE))
system.time(inla_mort_mental_40 <- fit_inla_spt(form.m40, mu.cd, dis = 'Mentalandsubstanceusedisorders', 
                                                standardize = TRUE,
                                                cores = core_config, logarithm = TRUE,
                                                verbose = TRUE, theta = thetas40))
system.time(inla_mort_mental_51 <- fit_inla_spt(form.m5I, mu.cd, dis = 'Mentalandsubstanceusedisorders', 
                                                standardize = TRUE,
                                                cores = core_config, logarithm = TRUE,
                                                verbose = TRUE, theta = thetas50))
system.time(inla_mort_mental_52 <- fit_inla_spt(form.m5II, mu.cd, dis = 'Mentalandsubstanceusedisorders', 
                                                cores = core_config, logarithm = TRUE,
                                                verbose = TRUE, theta = thetas50,
                                                standardize = TRUE))
system.time(inla_mort_mental_53 <- fit_inla_spt(form.m5III, mu.cd, dis = 'Mentalandsubstanceusedisorders', 
                                                cores = core_config, logarithm = TRUE,
                                                verbose = TRUE, theta = thetas50,
                                                standardize = TRUE))
# Final IV model
system.time(inla_mort_mental_54 <- fit_inla_spt(form.m5IV, mu.cd, dis = 'Mentalandsubstanceusedisorders', 
                                                lincomb = lcIV, 
                                                theta = thetas50, verbose = TRUE, 
                                                logarithm = TRUE, cores = core_config, 
                                                standardize = TRUE))


# Repeat model fitting until the sound mode status is obtained
j <- 1
while (inla_mort_mental_10$mode$mode.status %% 1000 != 0){
  inla_mort_mental_10 <- inla.rerun(inla_mort_mental_10)
  j <- j + 1
}
j <- 1
while (inla_mort_mental_20$mode$mode.status %% 1000 != 0){
  inla_mort_mental_20 <- inla.rerun(inla_mort_mental_20)
  j <- j + 1
}
j <- 1
while (inla_mort_mental_30$mode$mode.status %% 1000 != 0){
  inla_mort_mental_30 <- inla.rerun(inla_mort_mental_30)
  j <- j + 1
}
j <- 1
while (inla_mort_mental_40$mode$mode.status %% 1000 != 0){
  inla_mort_mental_40 <- inla.rerun(inla_mort_mental_40)
  j <- j + 1
}

j <- 1
while (inla_mort_mental_51$mode$mode.status %% 1000 != 0){
  inla_mort_mental_51 <- inla.rerun(inla_mort_mental_51)
  j <- j + 1
}
j <- 1
while (inla_mort_mental_52$mode$mode.status %% 1000 != 0){
  inla_mort_mental_52 <- inla.rerun(inla_mort_mental_52)
  j <- j + 1
}
j <- 1
while (inla_mort_mental_53$mode$mode.status %% 1000 != 0){
  inla_mort_mental_53 <- inla.rerun(inla_mort_mental_53)
  j <- j + 1
}
j <- 1
while (inla_mort_mental_54$mode$mode.status %% 1000 != 0){
  inla_mort_mental_54 <- inla.rerun(inla_mort_mental_54)
  j <- j + 1
}


## WAIC ####
inla_mort_mental_10$waic$waic
inla_mort_mental_20$waic$waic
inla_mort_mental_30$waic$waic
inla_mort_mental_40$waic$waic
inla_mort_mental_51$waic$waic
inla_mort_mental_52$waic$waic
inla_mort_mental_53$waic$waic
inla_mort_mental_54$waic$waic

## DIC ####
inla_mort_mental_10$dic$dic
inla_mort_mental_20$dic$dic
inla_mort_mental_30$dic$dic
inla_mort_mental_40$dic$dic
inla_mort_mental_51$dic$dic
inla_mort_mental_52$dic$dic
inla_mort_mental_53$dic$dic
inla_mort_mental_54$dic$dic
