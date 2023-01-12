## 12/26/2022

if (!require(pacman)) { install.packages('pacman'); library(pacman)}
if (!require(INLA)) { 
    install.packages('INLA', 
                     repos = c(getOption('repos'), "INLA"="https://inla.r-inla-download.org/R/testing"),
                     dependencies = TRUE)
}
p_load(tidyverse, spdep, tmap, rmapshaper,sf, dtplyr, data.table, CCA, FactoMineR, car)


inla.setOption(
        pardiso.license = '/home/felix/OneDrive/0_IS_Personal/pardiso_license.lic',
        mkl = FALSE,
        inla.mode = 'experimental')

inla.pardiso.check()

# Read file
orwa_cnty_0618 = read_rds("./pnw_reanalysis/ORWA_county_0618_covariates_122622.rds") %>%
    mutate(GEOID_int1 = as.integer(plyr::mapvalues(GEOID10, unique(GEOID10), seq(1, length(unique(GEOID10))))),
           GEOID_int2 = GEOID_int1,
           GEOID_int3 = GEOID_int2,
           year_int1 = year - 2005,
           year_int2 = year_int1,
           GEOID_year1 = seq(1, nrow(.))) %>%
    mutate(d_parks_nonpriv_10k = 1e4 * n_parks_nonpriv / n_pop_total)



### Spatial weight matrix (first-order) ####
scalemodel = TRUE
swm_orwa_q1 = 
    orwa_cnty_0618 %>%
    filter(year == 2010) %>%
    poly2nb %>%
    nb2mat(style = 'B')
graph = swm_orwa_q1


form_base_spt <- outcome ~ 
    offset(log(n_pop_total)) +
    f(GEOID_int1, model = 'bym2', graph = graph, constr = TRUE, scale.model = scalemodel) +
    f(year_int1, model = 'ar1', constr = TRUE) +
    State +
    green + 
    pm25_mean +
    log(p_water_pct_strict+ exp(-25)) +
    d_nursing_10k +
    SES + 
    p_elderly + p_singleheaded + p_bmarital + log(p_nonwhite + 0.1)

source("./pnw_reanalysis/1_functions.r")

form_base_spt1 = update(form_base_spt,
            ~.+ f(GEOID_year1, model = 'iid', constr = TRUE)
            )
form_base_spt2 = update(form_base_spt,
            ~.+ f(GEOID_int2, model = 'bym2', group = year_int2, control.group = list(model = 'ar1'), graph = graph, constr = TRUE)
            )
form_base_spt3 = update(form_base_spt,
            ~.+ f(GEOID_int2, model = 'bym2', group = year_int2, control.group = list(model = 'iid'), graph = graph, constr = TRUE)
            )
form_base_spt4 = update(form_base_spt,
            ~.+ f(GEOID_int2, model = 'bym2', group = year_int2, control.group = list(model = 'ar1'), graph = graph, constr = TRUE)
            )

hyper3 = list(prec = list(theta1 = 'loggamma', param = c(1, 0.01)))

year_length = length(unique(orwa_cnty_0618$year))
tract_length = length(unique(orwa_cnty_0618$GEOID10))

# Linear combination definition for the combined effect
# BYM2 random effect (structure + non-structure)
lcs_bym1 = inla.make.lincombs(GEOID_int3 = cbind(Diagonal(tract_length), Diagonal(tract_length)))#,
# BYM2 fixed + spatial random (structure + non-structure)
lcs_bym2 = inla.make.lincombs(green = Matrix(1, nrow = tract_length, ncol = 1),
                         GEOID_int3 = cbind(Diagonal(tract_length), Diagonal(tract_length)))#,

# spatial random slope
form_nsub_4ext = update(form_base_spt4,
                    ~ . + f(GEOID_int3, green, model = 'bym2', graph = graph, constr=TRUE, scale.model = scalemodel))
form_nsub_3ext = update(form_base_spt3,
                    ~ . + f(GEOID_int3, green, model = 'bym2', graph = graph, constr=TRUE, scale.model = scalemodel))
form_nsub_2ext = update(form_base_spt2,
                    ~ . + f(GEOID_int3, green, model = 'bym2', graph = graph, constr=TRUE, scale.model = scalemodel))
form_nsub_1ext = update(form_base_spt1,
                    ~ . + f(GEOID_int3, green, model = 'bym2', graph = graph, constr=TRUE, scale.model = scalemodel))

# determine the thread composition in your own discretion
threads = "8:1"

## Main model fitting: takes some time (3+ hours)
# Spatial only random slope models (Model 5)
mod_nsub_aq1_4pois = fit_model(form_nsub_4ext, orwa_cnty_0618, lincomb=lcs_bym1, outcome = 'n_nonsubstance2', family = 'poisson', graph = graph, greenvar = "d_parks_nonpriv_10k", threads = threads, thetas = NULL)
mod_nsub_aq1_4pois = repeat_refit(mod_nsub_aq1_4pois)

# Models 1-4
mod_nsub_aq1_0pois = fit_model(form_base_spt, orwa_cnty_0618, lincomb = NULL, outcome = 'n_nonsubstance2', family = 'poisson', graph = graph, greenvar = "d_parks_nonpriv_10k", threads = threads, thetas = NULL)
mod_nsub_aq1_1pois = fit_model(form_nsub_1ext, orwa_cnty_0618, lincomb=NULL, outcome = 'n_nonsubstance2', family = 'poisson', graph = graph, greenvar = "d_parks_nonpriv_10k", threads = threads, thetas = NULL)
mod_nsub_aq1_2pois = fit_model(form_nsub_2ext, orwa_cnty_0618, lincomb=NULL, outcome = 'n_nonsubstance2', family = 'poisson', graph = graph, greenvar = "d_parks_nonpriv_10k", threads = threads, thetas = NULL)
mod_nsub_aq1_3pois = fit_model(form_nsub_3ext, orwa_cnty_0618, lincomb=NULL, outcome = 'n_nonsubstance2', family = 'poisson', graph = graph, greenvar = "d_parks_nonpriv_10k", threads = threads, thetas = NULL)

mod_nsub_aq1_0pois = repeat_refit(mod_nsub_aq1_0pois)
mod_nsub_aq1_1pois = repeat_refit(mod_nsub_aq1_1pois)
mod_nsub_aq1_2pois = repeat_refit(mod_nsub_aq1_2pois)
mod_nsub_aq1_3pois = repeat_refit(mod_nsub_aq1_3pois)




# Export plots ####
vis_spt_effect(map = orwa_cnty_0618,
                inla_fit = mod_nsub_aq1_4pois, effect_id = NULL, file_export = TRUE,
                filepath = "fig_beta1i.png")
map_emarginal(mod_nsub_aq1_4pois, map = orwa_cnty_0618 %>% filter(year == 2010), rand.d = NULL,
              file.export = TRUE, pdir = "./",
              filename = "fig_pbeta1i_lt0.png")

## Collinearity test ####
form.subs <- n_substance ~ offset(log(tpop)) +
    urban + n_medincome_100k + EVI_1y + pm25_mean +
    p_nonwhite + p_youngadult + p_poverty + p_edubac + p_unemp + 
    p_renter + p_vacant + p_singleheaded + p_noncitizen + p_veteran + p_bmarital + p_fulltime + p_substancefac_100k

car::vif(glm(form.subs, poisson, data = orwa_tracts_0618_mdpp))
car::vif(glm(form.nsub, poisson, data = orwa_tracts_0618_mdpp))
