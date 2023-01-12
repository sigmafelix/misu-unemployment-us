## Miscellaneous plots

## Check normality (visually)
gg_normality = 
    mu.cd %>%
    filter(disease == 'Mentalandsubstanceusedisorders') %>%
    ggplot(data = .,
           mapping = aes(x = Mortrate)) +
        facet_wrap(~year) +
        geom_density()
gg_normality

gg_covars = 
    mu.cd %>%
    filter(disease == 'Mentalandsubstanceusedisorders') %>%
    dplyr::select(disease, year, Mortrate, UELAG1, PELDER, PNONWHITE, PUNIVABOVE, PPOVTBELOW, PRURAL, PRENTER, PHSINGLE) %>%
    pivot_longer(cols = 3:ncol(.)) %>%
    ggplot(data = .,
           mapping = aes(x = value)) +
        facet_wrap(~name) +
        geom_density()



### Miscellaneous figures
state = tigris::states(year = 2010)
st_crs(mort) <- 4269
mort.simp <- mort %>% 
  filter(GEOID %in% mu.cd$GEOID) %>% 
  dplyr::select(GEOID) %>% 
  arrange(GEOID) %>% 
  st_transform(crs = 2163) %>% 
  ms_simplify(keep = 0.1, keep_shapes = TRUE)
state.simp = state %>%
    filter(grepl('^(41|53)', GEOID10)) %>%
    ms_simplify(keep = 0.33, keep_shapes = TRUE) %>%
    st_transform(crs = 2163)


mean_54_min = min(inla_mort_mental_54$summary.lincomb.derived$mean)
mean_54_max = max(inla_mort_mental_54$summary.lincomb.derived$mean)
brks = c(mean_54_min, -0.2, -0.1, 0, 0.1, 0.2, mean_54_max)

if (!dir.exists("./output")) {
  dir.create("./output")
}

map_dexport(inla_mort_mental_54, rand.d = NULL, map_state = state.simp,
            pdir = './output/',
            #brks = brks,
            filename = 'Mortality_Mental_54_log_rwbesag_beta1it_011123.png', file.export = TRUE,
            title = expression(beta['1it']))


map_dexport(inla_mort_mental_54, rand.d = NULL, map_state = state.simp,
            pdir = './output/',
            filename = 'Mortality_Mental_54_log_beta1it_120421.png', file.export = TRUE,
            title = expression(beta['1it']), additional = TRUE)


map_dexport(inla_mort_mental_54, rand.d = NULL, map_state = state.simp,
            pdir = './output/',
            wd = 1280, hg = 800,
            filename = 'Mortality_Mental_54_beta1it_111121.gif', file.export = TRUE,
            title = expression(beta['1it']), animation = TRUE)



## Residual map
2001:2014 %>% 
  split(.,.) %>% 
  lapply(function(x) map_inla_rs(inla_mort_mental_54, mort.simp, year_select = x, rand.d = 5)) -> map_rs_all_l
map_rs_all_a <- do.call(rbind, map_rs_all_l) %>% 
      dplyr::select(GEOID, year, mean)


inla_mort_mental_54$.args$data %>%
  mutate(GEOID = mort.simp$GEOID[GEOID],
         year = year + 2000,
         resid_mental = Mortrate - inla_mort_mental_54$summary.fitted.values$mean) %>%
    dplyr::select(GEOID, year, resid_mental) -> inla_mort_mental_54r
colnames(inla_mort_mental_54r)[3] <- 'resid_mental'
map_rs_all_a %>%
    left_join(inla_mort_mental_54r) -> map_rs_all_aa

tms_a <- tm_shape(map_rs_all_aa %>% filter(!is.na(resid_mental))) +
    tm_facets(by = 'year') +
    tm_borders(col = 'light grey', lwd = 0.08) +
    tm_fill('resid_mental', palette = '-RdBu', 
            title = 'Residuals') +
    tm_layout(frame = FALSE, 
              frame.lwd = NA, 
              panel.label.bg.color = NA, 
              panel.label.size = 1.5,
              legend.outside = TRUE, 
              legend.outside.position = 'right',
              outer.margins = c(0.01,0.01,0.01,-0.1)) +
    tm_shape(state.simp) +
    tm_borders(col = 'black', lwd = 0.2)

tmap_save(tms_a, 
            filename = './output/Residuals_mental_log_54.png', 
            width = 50, height = 50,
            units = 'cm', dpi = 300, pointsize = 30)


### Lagged unemployment rate
mu.cdm = mu.cd %>%
  filter(grepl('^Mental', disease) & year >= 2001)
map_rs_all_lue = map_rs_all_a %>%
  left_join(mu.cdm, by = c('GEOID', 'year'))
tms_lue <- tm_shape(map_rs_all_lue %>% filter(!is.na(UELAG1))) +
    tm_facets(by = 'year') +
    tm_borders(col = 'light grey', lwd = 0.08) +
    tm_fill('UELAG1', palette = 'Reds', style = 'fisher', #breaks=c(0,0.2,0.5,0.8,1),
            title = 'One-year lagged \nUnemployment rate (%)') +
    tm_layout(frame = FALSE, frame.lwd = NA, panel.label.bg.color = NA, 
              panel.label.size = 1.5,
              legend.outside = TRUE, legend.outside.position = 'right',
              outer.margins = c(0.01,0.01,0.01,-0.08)) +
    tm_shape(state.simp) +
    tm_borders(col = 'black', lwd = 0.5)
tmap_save(tms_lue, 
            filename = './output/UELAG1_maps.png', 
            width = 50, height = 50,
            units = 'cm', dpi = 300, pointsize = 30)

### Mortality rate
tms_mtr <- tm_shape(map_rs_all_lue %>% filter(!is.na(Mortrate))) +
    tm_facets(by = 'year') +
    tm_borders(col = 'light grey', lwd = 0.08) +
    tm_fill('Mortrate', palette = 'Reds', style = 'fisher', 
            title = 'Mental illness and \nsubstance-use \nmortality rate (per 100,000)') +
    tm_layout(frame = FALSE, frame.lwd = NA, panel.label.bg.color = NA, 
              panel.label.size = 1.5,
              legend.outside = TRUE, legend.outside.position = 'right',
              outer.margins = c(0.01,0.01,0.01,-0.05)) +
    tm_shape(state.simp) +
    tm_borders(col = 'black', lwd = 0.5)
tmap_save(tms_mtr, filename = './output/Mortality_rates_maps.png', 
            width = 50, height = 50,
            units = 'cm', dpi = 300, pointsize = 30)






### Review: KStest
mu.cd_ment = mu.cd %>%
  filter(grepl('Mental.*', disease) & year >= 2001)
ks.test(x = mu.cd_ment$Mortrate, 'pnorm', mean(mu.cd_ment$Mortrate), sd(mu.cd_ment$Mortrate))$p.value
shapiro.test(mu.cd_ment$Mortrate)

kstest_simple = function(vec) {
  ks.test(x = vec, 'pnorm', mean(vec), sd(vec))$p.value
}

ks.test(mu.cd_ment$Mortrate, 'pnorm', mean(mu.cd_ment$Mortrate), sd(mu.cd_ment$Mortrate)) %>% str
tseries::jarque.bera.test(scale(mu.cd_ment$Mortrate))
e1071::skewness(mu.cd_ment$Mortrate %>% log)
hist(mu.cd_ment$Mortrate)

mu.cd_ment %>%
  mutate(Mortrate = scale(Mortrate)) %>%
  group_by(year) %>%
  nest() %>%
  mutate(ks_pval = map_dbl(data, ~kstest_simple(.x$Mortrate))) %>% .$ks_pval
  ungroup


mortrate_gg = ggplot(data = mu.cd_ment, mapping = aes(x = (Mortrate))) +
  facet_wrap(~year) +
  geom_histogram() +
  theme_bw()

