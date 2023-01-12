### Data cleaning for PNW reanalysis
### 12/25/2022 21:21 PST
### Insang Song
if (!require(pacman)){
  install.packages('pacman')
  library(pacman)
}
p_load(tmap, tmaptools, tidyverse, dtplyr, sf, spdep, rvest)
p_load(tidycensus, tigris, FactoMineR, exactextractr, terra)

sf_use_s2(FALSE)
# For paper 2 main data: ASMR ####
username = 'felix'
# Detects /mnt/c/ then assigns appropriate drive value
# currently supports WSL vs Linux
if (dir.exists("/mnt/c/")) {
  drive = sprintf('/mnt/c/Users/%s/OneDrive/', username)
} else {
  drive = sprintf('/home/%s/OneDrive/', username)
}
datapath = str_c(drive, "UO/Dissertation/")
# source(str_c(datapath, 'Code/Common/Spatial_OR_cleaning_012922.r'))
load(str_c(datapath, 'Data/Tracts_ACS_ORWA_0919_121621.RData'))
load(str_c(datapath, 'Data/ACS_Population_ORWA_012722.RData'))
# modify path_to to the proper directory path


# CPI ####
cpi = read_csv(str_c(datapath, "/Data/CPI_1913_2021.csv"))
cpi_0619 = cpi %>%
  mutate(Year = 1913:2021) %>%
  filter(Year >= 2006 & Year <= 2019) %>%
  mutate(CPI10 = `Annual Average`[5] / `Annual Average`) %>%
  dplyr::select(1,3)


## Data acquisition: ACS 5-yr ####
capi = "05db0646f9956eba295b88d5a86bf4fe455459f3"
tidycensus::census_api_key(capi)

cnty_or = tigris::counties(state = 'Oregon', year = 2010)
cnty_wa = tigris::counties(state = 'Washington', year = 2010)

cnty_orwa = rbind(cnty_or, cnty_wa) %>%
  st_transform(5070) %>%
  arrange(GEOID10)

zcta_or = tigris::zctas(state = 'Oregon', year = 2010)
zcta_wa = tigris::zctas(state = 'Washington', year = 2010)

zcta_orwa = rbind(zcta_or, zcta_wa) %>%
  st_transform(5070) %>%
  arrange(ZCTA5CE10)


lvars <- load_variables(2009, 'acs5')
lvars13 <- load_variables(2013, 'acs5')
# lvars12 <- load_variables(2012, 'acs5')
# lvars11 <- load_variables(2011, 'acs5')
# lvars17 <- load_variables(2017, 'acs5')
lvars10 <- load_variables(2010, 'acs5')
acs_income <- lvars %>% filter(grepl('*(H|h)ousehold*', label))
acs_poverty <- lvars %>% filter(grepl('*(P|p)overty*', label))
acs_unins <- lvars13 %>% filter(grepl('*(U|u)ninsur|insur*', label))

acs_variables <- lvars %>% filter(grepl('B01001(|[A-I])_[0-9]{3}', name)) %>% data.frame


codefilter <- '^(B01001|B17001|B19013|B19083|B23025|B25003|B25004|B15001|B15003|B27001|B05001|B08201)(|[A-Z])_.*'
acs_variables <- lvars13 %>% filter(grepl(codefilter, name)) %>% data.frame
# B01001[A-I]_001-099: total population (sex, age, ethnicity)
# B01001C_001-099: total population (AIAN)
# B19013 Median household income (past 12 months)
# B19083 Gini Index
# B15003 Educational attainment (25+ yo)
# B23001 Employment status (16+ yo)
# B25001 Total housing units
# B25003 Tenure
# B25004 Vacancy Status
# B25009 Tenure by household size
# B15001 Sex by age by Educational Attainment for the population 18+yo
# B17001 Poverty status in the past 12 months sex by age
# B27001 Insurance status by sex by age
# B08201 Household size by vehicles available
# B05001 Nativity and citizenship status
# B12002 Sex by marital status
# B06008 Place of birth by marital status//Number of times married by sex by marital status
# B20005 Sex by work experience in the past 12 months (full-time)
# B21001 Veteran status

# Minimal code lists
b23001_names = lvars10 %>%
    filter(grepl('B23001', name)) %>%
    .$name
b23001_label = lvars10 %>%
    filter(grepl('B23001', name)) %>%
    .$label
b23001_names_inlabor = b23001_names[grep('In labor force$', b23001_label)]
b23001_names_unemp = b23001_names[grep('Unemployed$', b23001_label)]


code_list_10 = c('B19013_001',
              # 'B19083_001',
              'B01001_020', 'B01001_021', 'B01001_022', 'B01001_023', 'B01001_024', 'B01001_025', 'B01001_044', 'B01001_045', 'B01001_046', 'B01001_047', 'B01001_048', 'B01001_049', 'B01001_001',
              'B01001_007', 'B01001_008', 'B01001_009', 'B01001_010', 'B01001_031', 'B01001_032', 'B01001_033', 'B01001_034', 
              'B01001A_001', # white alone
              'B17001_002', 'B17001_001',
              'B15002_001', 'B15002_015', 'B15002_016', 'B15002_017', 'B15002_018', 'B15002_032', 'B15002_033', 'B15002_034', 'B15002_035',
              #'B06009_001', 'B06009_005', 'B06009_006',
              'B23001_001',
              b23001_names_inlabor,
              b23001_names_unemp,
              'B25003_003', 'B25003_001',
              'B25009_001', 'B25009_003', 'B25009_011',
              #'B08201_001', 'B08201_007',
              'B21001_001', 'B21001_002',
              'B05001_001', 'B05001_006',
              'B20005_001', 'B20005_003', 'B20005_050',
              'B12002_001', 'B12002_035', 'B12002_050', 'B12002_065', 'B12002_080', 'B12002_129', 'B12002_143', 'B12002_158', 'B12002_173',
              #'B06008_001', 'B06008_004', 'B06008_005', 'B06008_006',
              'B25001_001',
              'B25004_003', 'B25004_004', 'B25004_005', 'B25004_001'#,
              #'B27001_005', 'B27001_008', 'B27001_011', 'B27001_014', 'B27001_017', 'B27001_020', 'B27001_023', 'B27001_026', 'B27001_029', 'B27001_033', 'B27001_036', 'B27001_039', 'B27001_042', 'B27001_045', 'B27001_048', 'B27001_051', 'B27001_054', 'B27001_057', 'B27001_001'
              )


code_list = c('B19013_001',
              # 'B19083_001',
              'B01001_020', 'B01001_021', 'B01001_022', 'B01001_023', 'B01001_024', 'B01001_025', 'B01001_044', 'B01001_045', 'B01001_046', 'B01001_047', 'B01001_048', 'B01001_049', 'B01001_001',
              'B01001_007', 'B01001_008', 'B01001_009', 'B01001_010', 'B01001_031', 'B01001_032', 'B01001_033', 'B01001_034', 
              'B01001A_001', # white alone
              'B17001_002', 'B17001_001',
              #'B15003_002', 'B15003_003', 'B15003_004', 'B15003_005', 'B15003_006', 'B15003_007', 'B15003_008', 'B15003_009', 'B15003_010', 'B15003_011', 'B15003_012', 'B15003_013', 'B15003_001',
              #'B15003_022', 'B15003_023', 'B15003_024', 'B15003_025', 'B15003_001',
              'B15002_001', 'B15002_015', 'B15002_016', 'B15002_017', 'B15002_018', 'B15002_032', 'B15002_033', 'B15002_034', 'B15002_035',
              #'B06009_001', 'B06009_005', 'B06009_006',
              'B23025_005', 'B23025_002',
              'B25003_003', 'B25003_001',
              'B08201_001', 'B08201_007',
              'B21001_001', 'B21001_002',
              'B05001_001', 'B05001_006',
              'B20005_001', 'B20005_003', 'B20005_050',
              'B12505_001', 'B12505_003', 'B12505_009',
              'B25001_001',
              'B25004_003', 'B25004_004', 'B25004_005', 'B25004_001',
              'B27001_005', 'B27001_008', 'B27001_011', 'B27001_014', 'B27001_017', 'B27001_020', 'B27001_023', 'B27001_026', 'B27001_029', 'B27001_033', 'B27001_036', 'B27001_039', 'B27001_042', 'B27001_045', 'B27001_048', 'B27001_051', 'B27001_054', 'B27001_057', 'B27001_001'
              )

code_list = c('B19013_001',
              # 'B19083_001',
              'B01001_020', 'B01001_021', 'B01001_022', 'B01001_023', 'B01001_024', 'B01001_025', 'B01001_044', 'B01001_045', 'B01001_046', 'B01001_047', 'B01001_048', 'B01001_049', 'B01001_001',
              'B01001_007', 'B01001_008', 'B01001_009', 'B01001_010', 'B01001_031', 'B01001_032', 'B01001_033', 'B01001_034', 
              'B01001A_001', # white alone
              'B17001_002', 'B17001_001',
              'B15002_001', 'B15002_015', 'B15002_016', 'B15002_017', 'B15002_018', 'B15002_032', 'B15002_033', 'B15002_034', 'B15002_035',
              #'B06009_001', 'B06009_005', 'B06009_006',
              'B23001_001',
              b23001_names_inlabor,
              b23001_names_unemp,
              'B25003_003', 'B25003_001',
              'B25009_001', 'B25009_003', 'B25009_011',
              #'B08201_001', 'B08201_007',
              'B21001_001', 'B21001_002',
              'B05001_001', 'B05001_006',
              'B20005_001', 'B20005_003', 'B20005_050',
              'B12002_001', 'B12002_035', 'B12002_050', 'B12002_065', 'B12002_080', 'B12002_129', 'B12002_143', 'B12002_158', 'B12002_173',
              #'B06008_001', 'B06008_004', 'B06008_005', 'B06008_006',
              'B25001_001',
              'B25004_003', 'B25004_004', 'B25004_005', 'B25004_001',
              'B27001_005', 'B27001_008', 'B27001_011', 'B27001_014', 'B27001_017', 'B27001_020', 'B27001_023', 'B27001_026', 'B27001_029', 'B27001_033', 'B27001_036', 'B27001_039', 'B27001_042', 'B27001_045', 'B27001_048', 'B27001_051', 'B27001_054', 'B27001_057', 'B27001_001'
              )



code_list_new = 
    c('B01001_020', 'B01001_021', 'B01001_022', 'B01001_023', 'B01001_024', 'B01001_025', 'B01001_044', 'B01001_045', 'B01001_046', 'B01001_047', 'B01001_048', 'B01001_049', 'B01001_001',
      'B01001_007', 'B01001_008', 'B01001_009', 'B01001_010', 'B01001_031', 'B01001_032', 'B01001_033', 'B01001_034'
        )

get_clear_acs_new = function(geography = 'county', year, state, survey = 'acs5', variable_list = acs_variables$name) {
    acs_5yr <- get_acs(geography = geography, 
                        year = year, 
                        survey = survey, 
                        variables = variable_list, 
                        state = state)
    if (any(nchar(acs_5yr$variable) > 10)) {
        acs_5yr = acs_5yr %>%
            mutate(variable = str_extract(variable, '[A-Z][0-9]{5}_[0-9]{3}'))
    }
    acs_5yr_c <- tryCatch(expr = {
        acs_5yr %>% 
        dplyr::select(-moe) %>% 
        pivot_wider(id_cols = c(GEOID, NAME), names_from = variable, values_from = c(estimate)) %>%
        transmute(GEOID = GEOID,
              p_elderly = 100 * (B01001_020+B01001_021+B01001_022+B01001_023+B01001_024+B01001_025+B01001_044+B01001_045+B01001_046+B01001_047+B01001_048+B01001_049)/B01001_001,
              p_youngadult = 100 * (B01001_007+B01001_008+B01001_009+B01001_010+B01001_031+B01001_032+B01001_033+B01001_034)/B01001_001
                )},
            error = function(e) {print(e) 
                                 return(acs_5yr)}
    )

    acs_5yr_moe <- acs_5yr %>% 
        dplyr::select(-estimate) %>% 
        pivot_wider(id_cols = c(GEOID, NAME), names_from = variable, values_from = c(moe)) %>%
        rename_at(.vars = vars(-1:-2), .funs = list(~str_c(., '_moe')))

    return(list(est = acs_5yr_c, moe = acs_5yr_moe))

}

years = 2009:2019 %>% split(.,.)


get_allyrs_acs_new = function(years, state, variable_list_old, variable_list) {
    yrlist = years %>% split(.,.)
    acs_yrs = yrlist %>%
        lapply(function(x) {
                get_clear_acs_new(year = x, state = state, variable_list = if (x <= 2011) variable_list_old else variable_list)
        })
    return(acs_yrs)
}

or_tract_0919 = get_allyrs_acs_new(2009:2019, 'Oregon', code_list_new, code_list_new)
wa_tract_0919 = get_allyrs_acs_new(2009:2019, 'Washington', code_list_new, code_list_new)

or_tract_0919d = mapply(function(x, y) x$est %>% mutate(year = y),
                        or_tract_0919, 2009:2019 %>% split(.,.), SIMPLIFY = FALSE) %>%
    Reduce(rbind, .)
wa_tract_0919d = mapply(function(x, y) x$est %>% mutate(year = y),
                        wa_tract_0919, 2009:2019 %>% split(.,.), SIMPLIFY = FALSE) %>%
    Reduce(rbind, .)
orwa_tract_0919d = bind_rows(or_tract_0919d, wa_tract_0919d)

cw_orwa = bind_rows(cw_or, cw_wa)
orwa_tract_09d = orwa_tract_0919d %>% filter(year == 2009)
orwa_tract_09df = orwa_tract_09d %>%
    full_join(cw_orwa %>% mutate(GEOID00 = as.character(GEOID00)), by = c('GEOID' = 'GEOID00')) %>%
    group_by(GEOID10) %>%
    summarize_at(.vars = vars(p_elderly, p_youngadult),
                .funs = list(~sum(. * ((POP00 * POPPCT00 / 100)/ sum(POP00 * POPPCT00 / 100))))) %>%
    ungroup %>%
    mutate(GEOID = as.character(GEOID10)) %>%
    dplyr::select(-GEOID10)

orwa_tract_0919df = bind_rows(
    orwa_tract_09df %>% mutate(year = 2006),
    orwa_tract_09df %>% mutate(year = 2007),
    orwa_tract_09df %>% mutate(year = 2008),
    orwa_tract_09df %>% mutate(year = 2009),
    orwa_tract_0919d %>% filter(year >= 2010)
)

## Gini Index ####
ginis = list.files(pattern = "ACSDT5Y[0-9]{,4}.B19083_data*.*.csv",
                   path = "./Data/Gini/",
                   full.names = TRUE) %>%
    split(.,.) %>%
    mapply(function(x, y) {
            xread = read_csv(x, skip = 2, col_names = FALSE) 
            colnames(xread) = c('i_gini', 'i_gini_margin', 'GEO_ID', 'name')
            xread %>% transmute(GEOID = str_sub(GEO_ID, 10, 20), year = y, i_gini = as.numeric(i_gini))},
            ., 2010:2018 %>% split(.,.), SIMPLIFY = FALSE)
ginis_df = Reduce(rbind, ginis)

orwa_tract_0919_add = orwa_tract_0919df %>%
    left_join(ginis_df)


# Original
get_clear_acs = function(geography = 'county', year, state, survey = 'acs5', variable_list = acs_variables$name) {
    acs_5yr <- get_acs(geography = geography, 
                        year = year, 
                        survey = survey, 
                        variables = variable_list, 
                        state = state)
    if (any(nchar(acs_5yr$variable) > 10)) {
        acs_5yr = acs_5yr %>%
            mutate(variable = str_extract(variable, '[A-Z][0-9]{5}[A-Z]{0,1}_[0-9]{3}'))
    }



    acs_5yr_c <- tryCatch(expr = {
        acs_5yr %>% 
        dplyr::select(-moe) %>% 
        pivot_wider(id_cols = c(GEOID, NAME), names_from = variable, values_from = c(estimate)) %>%
        transmute(GEOID = GEOID,
              n_pop_total = B01001_001,
              n_medincome = B19013_001,
              p_nonwhite = 100 * (B01001_001 - B01001A_001) / B01001_001,
              # i_gini = B19083_001,
              p_elderly = 100 * (B01001_020+B01001_021+B01001_022+B01001_023+B01001_024+B01001_025+B01001_044+B01001_045+B01001_046+B01001_047+B01001_048+B01001_049)/B01001_001,
              p_youngadult = 100 * (B01001_007+B01001_008+B01001_009+B01001_010+B01001_031+B01001_032+B01001_033+B01001_034)/B01001_001,
              #18-24
              p_poverty = 100 * B17001_002/B17001_001,
              #p_edub9 =  100 * (B15003_002+B15003_003+B15003_004+B15003_005+B15003_006+B15003_007+B15003_008+B15003_009+B15003_010+B15003_011+B15003_012+B15003_013)/B15003_001,
              p_edubac = 100 * (B15002_015 + B15002_016 + B15002_017 + B15002_018 + B15002_032 + B15002_033 + B15002_034 + B15002_035)/B15002_001,
              p_unemp := !!rlang::parse_quo(str_c('100*(',str_c(b23001_names_unemp, collapse = '+'), ')/(', str_c(b23001_names_inlabor, collapse = '+'), ')'), environment()),
              #p_unemp = 100 * B23025_005/B23025_002, # denominator is the population in labor force
              p_renter = 100 * (B25003_003/B25003_001), # denominator is the total number of occupied units
              p_vacant = 100 * (B25004_001)/B25001_001,
              p_singleheaded = 100 * (B25009_003+B25009_011)/B25009_001,
              p_noncitizen = 100 * (B05001_006/B05001_001),
              p_veteran = 100 * (B21001_002 / B21001_001),
              p_bmarital = 100 * (B12002_035+B12002_050+B12002_065+B12002_080+B12002_129+B12002_143+B12002_158+B12002_173)/B12002_001,
              #'B12002_001', 'B12002_035', 'B12002_050', 'B12002_065', 'B12002_080', 'B12002_129', 'B12002_143', 'B12002_158', 'B12002_173',
              p_fulltime = 100 * (B20005_003 + B20005_050)/B20005_001,
              p_uninsured = ifelse(any(grepl('B27001', acs_5yr$variable)>=1), 
                                100 * ((B27001_005+B27001_008+B27001_011+B27001_014+B27001_017+B27001_020+B27001_023+B27001_026+B27001_029+B27001_033+B27001_036+B27001_039+B27001_042+B27001_045+B27001_048+B27001_051+B27001_054+B27001_057)/B27001_001),
                                NA)
                )},
            error = function(e) {print(e) 
                                 return(acs_5yr)}
    )

    acs_5yr_moe <- acs_5yr %>% 
        dplyr::select(-estimate) %>% 
        pivot_wider(id_cols = c(GEOID, NAME), names_from = variable, values_from = c(moe)) %>%
        rename_at(.vars = vars(-1:-2), .funs = list(~str_c(., '_moe')))

    return(list(est = acs_5yr_c, moe = acs_5yr_moe))

}

# acs09_tr = get_clear_acs(year = 2009, state = 'Washington', variable_list = code_list_10)
# acs01_tr = get_clear_acs(year = 2009, state = 'Washington', variable_list = code_list_10)
# acs15_tr = get_clear_acs(year = 2015, state = 'Washington', variable_list = code_list)
# acs18_tr = get_clear_acs(year = 2018, state = 'Oregon', variable_list = code_list)

get_allyrs_acs = function(years, state, variable_list_old, variable_list, geography = 'county') {
    yrlist = years %>% split(.,.)
    acs_yrs = yrlist %>%
        lapply(function(x) {
                get_clear_acs(geography = geography, year = x, state = state, variable_list = if (x <= 2011) variable_list_old else variable_list)
        })
    return(acs_yrs)
}


years = 2009:2019 %>% split(.,.)

or_cnty_0919 = get_allyrs_acs(2009:2019, 'Oregon', code_list_10, code_list, geography = 'county')
wa_cnty_0919 = get_allyrs_acs(2009:2019, 'Washington', code_list_10, code_list, geography = 'county')

# or_cnty_pop_0919 = get_allyrs_acs(2009:2019, 'Oregon', acs_variables$name[1:49], acs_variables$name[1:49])
# wa_cnty_pop_0919 = get_allyrs_acs(2009:2019, 'Washington', acs_variables$name[1:49], acs_variables$name[1:49])

or_cnty_0919d = mapply(function(x, y) x$est %>% mutate(year = y),
                        or_cnty_0919, 2009:2019 %>% split(.,.), SIMPLIFY = FALSE) %>%
    Reduce(rbind, .)
wa_cnty_0919d = mapply(function(x, y) x$est %>% mutate(year = y),
                        wa_cnty_0919, 2009:2019 %>% split(.,.), SIMPLIFY = FALSE) %>%
    Reduce(rbind, .)
orwa_cnty_0919d = bind_rows(or_cnty_0919d, wa_cnty_0919d)

orwa_cnty_0919df = bind_rows(
    orwa_cnty_0919d %>% filter(year == 2009) %>% mutate(year = 2006),
    orwa_cnty_0919d %>% filter(year == 2009)%>% mutate(year = 2007),
    orwa_cnty_0919d %>% filter(year == 2009)%>% mutate(year = 2008),
    orwa_cnty_0919d %>% filter(year == 2009)%>% mutate(year = 2009),
    orwa_cnty_0919d %>% filter(year > 2009)
  ) %>%
  arrange(year, GEOID) %>%
  left_join(cpi_0619, by = c('year' = 'Year')) %>%
  group_by(GEOID) %>%
  mutate(n_medincome_adj = n_medincome * CPI10) %>%
  ungroup


# ZCTA (after 2011)
or_zcta_1119 = get_allyrs_acs(2011:2019, 'Oregon', code_list_10, code_list, geography = 'zcta')
wa_zcta_1119 = get_allyrs_acs(2011:2019, 'Washington', code_list_10, code_list, geography = 'zcta')

# or_zcta_pop_0919 = get_allyrs_acs(2009:2019, 'Oregon', acs_variables$name[1:49], acs_variables$name[1:49])
# wa_zcta_pop_0919 = get_allyrs_acs(2009:2019, 'Washington', acs_variables$name[1:49], acs_variables$name[1:49])

or_zcta_1119d = mapply(function(x, y) x$est %>% mutate(year = y),
                        or_zcta_1119, 2011:2019 %>% split(.,.), SIMPLIFY = FALSE) %>%
    Reduce(rbind, .)
wa_zcta_1119d = mapply(function(x, y) x$est %>% mutate(year = y),
                        wa_zcta_1119, 2011:2019 %>% split(.,.), SIMPLIFY = FALSE) %>%
    Reduce(rbind, .)
orwa_zcta_1119d = bind_rows(or_zcta_1119d, wa_zcta_1119d)

orwa_zcta_1119df = bind_rows(
    orwa_zcta_1119d %>% filter(year == 2011) %>% mutate(year = 2006),
    orwa_zcta_1119d %>% filter(year == 2011)%>% mutate(year = 2007),
    orwa_zcta_1119d %>% filter(year == 2011)%>% mutate(year = 2008),
    orwa_zcta_1119d %>% filter(year == 2011)%>% mutate(year = 2009),
    orwa_zcta_1119d %>% filter(year == 2011)%>% mutate(year = 2010),
    orwa_zcta_1119d %>% filter(year > 2011)
)

save(orwa_cnty_0919df, orwa_zcta_1119df,
     cnty_orwa, zcta_orwa,
     file = str_c(datapath, "Data/orwa_county_zcta_acs_122622.RData"),
     compress = 'xz', compression_level = 9)



## Oregon data cleaning ####
## Oregon ####
ord <- read_rds('/home/felix/Documents/registry/ORD_2006_2019_Raw.rds')

ord_s <- ord %>% 
#ord <- ord %>% 
  # Refactorize
  mutate(mental_underly = ifelse(grepl('F.*', dmcaACME), 1, 0)) %>% 
  mutate(ageunit = plyr::mapvalues(ddageunit, c(1,2,4,5,6,9), c(1,2,3,4,5,9)),
         married = plyr::mapvalues(ddmarital, c("A","D","M","S","W","U"), c("U","D","M","S","W","U")),
         educ = ddeduc,
         sex = ddsex,
         pregnancy = plyr::mapvalues(ddpregoutcome, c(1,2,3,4,7,8,9), c(1,2,3,4,7,8,9)),
         smoking = ddtobacco,
         dth_year = ddodyear,
         dth_month = ddodmonth %>% as.integer,
         dth_day = ddodday
  ) %>% 
  mutate(age = case_when(
    ageunit %in% c(5:6) ~ 0,
    ageunit == 4 ~ floor((1 * ddagenum)/365.24),
    ageunit == 2 ~ floor((30 * ddagenum)/365.24),
    ageunit == 1 ~ ddagenum,
    ageunit == 9 ~ 9999,
    ddagenum == 999 ~ 9999,
    TRUE ~ 9999
  )) %>% 
  mutate(race = case_when(
    ddracewh=='Y' & ddethnicmex!='H' & ddethnicpr!='H' & ddethniccuban!='H' & ddethnicoth!='H' ~ 1,
    ddracebl=='Y' ~ 2,
    ddethnicmex=='H' | ddethnicpr=='H' | ddethniccuban=='H' | ddethnicoth=='H' ~ 3,
    ddethnicmex=='U' | ddethnicpr =='U'| ddethniccuban=='U' | ddethnicoth=='U' | ddracewh=='U' | ddracebl=='U' | ddraceaian=='U' | ddraceasianind=='U' |
    ddracech=='U' | ddracefi=='U' | ddracekor=='U' | ddracevt=='U' | ddraceoasian=='U' | ddracenh =='U' | ddracegu=='U' | ddracesm=='U' |
    ddraceopi=='U' | ddraceospf == 'U' ~ 9,
    TRUE ~ 4
  )) %>% 
  rename_with(.cols = colnames(.)[grep('^dmca([1-9]|1[0-9]|2[0])$', colnames(.))],
              .fn = ~gsub('dmca', 'mltcse', .x, fixed = TRUE)) %>% 
  mutate(mental = grepl('F.*', mltcse1),
         mental3 = as.logical(grepl('F.*', mltcse1)+grepl('F.*', mltcse2)+grepl('F.*', mltcse3)),
         mental_new = ifelse(grepl('F[0-4][0-9]', mltcse1) + grepl('F[0-4][0-9]', mltcse2) + grepl('F[0-4][0-9]', mltcse3) + 
                               grepl('F[0-4][0-9]', mltcse4) + grepl('F[0-4][0-9]', mltcse5) + grepl('F[0-4][0-9]', mltcse6) +
                               grepl('F[0-4][0-9]', mltcse7) + grepl('F[0-4][0-9]', mltcse8) + grepl('F[0-4][0-9]', mltcse9) + grepl('F[0-4][0-9]', mltcse10), 1, 0),
         mental_new2 = ifelse(grepl('F[0-6][0-9]', mltcse1) + grepl('F[0-6][0-9]', mltcse2) + grepl('F[0-6][0-9]', mltcse3) + 
                               grepl('F[0-6][0-9]', mltcse4) + grepl('F[0-6][0-9]', mltcse5) + grepl('F[0-6][0-9]', mltcse6) +
                               grepl('F[0-6][0-9]', mltcse7) + grepl('F[0-6][0-9]', mltcse8) + grepl('F[0-6][0-9]', mltcse9) + grepl('F[0-6][0-9]', mltcse10), 1, 0),
         substance = ifelse(grepl('F1[0-9]', mltcse1) + grepl('F1[0-9]', mltcse2) + grepl('F1[0-9]', mltcse3) + 
                               grepl('F1[0-9]', mltcse4) + grepl('F1[0-9]', mltcse5) + grepl('F1[0-9]', mltcse6) +
                               grepl('F1[0-9]', mltcse7) + grepl('F1[0-9]', mltcse8) + grepl('F1[0-9]', mltcse9) + grepl('F1[0-9]', mltcse10), 1, 0),
         nonsubstance = ifelse(grepl('F[2-4][0-9]', mltcse1) + grepl('F[2-4][0-9]', mltcse2) + grepl('F[2-4][0-9]', mltcse3) + 
                               grepl('F[2-4][0-9]', mltcse4) + grepl('F[2-4][0-9]', mltcse5) + grepl('F[2-4][0-9]', mltcse6) +
                               grepl('F[2-4][0-9]', mltcse7) + grepl('F[2-4][0-9]', mltcse8) + grepl('F[2-4][0-9]', mltcse9) + grepl('F[2-4][0-9]', mltcse10), 1, 0),
         nonsubstance2 = ifelse(grepl('F(0|[2-6])[0-9]', mltcse1) + grepl('F(0|[2-6])[0-9]', mltcse2) + grepl('F(0|[2-6])[0-9]', mltcse3) + 
                               grepl('F(0|[2-6])[0-9]', mltcse4) + grepl('F(0|[2-6])[0-9]', mltcse5) + grepl('F(0|[2-6])[0-9]', mltcse6) +
                               grepl('F(0|[2-6])[0-9]', mltcse7) + grepl('F(0|[2-6])[0-9]', mltcse8) + grepl('F(0|[2-6])[0-9]', mltcse9) + grepl('F(0|[2-6])[0-9]', mltcse10), 1, 0),
         moodanxiety = ifelse(grepl('F[3-4][0-9]', mltcse1) + grepl('F[3-4][0-9]', mltcse2) + grepl('F[3-4][0-9]', mltcse3) + 
                               grepl('F[3-4][0-9]', mltcse4) + grepl('F[3-4][0-9]', mltcse5) + grepl('F[3-4][0-9]', mltcse6) +
                               grepl('F[3-4][0-9]', mltcse7) + grepl('F[3-4][0-9]', mltcse8) + grepl('F[3-4][0-9]', mltcse9) + grepl('F[3-4][0-9]', mltcse10), 1, 0),
         parkinson = ifelse(grepl('^(A521|F023|G20|G21|G22|G903).*', mltcse1) | grepl('^(A521|F023|G20|G21|G22|G903).*', mltcse2) | grepl('^(A521|F023|G20|G21|G22|G903).*', mltcse3), 1, 0),
         dementia = ifelse(grepl('^(F00|F01|F02|F03).*', mltcse1) | grepl('^(F00|F01|F02|F03).*', mltcse2) | grepl('^(F00|F01|F02|F03).*', mltcse3), 1, 0),
         dementia_vascular = ifelse(grepl('F01', mltcse1) + grepl('F01', mltcse2) + grepl('F01', mltcse3) + 
                               grepl('F01', mltcse4) + grepl('F01', mltcse5) + grepl('F01', mltcse6) +
                               grepl('F01', mltcse7) + grepl('F01', mltcse8) + grepl('F01', mltcse9) + grepl('F01', mltcse10), 1, 0),
         dementia_unspecified = ifelse(grepl('F03', mltcse1) + grepl('F03', mltcse2) + grepl('F03', mltcse3) + 
                               grepl('F03', mltcse4) + grepl('F03', mltcse5) + grepl('F03', mltcse6) +
                               grepl('F03', mltcse7) + grepl('F03', mltcse8) + grepl('F03', mltcse9) + grepl('F03', mltcse10), 1, 0),
         alzheimer = ifelse(grepl('^(G30).*', mltcse1) | grepl('^(G30).*', mltcse2) | grepl('^(G30).*', mltcse3), 1, 0),
         educ_f = factor(educ),
         married_f = factor(married),
         smoking_f = factor(smoking),
         race_f = factor(plyr::mapvalues(race, LETTERS[1:8], rep('Others',8)))) %>% 
  mutate(mental_u = as.integer(ifelse(grepl('F.*', dmcaACME), 1, 0)),
         parkinson_u = as.integer(ifelse(grepl('^(A521|F023|G20|G21|G22|G903).*', dmcaACME), 1, 0)),
         dementia_u = as.integer(ifelse(grepl('^(F00|F01|F02|F03).*', dmcaACME), 1, 0)),
         alzheimer_u = as.integer(ifelse(grepl('^(G30).*', dmcaACME), 1, 0)))




## Clear coordinates
ord_sweird <- ord_s %>% 
  filter(ddreslong >= 180) %>% 
  st_as_sf(coords = c('ddreslong', 'ddreslat'), crs = 2992) %>% 
  st_transform(4326) %>% 
  mutate(ddreslong = st_coordinates(.)[,1],
         ddreslat = st_coordinates(.)[,2]) %>% 
  st_set_geometry(NULL)

ord_snormal <- ord_s %>% 
  filter(ddreslong <= 180)

ord_sre <- ord_sweird %>% 
    bind_rows(ord_snormal) %>% 
    filter(ddreslat >= 40 & ddreslong < -80)

ord_sf <- ord_sre %>% 
    st_as_sf(coords = c('ddreslong', 'ddreslat'), crs = 4326)



# or_tract_pop_0919_n = or_tract_pop_0919 %>%
#   lapply(function(x) { x$est %>% 
#       dplyr::select(1, 3, 4) %>% 
#       pivot_wider(names_from = variable, 
#                   values_from = estimate) %>%
#       transmute(GEOID = GEOID,
#              n_00_05 = B01001_003 + B01001_027,
#              n_05_10 = B01001_004 + B01001_028,
#              n_10_15 = B01001_005 + B01001_029,
#              n_15_20 = B01001_006 + B01001_007 + B01001_030 + B01001_031,
#              n_20_25 = B01001_008 + B01001_009 + B01001_010 + B01001_032 + B01001_033 + B01001_034,
#              n_25_30 = B01001_011 + B01001_035,
#              n_30_35 = B01001_012 + B01001_036,
#              n_35_40 = B01001_013 + B01001_037,
#              n_40_45 = B01001_014 + B01001_038,
#              n_45_50 = B01001_015 + B01001_039,
#              n_50_55 = B01001_016 + B01001_040,
#              n_55_60 = B01001_017 + B01001_041,
#              n_60_65 = B01001_018 + B01001_019 + B01001_042 + B01001_043,
#              n_65_70 = B01001_020 + B01001_021 + B01001_044 + B01001_045,
#              n_70_75 = B01001_022 + B01001_046,
#              n_75_80 = B01001_023 + B01001_047,
#              n_80_85 = B01001_024 + B01001_048,
#              n_85_Inf = B01001_025 + B01001_049)}) %>%
#   mapply(function(x, y) { x %>% mutate(year = y) %>%
#               pivot_longer(cols = 2:(ncol(.) - 1))}, 
#               ., 2009:2019 %>% split(.,.), SIMPLIFY = FALSE) %>%
#   do.call(bind_rows, .)

# # 10 year basis

# or_tract_pop_0919_n = or_tract_pop_0919 %>%
#   lapply(function(x) { x$est %>% 
#       dplyr::select(1, 3, 4) %>% 
#       pivot_wider(names_from = variable, 
#                   values_from = estimate) %>%
#       transmute(GEOID = GEOID,
#              n_00_10 = B01001_003 + B01001_027 + B01001_004 + B01001_028,
#              n_10_20 = B01001_005 + B01001_029 + B01001_006 + B01001_007 + B01001_030 + B01001_031,
#              n_20_30 = B01001_008 + B01001_009 + B01001_010 + B01001_032 + B01001_033 + B01001_034 + B01001_011 + B01001_035,
#              n_30_40 = B01001_012 + B01001_036 + B01001_013 + B01001_037,
#              n_40_50 = B01001_014 + B01001_038 + B01001_015 + B01001_039,
#              n_50_60 = B01001_016 + B01001_040 + B01001_017 + B01001_041,
#              n_60_70 = B01001_018 + B01001_019 + B01001_042 + B01001_043 + B01001_020 + B01001_021 + B01001_044 + B01001_045,
#              n_70_80 = B01001_022 + B01001_046 + B01001_023 + B01001_047,
#              n_80_Inf = B01001_024 + B01001_048 + B01001_025 + B01001_049)})
# or_tract_pop_0919_n[[1]] = or_tract_pop_0919_n[[1]] %>%
#   full_join(cw_or %>% mutate(GEOID00 = as.character(GEOID00)), by = c('GEOID' = 'GEOID00')) %>%
#   group_by(GEOID10) %>%
#   summarize_at(.vars = vars(n_00_10:n_80_Inf),
#                .funs = list(~floor(sum(. * ((POP00 * POPPCT00 / 100)/ sum(POP00 * POPPCT00 / 100)))))) %>%
#   ungroup %>%
#   mutate(GEOID = as.character(GEOID10)) %>%
#   dplyr::select(-GEOID10) %>%
#   dplyr::select(GEOID, n_00_10:n_80_Inf)


# or_tract_pop_0919_df = or_tract_pop_0919_n %>%
#   mapply(function(x, y) { x %>% mutate(year = y) %>%
#               pivot_longer(cols = 2:(ncol(.) - 1))}, 
#               ., 2009:2019 %>% split(.,.), SIMPLIFY = FALSE) %>%
#   do.call(bind_rows, .)

# or_tract_pop_0619_df =
#   bind_rows(
#       or_tract_pop_0919_df %>% filter(year == 2009) %>% mutate(year = 2006),
#       or_tract_pop_0919_df %>% filter(year == 2009) %>% mutate(year = 2007),
#       or_tract_pop_0919_df %>% filter(year == 2009) %>% mutate(year = 2008),
#       or_tract_pop_0919_df) %>%
#   rename(agegroup = name,
#          tpop = value) %>%
#   mutate(agegroup = str_replace_all(agegroup, 'n_', ''))

# # total population
# or_tract_pop_0619_total = 
#   or_tract_pop_0619_df %>%
#   group_by(GEOID, year) %>%
#   summarize(tpop = sum(tpop)) %>%
#   # NA filling with 2009 population (tract 9705 in Malheur county)
#   mutate(tpop = ifelse(is.na(tpop), tpop[which(!is.na(tpop))[1]], tpop)) %>%
#   ungroup

  
# tract data (assumption: the working directory points the dissertation base directory)
or_tract = 
  st_read("./Data/Basemap/ORWA_Tracts_2010.gpkg") %>%
  filter(grepl('^41', GEOID10)) %>%
  rmapshaper::ms_simplify(keep = 0.33, keep_shapes = TRUE)
or_tract_lite = or_tract %>%
  dplyr::select(GEOID10) %>%
  st_transform(4326)
ord_sf_tr = ord_sf %>%
  st_join(or_tract_lite)
ord_sf_tr_agecounts = ord_sf_tr %>%
  st_drop_geometry %>%
  mutate(agegroup = cut(age, c(seq(0, 80, 10), Inf), right = FALSE),
         agegroup = str_replace_all(agegroup, '[[:punct:]]', '_')) %>%
  group_by(GEOID10, dth_year, agegroup) %>%
  summarize(n_mental1 = sum(mental, na.rm = TRUE),
            n_mental3 = sum(mental3, na.rm = TRUE),
            n_mentaln = sum(mental_new, na.rm = TRUE),
            n_mentaln2 = sum(mental_new2, na.rm = TRUE),
            n_moodanxiety = sum(moodanxiety, na.rm = TRUE),
            n_substance = sum(substance, na.rm = TRUE),
            n_nonsubstance = sum(nonsubstance, na.rm = TRUE),
            n_nonsubstance2 = sum(nonsubstance2, na.rm = TRUE),
            n_dementia_vas = sum(dementia_vascular, na.rm = TRUE),
            n_dementia_uns = sum(dementia_unspecified, na.rm = TRUE)) %>%
  ungroup %>%
  mutate(agegroup = str_sub(agegroup, 2, nchar(agegroup)-1))
ord_sf_tr_summary = ord_sf_tr %>%
  st_drop_geometry %>%
  group_by(GEOID10, dth_year) %>%
  summarize(n_mental1 = sum(mental, na.rm = TRUE),
            n_mental3 = sum(mental3, na.rm = TRUE),
            n_mentaln = sum(mental_new, na.rm = TRUE),
            n_mentaln2 = sum(mental_new2, na.rm = TRUE),
            n_moodanxiety = sum(moodanxiety, na.rm = TRUE),
            n_substance = sum(substance, na.rm = TRUE),
            n_nonsubstance = sum(nonsubstance, na.rm = TRUE),
            n_nonsubstance2 = sum(nonsubstance2, na.rm = TRUE),
            n_dementia_vas = sum(dementia_vascular, na.rm = TRUE),
            n_dementia_uns = sum(dementia_unspecified, na.rm = TRUE),
            n_meanage_nonsubstance2 = sum(age * nonsubstance2, na.rm = T) / sum(nonsubstance2, na.rm = T),
            n_meanage_dementia_vas = sum(age * dementia_vascular, na.rm = T) / sum(dementia_vascular, na.rm = T),
            n_meanage_dementia_uns = sum(age * dementia_unspecified, na.rm = T) / sum(dementia_unspecified, na.rm = T)) %>%
  ungroup

pop_agestr = tribble(
  ~agegroup,  ~n_all, ~n_male,  ~n_female,  ~newgroup,
  "00_05",  20201362, 10319427, 9881935,  "00_10",
  "05_10", 20348657, 10389638, 9959019, "00_10",
  "10_15",  20677194, 10579862, 10097332, "10_20",
  "15_20",  22040343, 11303666, 10736677, "10_20",
  "20_25",  21585999, 11014176, 10571823, "20_30",
  "25_30",  21101849, 10635591, 10466258, "20_30",
  "30_35",  19962099, 9996500,  9965599,  "30_40",
  "35_40",  20179642, 10042022, 10137620, "30_40",
  "40_45",  20890964, 10393977, 10496987, "40_50",
  "45_50",  22708591, 11209085, 11499506, "40_50",
  "50_55",  22298125, 10933274, 11364851, "50_60",
  "55_60",  19664805, 9523648,  10141157, "50_60",
  "60_65",  16817924, 8077500,  8740424,  "60_70",
  "65_70",  12435263, 5852547,  6582716,  "60_70",
  "70_75",  9278166,  4243972,  5034194,  "70_80",
  "75_80",  7317795,  3182388,  4135407,  "70_80",
  "80_85",  5743327,  2294374,  3448953,  "80_Inf",
  "85_Inf", 5493433,  1789679,  3703754,  "80_Inf",
)

# CPI (from Minneapolis Fed; base 1982-84)
# cpi = rvest::read_html("https://www.minneapolisfed.org/about-us/monetary-policy/inflation-calculator/consumer-price-index-1913-") %>%
#   rvest::html_table() %>%
#   .[[1]] %>%
#   .[,-3]
cpi = read_csv(str_c(datapath, "/Data/CPI_1913_2021.csv"))
cpi_0619 = cpi %>%
  filter(Year %in% 2006:2019) %>%
  mutate(Year = as.integer(Year),
         CPI10 = `Annual Average`[5] / `Annual Average`) %>%
  dplyr::select(1,3)

pop_agestr = pop_agestr %>%
  # regroup
  group_by(newgroup) %>%
  summarize_at(.vars = vars(n_all:n_female),
               .funs = list(~sum(.))) %>%
  ungroup %>%
  rename(agegroup = newgroup) %>%
  # regroup end
  pivot_longer(cols = 2:4) %>%
  group_by(name) %>% 
  mutate_at(vars(value), list(~(./sum(.)))) %>%
  ungroup %>%
  pivot_wider(names_from = name,
              values_from = value)
or_tract_asmr_base = or_tract_pop_0619_df %>%
  left_join(ord_sf_tr_agecounts, by = c('GEOID' = 'GEOID10', 'year' = 'dth_year', 'agegroup' = 'agegroup')) %>%
  left_join(pop_agestr, by = c('agegroup' = 'agegroup')) %>%
  mutate_at(.vars = vars(n_mental1:n_substance),
            .funs = list(~ifelse(is.na(.), 0, .))) %>%
  mutate(r_mental1_100k = 1e5 * (n_mental1 / tpop),
         r_mental3_100k = 1e5 * (n_mental3 / tpop),
         r_mentaln_100k = 1e5 * (n_mentaln / tpop))
or_tract_asmr = or_tract_asmr_base %>%
  group_by(GEOID, year) %>%
  summarize_at(.vars = vars(r_mental1_100k:r_mentaln_100k),
               .funs = list(~sum(.*n_all, na.rm = TRUE))) %>%
  ungroup
or_tract_n = or_tract_asmr_base %>%
  group_by(GEOID, year) %>%
  summarize_at(.vars = vars(n_mental1:n_substance),
               .funs = list(~sum(., na.rm = T))) %>%
  ungroup

or_tract_asmr_gg = ggplot(or_tract_asmr,
                          mapping = aes(x = year, y = r_mental3_100k)) +
                  geom_line(alpha = 0.3, lwd = 0.66, mapping = aes(group = GEOID)) +
                  stat_summary(geom = 'line', fun = mean, color = 'red', lwd = 1.5)
or_tract_n_gg = ggplot(or_tract_n,
                          mapping = aes(x = year, y = d_mental3)) +
                  geom_line(alpha = 0.3, lwd = 0.66, mapping = aes(group = GEOID)) +
                  stat_summary(geom = 'line', fun = mean, color = 'red', lwd = 1.5)

###
# composite map
2006:2019 %>% 
  split(.,.) %>%
  lapply(function(x) or_tract_lite %>% st_transform(3857) %>% mutate(year = x)) %>%
  do.call(rbind, .) -> or_tracts_0619
or_tracts_0619_md = or_tracts_0619 %>%
  left_join(or_tract_n, by = c('GEOID10' = 'GEOID', 'year' = 'year'))


## ready-made ACS covariates (00-10 population weighted mean) ####
or_tract_0919_d = 2009:2019 %>%
  split(.,.) %>%
  mapply(function(x, y) x$est %>% mutate(year = y),
         or_tract_0919, ., SIMPLIFY = FALSE)
or_tract_0919_d[[1]] = or_tract_0919_d[[1]] %>%
  full_join(cw_or %>% mutate(GEOID00 = as.character(GEOID00)), by = c('GEOID' = 'GEOID00')) %>%
  group_by(GEOID10) %>%
  summarize_at(.vars = vars(n_medincome:p_uninsured),
               .funs = list(~sum(sum(. * ((POP00 * POPPCT00 / 100)/ sum(POP00 * POPPCT00 / 100)))))) %>%
  ungroup %>%
  mutate(GEOID = as.character(GEOID10)) %>%
  dplyr::select(-GEOID10) %>%
  dplyr::select(GEOID, n_medincome:p_uninsured)
or_tract_0919_df = or_tract_0919_d %>%
  do.call(bind_rows, .) %>%
  mutate(year = ifelse(is.na(year), 2009, year))
or_tract_0619_df = 
  bind_rows(
    or_tract_0919_d[[1]] %>% mutate(year = 2006),
    or_tract_0919_d[[1]] %>% mutate(year = 2007),
    or_tract_0919_d[[1]] %>% mutate(year = 2008),
    or_tract_0919_df)    




## Washington data cleaning ####

## Location filtering
wad = read_rds('/home/felix/Documents/registry/WAD_DF_111420.rds')
wad_sf = wad %>% 
  filter(!is.na(long) & !is.na(lat)) %>% 
  filter(long < 0 & lat > 40) %>% 
  st_as_sf(coords = c('long', 'lat'), crs = 4326) %>%
  mutate(
        record_year = as.integer(str_sub(certno, 1, 4)),
        fac_type = as.numeric(as.character(fac_type)),
        educ = as.numeric(as.character(educ)),
        #ageunit = plyr::mapvalues(ageunit, c(0,1,2,4,5,6,9), c(1,1,2,3,4,5,9)),
         fac_type = 
           case_when(
             record_year > 2003 & record_year <= 2015 ~ plyr::mapvalues(fac_type, c(0,1,2,3,4,5,6,7,9), c(0,1,1,3,4,5,4,5,9)),
             record_year <= 2003 ~ plyr::mapvalues(fac_type, c(0,1,2,3,4,5,6,7,8,9), c(0,1,2,3,4,5,4,5,1,9)),
             TRUE ~ fac_type
           ),
         armforce = plyr::mapvalues(armforce, c(1,2,9), c('Y', 'N', 'U')),
         married = 
           case_when(
             record_year > 2003 & record_year <= 2015 ~ plyr::mapvalues(married, c(1:6,9), c('S', 'M', 'D', 'W', 'M', 'M', 'U')),
             record_year <= 2003 ~ plyr::mapvalues(married, c(1,2,3,4,9), c('S', 'M', 'D', 'W', 'U')),
             TRUE ~ married
             ),
         educ =
           case_when(
             record_year > 2003 & record_year <= 2015 ~ plyr::mapvalues(educ, 1:9, c(1,2,3,4,4,5,7,7,9)),
             record_year <= 2003 ~ plyr::mapvalues(educ, c(0:17,99), c(1,1,1,1,1,1,1,1,1,2,2,2,3,4,4,4,6,7,9)),
             TRUE ~ educ
             ),
         resunit = ifelse(is.na(resunit), res_lena, resunit),
         resunit = 
           ifelse(record_year <= 2015, plyr::mapvalues(resunit, c(0,1,2,3,4,5,6,9), c('Y','Y','M', 'W', 'D', 'H', 'M', 'U')), resunit),
         #referred = ifelse(record_year <= 2015, plyr::mapvalues(referred, c(1,2,9), c('Y', 'N', 'U')), referred),
         #smoking = ifelse(is.na(smoking), tbcontri, smoking),
         smoking = case_when(
           record_year <= 2003 ~ plyr::mapvalues(smoking, c(1,2,9), c('Y', 'N', 'U')),
           record_year > 2003 & record_year <= 2015 ~ plyr::mapvalues(smoking, c(1,2,7,8,9), c('Y', 'N', 'P', 'U','U')),
           TRUE ~ smoking
         )
         )

wad_sf <- wad_sf %>%
  filter(dth_yr >= 2006) %>% 
#   filter(!is.na(long) & !is.na(lat)) %>% 
#   filter(long < 0 & lat > 40) %>% 
  mutate(resunum = as.integer(resunum)) %>% 
  mutate(parkinson = ifelse(grepl('^(A521|F023|G20|G21|G22|G903).*', mltcse1), 'Parkinson', 'Others'),
         mental = ifelse(grepl('F.*', mltcse1), 1, 0),
         mental2 = ifelse(grepl('F.*', mltcse1) + grepl('F.*', mltcse2) + grepl('F.*', mltcse3), 1, 0),
         mental3 = ifelse(grepl('F.*', mltcse1) + grepl('F.*', mltcse2) + grepl('F.*', mltcse3) + 
                        grepl('F.*', mltcse4) + grepl('F.*', mltcse5) + grepl('F.*', mltcse6) +
                        grepl('F.*', mltcse7) + grepl('F.*', mltcse8) + grepl('F.*', mltcse9) + grepl('F.*', mltcse10), 1, 0),
         substance = ifelse(grepl('F1[0-9]', mltcse1) + grepl('F1[0-9]', mltcse2) + grepl('F1[0-9]', mltcse3) + 
                               grepl('F1[0-9]', mltcse4) + grepl('F1[0-9]', mltcse5) + grepl('F1[0-9]', mltcse6) +
                               grepl('F1[0-9]', mltcse7) + grepl('F1[0-9]', mltcse8) + grepl('F1[0-9]', mltcse9) + grepl('F1[0-9]', mltcse10), 1, 0),
         mental_new2 = ifelse(grepl('F[0-6][0-9]', mltcse1) + grepl('F[0-6][0-9]', mltcse2) + grepl('F[0-6][0-9]', mltcse3) + 
                               grepl('F[0-6][0-9]', mltcse4) + grepl('F[0-6][0-9]', mltcse5) + grepl('F[0-6][0-9]', mltcse6) +
                               grepl('F[0-6][0-9]', mltcse7) + grepl('F[0-6][0-9]', mltcse8) + grepl('F[0-6][0-9]', mltcse9) + grepl('F[0-6][0-9]', mltcse10), 1, 0),
         nonsubstance = ifelse(grepl('F[2-4][0-9]', mltcse1) + grepl('F[2-4][0-9]', mltcse2) + grepl('F[2-4][0-9]', mltcse3) + 
                               grepl('F[2-4][0-9]', mltcse4) + grepl('F[2-4][0-9]', mltcse5) + grepl('F[2-4][0-9]', mltcse6) +
                               grepl('F[2-4][0-9]', mltcse7) + grepl('F[2-4][0-9]', mltcse8) + grepl('F[2-4][0-9]', mltcse9) + grepl('F[2-4][0-9]', mltcse10), 1, 0),
         nonsubstance2 = ifelse(grepl('F(0|[2-6])[0-9]', mltcse1) + grepl('F(0|[2-6])[0-9]', mltcse2) + grepl('F(0|[2-6])[0-9]', mltcse3) + 
                               grepl('F(0|[2-6])[0-9]', mltcse4) + grepl('F(0|[2-6])[0-9]', mltcse5) + grepl('F(0|[2-6])[0-9]', mltcse6) +
                               grepl('F(0|[2-6])[0-9]', mltcse7) + grepl('F(0|[2-6])[0-9]', mltcse8) + grepl('F(0|[2-6])[0-9]', mltcse9) + grepl('F(0|[2-6])[0-9]', mltcse10), 1, 0),
         moodanxiety = ifelse(grepl('F[3-4][0-9]', mltcse1) + grepl('F[3-4][0-9]', mltcse2) + grepl('F[3-4][0-9]', mltcse3) + 
                               grepl('F[3-4][0-9]', mltcse4) + grepl('F[3-4][0-9]', mltcse5) + grepl('F[3-4][0-9]', mltcse6) +
                               grepl('F[3-4][0-9]', mltcse7) + grepl('F[3-4][0-9]', mltcse8) + grepl('F[3-4][0-9]', mltcse9) + grepl('F[3-4][0-9]', mltcse10), 1, 0),
         mental_new = ifelse(grepl('F[0-4][0-9]', mltcse1) + grepl('F[0-4][0-9]', mltcse2) + grepl('F[0-4][0-9]', mltcse3) + 
                               grepl('F[0-4][0-9]', mltcse4) + grepl('F[0-4][0-9]', mltcse5) + grepl('F[0-4][0-9]', mltcse6) +
                               grepl('F[0-4][0-9]', mltcse7) + grepl('F[0-4][0-9]', mltcse8) + grepl('F[0-4][0-9]', mltcse9) + grepl('F[0-4][0-9]', mltcse10), 1, 0),
         dementia_vascular = ifelse(grepl('F01', mltcse1) + grepl('F01', mltcse2) + grepl('F01', mltcse3) + 
                               grepl('F01', mltcse4) + grepl('F01', mltcse5) + grepl('F01', mltcse6) +
                               grepl('F01', mltcse7) + grepl('F01', mltcse8) + grepl('F01', mltcse9) + grepl('F01', mltcse10), 1, 0),
         dementia_unspecified = ifelse(grepl('F03', mltcse1) + grepl('F03', mltcse2) + grepl('F03', mltcse3) + 
                               grepl('F03', mltcse4) + grepl('F03', mltcse5) + grepl('F03', mltcse6) +
                               grepl('F03', mltcse7) + grepl('F03', mltcse8) + grepl('F03', mltcse9) + grepl('F03', mltcse10), 1, 0)) %>%
  mutate(mental_underly = ifelse(grepl('F.*', underly), 1, 0)) %>% 
  mutate(resdays = case_when(
        resunit %in% c('H', 'N') ~ 1,
        resunit == 'D' ~ (1 * resunum),
        resunit == 'W' ~ (7 * resunum),
        resunit == 'M' ~ (30 * resunum),
        resunit == 'Y' ~ (365 * resunum),
        resunit == 'U' ~ as.numeric(resunum)
  )) %>% 
  mutate(race = factor(plyr::mapvalues(race, c(1:8, LETTERS[1:8]), c(1,2,rep(4,8), 3, rep(4,5))))) %>% 
  mutate(educ = factor(educ),
         race = factor(race),
         sex = factor(sex)) %>%
  mutate_if(is.factor, droplevels) 




# tract data
wa_tract = 
  st_read("./Data/Basemap/ORWA_Tracts_2010.gpkg") %>%
  filter(grepl('^53', GEOID10)) %>%
  rmapshaper::ms_simplify(keep = 0.33, keep_shapes = TRUE)
wa_tract_lite = wa_tract %>%
  dplyr::select(GEOID10) %>%
  st_transform(4326)
wad_sf_tr = wad_sf %>%
  st_join(wa_tract_lite)
wad_sf_tr_agecounts = wad_sf_tr %>%
  st_drop_geometry %>%
  mutate(agegroup = cut(age, c(seq(0, 80, 10), Inf), right = FALSE),
         agegroup = str_replace_all(agegroup, '[[:punct:]]', '_'),
         dth_year = dth_yr) %>%
  group_by(GEOID10, dth_year, agegroup) %>%
  summarize(n_mental1 = sum(mental, na.rm = TRUE),
            n_mental3 = sum(mental3, na.rm = TRUE),
            n_mentaln = sum(mental_new, na.rm = TRUE),
            n_mentaln2 = sum(mental_new2, na.rm = TRUE),
            n_moodanxiety = sum(moodanxiety, na.rm = TRUE),
            n_substance = sum(substance, na.rm = TRUE),
            n_nonsubstance = sum(nonsubstance, na.rm = TRUE),
            n_nonsubstance2 = sum(nonsubstance2, na.rm = TRUE),
            n_dementia_vas = sum(dementia_vascular, na.rm = TRUE),
            n_dementia_uns = sum(dementia_unspecified, na.rm = TRUE)) %>%
  ungroup %>%
  mutate(agegroup = str_sub(agegroup, 2, nchar(agegroup)-1))
wad_sf_tr_summary = wad_sf_tr %>%
  st_drop_geometry %>%
  mutate(dth_year = dth_yr) %>%
  group_by(GEOID10, dth_year) %>%
  summarize(n_mental1 = sum(mental, na.rm = TRUE),
            n_mental3 = sum(mental3, na.rm = TRUE),
            n_mentaln = sum(mental_new, na.rm = TRUE),
            n_mentaln2 = sum(mental_new2, na.rm = TRUE),
            n_moodanxiety = sum(moodanxiety, na.rm = TRUE),
            n_substance = sum(substance, na.rm = TRUE),
            n_nonsubstance = sum(nonsubstance, na.rm = TRUE),
            n_nonsubstance2 = sum(nonsubstance2, na.rm = TRUE),
            n_dementia_vas = sum(dementia_vascular, na.rm = TRUE),
            n_dementia_uns = sum(dementia_unspecified, na.rm = TRUE),
            n_meanage_nonsubstance2 = sum(age * nonsubstance2, na.rm = T) / sum(nonsubstance2, na.rm = T),
            n_meanage_dementia_vas = sum(age * dementia_vascular, na.rm = T) / sum(dementia_vascular, na.rm = T),
            n_meanage_dementia_uns = sum(age * dementia_unspecified, na.rm = T) / sum(dementia_unspecified, na.rm = T)) %>%
  ungroup



# One-shot processing ####
load(str_c(datapath, "Data/orwa_county_zcta_acs_122622.RData"))
or_cnty = 
  cnty_orwa %>%
  filter(grepl('^41', GEOID10)) %>%
  rmapshaper::ms_simplify(keep = 0.2, keep_shapes = TRUE)
or_cnty_lite = or_cnty %>%
  dplyr::select(GEOID10) %>%
  st_transform(4326)
ord_sf_tr = ord_sf %>%
  st_join(or_cnty_lite)
ord_sf_tr_summary = ord_sf_tr %>%
  st_drop_geometry %>%
  group_by(GEOID10, dth_year) %>%
  summarize(n_mental1 = sum(mental, na.rm = TRUE),
            n_mental3 = sum(mental3, na.rm = TRUE),
            n_mentaln = sum(mental_new, na.rm = TRUE),
            n_mentaln2 = sum(mental_new2, na.rm = TRUE),
            n_moodanxiety = sum(moodanxiety, na.rm = TRUE),
            n_substance = sum(substance, na.rm = TRUE),
            n_nonsubstance = sum(nonsubstance, na.rm = TRUE),
            n_nonsubstance2 = sum(nonsubstance2, na.rm = TRUE),
            n_dementia_vas = sum(dementia_vascular, na.rm = TRUE),
            n_dementia_uns = sum(dementia_unspecified, na.rm = TRUE),
            n_meanage_nonsubstance2 = sum(age * nonsubstance2, na.rm = T) / sum(nonsubstance2, na.rm = T),
            n_meanage_dementia_vas = sum(age * dementia_vascular, na.rm = T) / sum(dementia_vascular, na.rm = T),
            n_meanage_dementia_uns = sum(age * dementia_unspecified, na.rm = T) / sum(dementia_unspecified, na.rm = T)) %>%
  ungroup

wa_cnty_lite = cnty_orwa %>%
  filter(grepl('^53', GEOID10)) %>%
  rmapshaper::ms_simplify(keep = 0.2, keep_shapes = TRUE) %>%
  dplyr::select(GEOID10) %>%
  st_transform(4326)
wad_sf_tr = wad_sf %>%
  st_join(wa_cnty_lite)
wad_sf_tr_summary = wad_sf_tr %>%
  st_drop_geometry %>%
  mutate(dth_year = dth_yr) %>%
  group_by(GEOID10, dth_year) %>%
  summarize(n_mental1 = sum(mental, na.rm = TRUE),
            n_mental3 = sum(mental3, na.rm = TRUE),
            n_mentaln = sum(mental_new, na.rm = TRUE),
            n_mentaln2 = sum(mental_new2, na.rm = TRUE),
            n_moodanxiety = sum(moodanxiety, na.rm = TRUE),
            n_substance = sum(substance, na.rm = TRUE),
            n_nonsubstance = sum(nonsubstance, na.rm = TRUE),
            n_nonsubstance2 = sum(nonsubstance2, na.rm = TRUE),
            n_dementia_vas = sum(dementia_vascular, na.rm = TRUE),
            n_dementia_uns = sum(dementia_unspecified, na.rm = TRUE),
            n_meanage_nonsubstance2 = sum(age * nonsubstance2, na.rm = T) / sum(nonsubstance2, na.rm = T),
            n_meanage_dementia_vas = sum(age * dementia_vascular, na.rm = T) / sum(dementia_vascular, na.rm = T),
            n_meanage_dementia_uns = sum(age * dementia_unspecified, na.rm = T) / sum(dementia_unspecified, na.rm = T)) %>%
  ungroup


# sf features with ACS data
orwa_cnty_0618df = 
  cnty_orwa %>%
  group_by(STATEFP) %>%
  nest() %>%
  mutate(data = map(.x = data, .f = ~right_join(.x, orwa_cnty_0919df %>% filter(grepl(str_c('^', unique(STATEFP)), GEOID)), by = c('GEOID10'='GEOID')))) %>%
  ungroup %>%
  mutate(data = case_when(
    STATEFP == 41 ~ map(.x = data, .f = ~left_join(.x, ord_sf_tr_summary, by = c('GEOID10' = 'GEOID10', 'year' = 'dth_year'))),
    STATEFP == 53 ~ map(.x = data, .f = ~left_join(.x, wad_sf_tr_summary, by = c('GEOID10' = 'GEOID10', 'year' = 'dth_year')))
  )) %>%
  unnest(cols = 'data') %>%
  ungroup %>%
  filter(year %in% 2006:2018) %>%
  mutate(State = ifelse(STATEFP == 41, "Oregon", 'Washington')) %>%
  arrange(year, GEOID10)



## Nursing home (output: nursing_join) ####
nursing = st_read(str_c(datapath, "/Data/NursingHomes_BusinessAnalyst.geojson")) %>%
  st_transform(5070)
cenpop = read.csv(str_c(datapath, "/Data/TotalPopulation_2010.csv")) %>%
    transmute(GEOID10 = str_sub(GEO_ID, 10, 14),
              pop_ref = as.integer(P001001)) %>%
    group_by(GEOID10) %>%
    summarize(pop_ref = sum(pop_ref, na.rm = TRUE)) %>%
    ungroup

nursing_join = nursing %>%
    filter(STATE %in% c('OR', 'WA') & STATUS != "T") %>%
    st_join(cnty_orwa) %>%
    st_drop_geometry %>%
    group_by(GEOID10) %>%
    summarize(n_nursing = n(),
              n_nursing_emp = sum(EMPNUM)) %>%
    ungroup %>%
    full_join(cenpop) %>%
    mutate_at(.vars = vars(starts_with('n_nursing')),
              .funs = list(~ifelse(is.na(.), 0, .))) %>%
    transmute(GEOID10 = GEOID10, 
              d_nursing_10k = 1e4 * n_nursing / pop_ref,
              d_nursing_emp_10k = 1e4 * n_nursing_emp / pop_ref)


## Waterbody and greenspace
# Waterbody (output: orwa_cnty_0618_nws) ####
orwa_tracts_0618_nw = st_read(str_c(datapath, "/Data/Water/ORWATracts2010.gpkg"), layer = "Tracts_WaterErased")
orwa_tracts_0618_ns = st_read(str_c(datapath, "/Data/Water/ORWATracts2010.gpkg"), layer = "Tracts_OceanSubtracted")
orwa_tracts_0618 = st_read(str_c(datapath, "/Data/Water/ORWATracts2010.gpkg"), layer = "ORWATracts2010")

library(units)
## area(ns) - area(nw) = area(freshwater)
orwa_cnty_0618_nw = orwa_tracts_0618_nw %>%
    mutate(GEOID10_cnty = str_sub(GEOID10, 1, 5)) %>%
    mutate(area_nonwater = set_units(st_area(geom), 'mi^2')) %>%
    st_drop_geometry %>%
    group_by(GEOID10_cnty) %>%
    summarize(area_nonwater = as.numeric(sum(area_nonwater, na.rm = TRUE))) %>%
    ungroup
orwa_cnty_0618_ns = orwa_tracts_0618_ns %>%
    mutate(GEOID10_cnty = str_sub(GEOID10, 1, 5)) %>%
    mutate(area_nonwater = set_units(st_area(geom), 'mi^2')) %>%
    st_drop_geometry %>%
    group_by(GEOID10_cnty) %>%
    summarize(area_nonsea = as.numeric(sum(area_nonwater, na.rm = TRUE))) %>%
    ungroup
orwa_cnty_0618_nws = orwa_cnty_0618_nw %>%
    left_join(orwa_cnty_0618_ns) %>%
    mutate(p_water_pct_strict = as.numeric(100 * (area_nonsea - area_nonwater)/area_nonsea)) %>%
    mutate(p_water_pct_strict = if_else(p_water_pct_strict < 0, 0, p_water_pct_strict)) %>%
    rename(GEOID10 = GEOID10_cnty)


### Objective: to calculate the proportion of parks
parks = st_read(str_c(datapath, "/Data/Parks/ORWA_Parks_ParkServe_Parks_Playgrounds_020722.shp")) %>%
  st_transform(5070)

# Nonprivate (output: orwa_cnty_pks) ####
# parks_proj = parks %>% st_transform(st_crs(orwa_tracts))
parks_proj = parks %>% 
    filter(!Park_Manag %in% c('PRIV')) %>%
    dplyr::select(ParkID)
orwa_cnty_parks = st_intersection(cnty_orwa, parks_proj)
orwa_cnty_parksd = orwa_cnty_parks %>%
    filter(!is.na(ParkID)) %>%
    mutate(park_nonpriv_area_km2 = as.numeric(units::set_units(st_area(geometry), "km^2"))) %>%
    st_drop_geometry %>%
    group_by(GEOID10) %>%
    summarize(n_parks_nonpriv = n(),
              park_nonpriv_area_km2 = sum(park_nonpriv_area_km2, na.rm = T)) %>%
    ungroup
orwa_cnty_pks = cnty_orwa %>% 
    #dplyr::select(-all_of(colnames(orwa_tracts_parksd)[-1])) %>%
    left_join(orwa_cnty_parksd) %>%
    transmute(GEOID10 = GEOID10,
              park_nonpriv_area_km2 = ifelse(is.na(park_nonpriv_area_km2), 0, park_nonpriv_area_km2),
              n_parks_nonpriv = ifelse(is.na(n_parks_nonpriv), 0, n_parks_nonpriv),
              #  d_parks_nonpriv_10k = 1e4 * n_parks_nonpriv / n_pop_total,
              d_parks_nonpriv_km2 = n_parks_nonpriv / as.numeric(units::set_units(st_area(geometry), "km^2")),
              p_park_nonpriv = 100* park_nonpriv_area_km2 / as.numeric(units::set_units(st_area(geometry), "km^2"))) %>%
    st_drop_geometry


## SES composite (output: orwa_cnty_pca$ind$coord[,1])
orwa_cnty_ess = orwa_cnty_0618df %>%
    mutate(n_medincome_10k = n_medincome_adj / 1e4) %>%
    dplyr::select(n_medincome_10k,
                  p_poverty, p_unemp, p_edubac,
                  year,
                  State) %>%
    st_drop_geometry

orwa_cnty_pca = FactoMineR::PCA(orwa_cnty_ess, ncp = 3, quanti.sup = 5, quali.sup = 6, graph= FALSE)
plot(orwa_cnty_pca, choix = "var", label = c("var", "quali"), graph.type = "ggplot", col.quali = 'red')

# write.csv(orwa_cnty_pca$var$cor, "/mnt/c/Users/isong/OneDrive/UO/Paper/MentalMortality_Tract/SESComposite_loading_061922.csv")
# write.csv(orwa_cnty_pca$eig, "/mnt/c/Users/isong/OneDrive/UO/Paper/MentalMortality_Tract/SESComposite_eigenvalues_061922.csv")

# orwa_cnty_0618_ses = orwa_cnty_0618_covars %>%
#     dplyr::select(-starts_with('SES')) %>%
#     mutate(SES = orwa_cnty_pca$ind$coord[,1])



## PM2.5 (output: pm25_vec)
mean_extract = function(ras, pol, foo = 'mean') {
    orwa_tracts_evi = exactextractr::exact_extract(ras, pol, foo, full_colnames = TRUE)
    return(orwa_tracts_evi)
}

pm25_list = list.files("/media/felix/Processing/PM25-selected/",
                       pattern = "*.asc$",
                       full.names = TRUE)
pm25_list_l = pm25_list[-length(pm25_list)] %>%
    split(.,.) %>%
    lapply(function(x) mean_extract(terra::rast(x), cnty_orwa %>% arrange(GEOID10) %>% st_transform(4326)))
pm25_vec = Reduce(c, pm25_list_l)


## Join all
orwa_cnty_0618_covars = orwa_cnty_0618df %>%
  left_join(nursing_join) %>%
  left_join(orwa_cnty_0618_nws) %>%
  left_join(orwa_cnty_pks) %>%
  mutate(SES = orwa_cnty_pca$ind$coord[,1],
         pm25_mean = pm25_vec) %>%
  st_as_sf(sf_column_name = 'geometry')

write_rds(orwa_cnty_0618_covars, "./pnw_reanalysis/ORWA_county_0618_covariates_122622.rds",
          compress = 'xz')


form_base_spt <- outcome ~ 
    offset(log(tpop)) +
    f(GEOID_int1, model = 'bym2', graph = graph, constr = TRUE, scale.model = scalemodel) +
    f(year_int1, model = 'ar1', constr = TRUE) +
    State +
    green + 
    pm25_mean +
    log(p_water_pct_strict+ exp(-25)) +
    d_nursing_10k +
    SES + 
    p_elderly + p_singleheaded + p_bmarital + log(p_nonwhite + 0.1)





### End of file ####