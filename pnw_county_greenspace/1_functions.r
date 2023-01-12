## Base file
## Main function ####
fit_model = function(formula, data, greenvar = "EVI_1y", outcome = 'n_nonsubstance', lincomb, thetas = NULL, 
                     graph = graph, family = 'poisson', threads = "12:1",
                     mode = "plain") {
    data_norm = data %>%
        mutate(green = !!sym(greenvar)) %>%
        group_by(year) %>%
        mutate(Expected = !!sym(outcome) * (sum(!!sym(outcome))/sum(n_pop_total)),
               outcome = !!sym(outcome))
    #print(fivenum(data_norm$Expected))
    inla(
        formula = formula,
            data = data_norm,
            family = family, verbose = TRUE,
            lincomb = lincomb,
            #E = ifelse(mode != 'plain', Expected, NULL),
            control.lincomb = list(verbose = FALSE),
            control.inla = control.inla(#b.strategy = 'keep', 
                                        strategy = "simplified.laplace", 
                                        int.strategy = 'eb',
                                        restart = 5,
                                        optimiser = 'default',
                                        stupid.search.factor = 1.05,
                                        #tolerance = 1e-6, 
                                        npoints = 9,
                                        #dz = 0.75,
                                        #cmin = -Inf,
                                        improved.simplified.laplace = TRUE,
                                        use.directions = TRUE),
            control.predictor = list(compute = TRUE, link = 1), #
            control.compute = list(cpo = TRUE, waic = TRUE, dic = TRUE,
                                   return.marginals.predictor = TRUE,
                                   openmp.strategy = 'huge',
                                   smtp = 'pardiso'),
            control.mode = list(restart = TRUE, fixed = FALSE,
                                #x = 0,
                                theta = thetas),
            num.threads = threads,
            safe = TRUE)
}


## Refit until getting the mode.status==0
repeat_refit = function(inla_result, criterion = 'mode', crit.value = 0) {
    i = 1
    stopifnot(criterion %in% c('mode', 'dic', 'hessian'))
    if (criterion == 'mode') {
        while (inla_result$mode$mode.status != crit.value) {
            cat(str_c("Refit INLA model (criterion: mode.status), iteration: ", i, " \n"))
            inla_result = inla.rerun(inla_result)
            i = i + 1
        }
    } else if (criterion == "dic") {
        while (inla_result$dic$dic < crit.value) {
            cat(str_c("Refit INLA model (criterion: DIC), iteration: ", i, " \n"))
            inla_result = inla.rerun(inla_result)
            i = i + 1
        }
    } else if (criterion == "hessian") {
        while (sum(duplicated(inla_result$misc$cov.intern.eigenvalues)) != 0) {
            cat(str_c("Refit INLA model (criterion: DIC), iteration: ", i, " \n"))
            inla_result = inla.rerun(inla_result)
            i = i + 1
        }
    }
    cat(str_c("Finish refitting INLA models (total iteration: ", i, " times)\n"))
    return(inla_result)
}

## Effect classification and filtering
detect_effect = function(inla_coef, N_data, effect_id, add = FALSE) {
    if (sum(grepl("(GEOID|year)", effect_id) == 0)) {
        stop ("Please check your effect id. They are supposed to include GEOID or year.")
    }
    ncoef = nrow(inla_coef)
    if (N_data < ncoef) {
        if (add) {
            coefs = inla_coef[seq.int(1, N_data), "mean"] + inla_coef[seq.int(N_data + 1, 2* N_data), "mean"]
        } else {
            coefs = inla_coef[seq.int(1, N_data), "mean"]
        }
        
    } else {
        coefs = inla_coef[,"mean"]
    }
    coefs = unlist(coefs, use.names = FALSE)
    return(coefs)
}



## Mapping mode
vis_spt_effect = function(map = orwa_tracts_0618_mdpp,
                          inla_fit,
                          map_state = NULL,
                          time_col = "year",
                          sp_col = "GEOID10",
                          draw_col = "beta_1it",
                          effect_id = NULL,
                          outcome = "n_nonsubstance",
                          mode = "effect",
                          title = expression(beta, '[1i]'),
                          file_export = FALSE,
                          filepath = NULL) {
        
        if (!is.null(effect_id)) {
            inla_fit_re = inla_fit$summary.random[[effect_id]]
        } else {
            inla_fit_re = inla_fit$summary.lincomb.derived
        }

        N = nrow(unique(st_drop_geometry(map[,sp_col])))
        if (nrow(inla_fit_re) < N) {
            stop ("It seems like you want to visualize temporal effects.")
        } else if (nrow(inla_fit_re) == N | nrow(inla_fit_re) == 2 * N) {
            map = map %>%
                filter(year == 2010)
        } else {
            map = map
        }

        coefs = detect_effect(inla_fit_re, nrow(map), effect_id)

        if (mode == "effect") {
            map_ext = map %>%
                mutate(beta_1it = coefs)
        } else if (mode == "residuals") {
            map_ext = map %>%
                mutate(prediction = inla_fit$summary.fitted.values[,1] %>% unlist,
                    #delta_1it_park = mod_moodanxiety_orwa_inla_pois_aq1$summary.random[[2]][,'mean'] %>% unlist,
                    beta_1it = coefs) %>%
                mutate(resid_pred = !!sym(outcome) - prediction)
        } else {
            stop("The argument mode should be one of 'effect' or 'residuals.'")
        }

        if (is.null(effect_id)) {
            tm_inla = tm_shape(map_ext) +
                #tm_style('classic') +
                tm_fill(draw_col, pal = '-RdBu', n = 6, style = 'cont', 
                        aes.palette = 'seq', midpoint = 0,
                        title = title) +
                tm_borders(lwd = 0.15, col = 'dark grey') +
                tm_layout(inner.margins = 0.01, outer.margins = c(0.01, 0.01, 0.01, 0.01), frame = FALSE)
            if (length(coef) > N) {
                tm_inla = tm_inla +
                    tm_facets(time_col) +
                    tm_layout(inner.margins = 0.01, outer.margins = c(0.01, 0.01, 0.01, -0.12))
            }
        } else {
            tm_inla = tm_shape(map_ext) +
                #tm_style('classic') +
                tm_fill(draw_col, pal = '-RdBu', n = 6, style = 'cont', 
                        aes.palette = 'seq', midpoint = 0,
                        title = title) +
                tm_borders(lwd = 0.15, col = 'dark grey') +
                tm_layout(inner.margins = 0.01, outer.margins = c(0.01, 0.01, 0.01, 0.01), frame = FALSE)
        }

        if (!is.null(map_state)) {
            tm_inla = tm_inla +
                tm_shape(map_state) +
                tm_borders(col = 'black', lwd = 0.5)
        }

        tm_inla

        if (file_export) {
            tmap_save(tm_inla,
                      filename = filepath,
                      width = 10, height = 12.5,
                      #outer.margins = 0.02,
                      units = "in",
                      dpi = 508, scale = 1.33)
        }
                          }


## 
inla_predacc = function(inla_result, data, target_var, pop_var = NULL, acc_measure) {
    basicdiff = function(x, xhat) { x - xhat}
    mape = function(x, xhat) {
        100 * mean ( abs( basicdiff(x, xhat) / x ) )
    }
    mse = function(x, xhat) {
        mean(basicdiff(x, xhat) ^ 2)
    }
    rmse = function(x, xhat) {
        sqrt(mse(x, xhat))
    }
    mae = function(x, xhat) {
        mean(abs(basicdiff(x, xhat)))
    }
    data = st_drop_geometry(data)

    truth = unlist(data[,target_var])
    preded = inla_result$summary.fitted.values$mean * unlist(data[,pop_var])

    acc = 
    switch(acc_measure,
            mape = mape(truth, preded),
            mse = mse(truth, preded),
            rmse = rmse(truth, preded))
    return(acc)
}





### Spatiotemporal interaction based on generic definitions ####
## Ref: https://github.com/spatialstatisticsupna/Mmodels_SERRA_article/blob/master/R/functions/MCAR_INLA_st.R

# data: data frame
# graph: matrix or Matrix. Including binary connectivity
# formula
# ID.main: to identify the identifier for the spt structure matrix
# ID.area: the area identifier
# ID.year: the time identifier
# prior.interaction: 1-4. the type of spatiotemporal interaction
generate_spt_prec_formula = function(data, graph, formula, ID.main = 'GEOID_year1', ID.area = 'GEOID10', ID.year = 'year', prior.interaction = 1) {
    if (any(grepl('sf', class(data)))) {
    data = as.data.frame(st_drop_geometry(data))
    }
    n <- length(unique(data[,ID.area]))
    t <- length(unique(data[,ID.year]))

    ## Adjacency matrices                   ##
    ## ---------------------------------------
    ## W: adjacency matrix (spatial)
    Ws <- graph
                    
    ## W.t: adjacency matrix (temporal)
    Pt <- crossprod(diff(diag(t),differences=1))
    Wt <- INLA::inla.as.sparse(Matrix::Diagonal(n=nrow(Pt),diag(Pt))- Pt) ## temp.corre = TRUE

    ## Precision matrices                   ##
    ## ---------------------------------------
    Rs <- Matrix::Diagonal(n,colSums(Ws))-Ws
    Rt <- Matrix::Diagonal(t,colSums(Wt))-Wt

    # R.st: Precision matrices for spatio-temporal interaction
    if(prior.interaction %in% c(1)){
        R.st <- diag(n*t)
        r.def.st <- 0
    }
    if(prior.interaction %in% c(2)){
        R.st <- kronecker(Rt, Matrix::Diagonal(n,1))
        r.def.st <- n
    }
    if(prior.interaction %in% c(3)){
        R.st <- kronecker(Matrix::Diagonal(t,1),Rs)
        r.def.st <- t
    }
    if(prior.interaction %in% c(4)){
        R.st <- kronecker(Rt, Rs)
        r.def.st <- n+t-1
    }
    assign('R.st', R.st, envir = .GlobalEnv)
    assign('r.def.st', r.def.st, envir = .GlobalEnv)

    ## Define appropriate constraints matrices ##
    ## ---------------------------------------
    A.constr.s<- kronecker(diag(1), matrix(1,1,n))
    A.constr.t<- kronecker(diag(1), matrix(1,1,t))

    if(prior.interaction %in% c(1)){
        A.constr.st <- matrix(1, nrow=1, ncol=n*t)
    }
    if(prior.interaction %in% c(2)){
        A.constr.st <- kronecker(matrix(1,nrow=1, ncol=t), diag(1,n))
    }
    if(prior.interaction %in% c(3)){
        A.constr.st <- kronecker(diag(1,t), matrix(1,nrow=1, ncol=n))
    }
    if(prior.interaction %in% c(4)){
        A1 <- kronecker(matrix(1,nrow=1, ncol=t), diag(1,n))
        A2 <- kronecker(diag(1,t), matrix(1,nrow=1, ncol=n))
        A.constr.st <- rbind(A1[-1,], A2)
    }
    assign('A.constr.st', A.constr.st, envir = .GlobalEnv)

    txt_form = as.character(formula)
    form <- as.formula(
        paste(txt_form[[2]], "~", 
              paste(txt_form[[3]], collapse = "+"), "+",
                    paste0("f(", ID.main, ", model='generic0', Cmatrix=R.st, rankdef=r.def.st, constr=FALSE, extraconstr=list(A=A.constr.st, e=rep(0,dim(A.constr.st)[1])))")))
                    #paste0("f(", ID.main, ", model='generic0', Cmatrix=R.st, rankdef=r.def.st, constr=FALSE)")))
    # , hyper=list(prec=list(prior=sdunif))
    return(form)
}



##
# Posterior marginals for computing the probability of random effects larger than a specific value
map_emarginal <- function(idata, map = mort.simp, rand.d = 5, 
                        map_state = NULL,
                        file.export = FALSE,
                        pdir = './output/', 
                        title = expression(paste('p(', beta['1it'], '<0)', sep = '')),
                        filename){
  cat(paste("Calculating p(beta_{1i}<0)....\n"))
  if (is.null(rand.d)){
      emarginal <- lapply(idata$marginals.lincomb.derived,
                      function(x) inla.pmarginal(0, x)) # 1-p (x>0) or p (x<0)
  } else {
      emarginal <- lapply(idata$marginals.random[[rand.d]],#[indx],
                      function(x) inla.pmarginal(0, x))
  }
  emarginal <- do.call(c, emarginal)

  map_rs_all <- map %>%
    dplyr::select(year, GEOID10) %>% 
    mutate(pmarginal = emarginal)
  
  brks = c(0, 0.05, 0.2, 0.5, 0.8, 0.95, 1)
  tms <- tm_shape(map_rs_all) +
    tm_borders(col = 'light grey', lwd = 0.08) +
    tm_fill('pmarginal', palette = 'Blues', breaks=brks,
            title = title) +
    tm_layout(frame = FALSE, frame.lwd = NA, panel.label.bg.color = NA, 
              panel.label.size = 1.5,
              legend.outside = TRUE, legend.outside.position = 'right',
              outer.margins = c(0.01,0.01,0.01,-0.1))
    if (!is.null(map_state)) {
        tms <- tms + 
            tm_shape(map_state) +
            tm_borders(col = 'black', lwd = 0.88)
    }

  if (file.export){
    target_file <- str_c(pdir, filename)
  
    tmap_save(tms, filename = target_file, width = 50, height = 50,
              units = 'cm', dpi = 300, pointsize = 30)

  } else {
      tms
    }
  
}


##
recalculate_ses = function(data,
                           incols = c('n_medincome_10k', 'p_poverty', 'p_edubac', 'p_unemp')) {
    library(FactoMineR)
    data_ses = data %>%
        st_drop_geometry %>%
        dplyr::select(all_of(incols))
    data_ses_pca = FactoMineR::PCA(data_ses, ncp = length(incols), graph = FALSE)
    data = data %>%
        dplyr::select(-starts_with('SES')) %>%
        mutate(SES = data_ses_pca$ind$coord[,1])
    return(data)
                           }
