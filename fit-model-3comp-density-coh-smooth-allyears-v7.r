#### SETUP ####
source("Scripts/setup.R")

#### IMPORT LIFE TABLE CSV ####
source("Scripts/import_mortality.R")


#### IMPORT LIFE TABLE DATA AS RDS OBJECTS ####

prepare_data<-function(l, y, cod, min.age, max.age, s, years){
  m_r<-dplyr::filter(l[years==y][[1]], codreg==cod, sex==s)
  m_data=list(N=with(dplyr::filter(m_r, age>=min.age, age<=max.age),
                     sum(dx)),
              x1=min.age,
              x2=max.age,
              dx=with(dplyr::filter(m_r, age>=min.age, age<=max.age),
                      dx)
  )
  m_data
}

m_data<-readRDS("Life table data FVG.rds")
m_data_lombardia<-readRDS("Life table data Lombardia.rds")
m_data_lazio<-readRDS("Life table data Lazio.rds")
m_data_sicilia<-readRDS("Life table data Sicilia.rds")
m_data_sardegna<-readRDS("Life table data Sardegna.rds")


#### FIT ####

m_3comp_sm_coh_ay <- stan_model('single_region_mixture_Dx all_years density-like smooth cohort effect v7.stan', model_name="Single region density-like smooth cohort effect all years")

data_3comp=m_data_lombardia 
years_3comp=m_years
name_3comp="lombardia"
# data_3comp=m_data_sicilia
# years_3comp=m_years
# name_3comp="sicilia"
# data_3comp=m_data_sardegna 
# years_3comp=m_years
# name_3comp="sardegna"
# data_3comp=m_data_lazio
# years_3comp=m_years
# name_3comp="lazio"
# data_3comp=m_data #slightly less than 24 hours
# years_3comp=m_years
# name_3comp="fvg"

cohort_years_model=c(1915:1925)
coh_effect<-list(N_coh=length(cohort_years_model),
                 cohorts=cohort_years_model #etc., all cohorts for which the relative cohort effect exceeds a threshold (or even all!)
)
#prepare data
region_data_allyears<-list(x1=data_3comp[[1]]$x1,
                           x2=data_3comp[[1]]$x2,
                           t1=years_3comp[1],
                           t2=years_3comp[length(years_3comp)],
                           N=sapply(data_3comp, function(x) x$N),
                           dx=t(sapply(data_3comp, function(x) x$dx)))

init_args_3comp_c_z<-function(){
  zeta3<-runif(1, 0.6, 0.8)
  zeta2<-runif(1, 0.2, 0.3)
  zeta1<-runif(1, 0.05, 0.15)
  zeta<-c(zeta1/sum(zeta1, zeta2, zeta3), zeta2/sum(zeta1, zeta2, zeta3), zeta3/sum(zeta1, zeta2, zeta3))
  mu=c(runif(1, 20, 35), runif(1, 75, 85))
  sigma_i=runif(1, 0.6, 0.8)
  sigma_m=runif(1, 5, 12)
  gamma_m=runif(1, 0.15, 0.5)
  sigma_M=runif(1, 5, 12)
  gamma_M=runif(1, -0.9, -0.15)
  list(zeta=zeta, sigma_i=sigma_i, mu=mu, 
       sigma_m=sigma_m, gamma_m=gamma_m,
       sigma_M=sigma_M, gamma_M=gamma_M)
  
}

init_args_3comp_c_z_ay<-function(n, n_alpha){
  res<-lapply(1:n, function(x) init_args_3comp_c_z())
  zeta<-t(sapply(res, function(x) x$zeta))
  mu<-t(sapply(res, function(x) x$mu))
  sigma_m<-sapply(res, function(x) x$sigma_m)
  gamma_m<-sapply(res, function(x) x$gamma_m)
  sigma_M<-sapply(res, function(x) x$sigma_M)
  gamma_M<-sapply(res, function(x) x$gamma_M)
  sigma_i<-runif(n, 0.8, 1.2)
  log_alpha_raw=runif(n_alpha, -0.22, 0.18) #multiplier between 0.80 and 1.20
  return(list(zeta=zeta, mu=mu, sigma_i=sigma_i, sigma_m=sigma_m, gamma_m=gamma_m,
              sigma_M=sigma_M, gamma_M=gamma_M, log_alpha_raw=log_alpha_raw))
}


print(paste0("beginning ", name_3comp, ", ", now()))
#fit model
fit_3comp_sm_coh_ay<- sampling(m_3comp_sm_coh_ay, 
                               data = c(region_data_allyears, coh_effect),
                               refresh = 10, seed=1,
                               init=function() {init_args_3comp_c_z_ay(length(years_3comp), length(cohort_years_model))},
                               iter=2000,
                               warmup=1000,
                               chains=4,
                               cores=2)

print(paste0("finished, ", now()))


write_rds(fit_3comp_sm_coh_ay, paste0("3component centered ", name_3comp, " density model ",years_3comp[1], "-", years_3comp[length(years_3comp)], " smooth cohort effect all years v7 cauchy ", format(Sys.time(), format="%Y%m%d"), ".rds"))

#output to screen
print(fit_3comp_sm_coh_ay, pars = c(paste0('zeta[', 1:length(years_3comp), ',1]'), paste0('mu[', 1:length(years_3comp), ',1]')),
      probs = c(0.1, 0.5, 0.9))
print(fit_3comp_sm_coh_ay, c(paste0('sigma_m[', 1:length(years_3comp), ']'), paste0('gamma_m[', 1:length(years_3comp), ']')),
      probs = c(0.1, 0.5, 0.9))
print(fit_3comp_sm_coh_ay, pars = c(paste0('zeta[', 1:length(years_3comp), ',2]'), paste0('mu[', 1:length(years_3comp), ',2]')),
      probs = c(0.1, 0.5, 0.9))
print(fit_3comp_sm_coh_ay, pars = c(paste0('xi_M[', 1:length(years_3comp), ']'), paste0('sigma_M[', 1:length(years_3comp), ']')),
      probs = c(0.1, 0.5, 0.9))
print(fit_3comp_sm_coh_ay, pars = c(paste0('gamma_M[', 1:length(years_3comp), ']'), paste0('zeta[', 1:length(years_3comp), ',3]')),
      probs = c(0.1, 0.5, 0.9))
print(fit_3comp_sm_coh_ay, pars = c(paste0('log_alpha[', 1:length(cohort_years_model), ']')),
      probs = c(0.1, 0.5, 0.9))
print(fit_3comp_sm_coh_ay, pars = c("denominator"),
      probs = c(0.1, 0.5, 0.9))

#csv
par_neff_tseries<-lapply(1:length(years_3comp), function(i) {
  res<-data.frame(year=years_3comp[i],
                  as.data.frame(t(summary(fit_3comp_sm_coh_ay)$summary[c(paste0('zeta[',i,',1]'), paste0('zeta[',i,',2]'), paste0('mu[',i,',1]'), paste0('mu[',i,',2]'), paste0('sigma_m[',i,']'), paste0('gamma_m[',i,']'), paste0('xi_m[',i,']'), paste0('xi_M[',i,']'), paste0('sigma_M[',i,']'), paste0('gamma_M[',i,']'), paste0("log_alpha[", 1:length(cohort_years_model),"]"), paste0('sigma_i[',i,']')),"n_eff"])))
  names(res)<-c('year', 'zeta[1]', 'zeta[2]', 'mu[1]', 'mu[2]', 'sigma_m', 'gamma_m', 'xi_m', 'xi_M', 'sigma_M', 'gamma_M', paste0("alpha[", 1:length(cohort_years_model),"]"), 'sigma_i')
  res
}
)
par_neff_tseries<-as.data.frame(do.call(rbind, par_neff_tseries))

par_median_tseries<-lapply(1:length(years_3comp), function(i) {
  res<-data.frame(year=years_3comp[i],
                  as.data.frame(t(summary(fit_3comp_sm_coh_ay)$summary[c(paste0('zeta[',i,',1]'), paste0('zeta[',i,',2]'), paste0('mu[',i,',1]'), paste0('mu[',i,',2]'), paste0('sigma_m[',i,']'), paste0('gamma_m[',i,']'), paste0('xi_m[',i,']'), paste0('xi_M[',i,']'), paste0('sigma_M[',i,']'), paste0('gamma_M[',i,']'), paste0("log_alpha[", 1:length(cohort_years_model),"]"), paste0('sigma_i[',i,']')),"50%"])))
  names(res)<-c('year', 'zeta[1]', 'zeta[2]', 'mu[1]', 'mu[2]', 'sigma_m', 'gamma_m', 'xi_m', 'xi_M', 'sigma_M', 'gamma_M', paste0("log_alpha[", 1:length(cohort_years_model),"]"), 'sigma_i')
  res
}
)
par_median_tseries<-as.data.frame(do.call(rbind, par_median_tseries))

par_Rhat_tseries<-lapply(1:length(years_3comp), function(i) {
  res<-data.frame(year=years_3comp[i],
                  as.data.frame(t(summary(fit_3comp_sm_coh_ay)$summary[c(paste0('zeta[',i,',1]'), paste0('zeta[',i,',2]'), paste0('mu[',i,',1]'), paste0('mu[',i,',2]'), paste0('sigma_m[',i,']'), paste0('gamma_m[',i,']'), paste0('xi_m[',i,']'), paste0('xi_M[',i,']'), paste0('sigma_M[',i,']'), paste0('gamma_M[',i,']'), paste0("log_alpha[", 1:length(cohort_years_model),"]"), paste0('sigma_i[',i,']')), "Rhat"])))
  names(res)<-c('year', 'zeta[1]', 'zeta[2]', 'mu[1]', 'mu[2]', 'sigma_m', 'gamma_m', 'xi_m', 'xi_M', 'sigma_M', 'gamma_M', paste0("log_alpha[", 1:length(cohort_years_model),"]"), 'sigma_i')
  res
}
)
par_Rhat_tseries<-as.data.frame(do.call(rbind, par_Rhat_tseries))

par_10pc_tseries<-lapply(1:length(years_3comp), function(i) {
  res<-data.frame(year=years_3comp[i],
                  as.data.frame(t(summary(fit_3comp_sm_coh_ay, probs = c(0.1, 0.5, 0.9))$summary[c(paste0('zeta[',i,',1]'), paste0('zeta[',i,',2]'), paste0('mu[',i,',1]'), paste0('mu[',i,',2]'), paste0('sigma_m[',i,']'), paste0('gamma_m[',i,']'), paste0('xi_m[',i,']'), paste0('xi_M[',i,']'), paste0('sigma_M[',i,']'), paste0('gamma_M[',i,']'), paste0("log_alpha[", 1:length(cohort_years_model),"]"), paste0('sigma_i[',i,']')), "10%"])))
  names(res)<-c('year', 'zeta[1]', 'zeta[2]', 'mu[1]', 'mu[2]', 'sigma_m', 'gamma_m', 'xi_m', 'xi_M', 'sigma_M', 'gamma_M', paste0("log_alpha[", 1:length(cohort_years_model),"]"), 'sigma_i')
  res
}
)
par_10pc_tseries<-as.data.frame(do.call(rbind, par_10pc_tseries))

par_90pc_tseries<-lapply(1:length(years_3comp), function(i) {
  res<-data.frame(year=years_3comp[i],
                  as.data.frame(t(summary(fit_3comp_sm_coh_ay, probs = c(0.1, 0.5, 0.9))$summary[c(paste0('zeta[',i,',1]'), paste0('zeta[',i,',2]'), paste0('mu[',i,',1]'), paste0('mu[',i,',2]'), paste0('sigma_m[',i,']'), paste0('gamma_m[',i,']'), paste0('xi_m[',i,']'), paste0('xi_M[',i,']'), paste0('sigma_M[',i,']'), paste0('gamma_M[',i,']'), paste0("log_alpha[", 1:length(cohort_years_model),"]"), paste0('sigma_i[',i,']')), "90%"])))
  names(res)<-c('year', 'zeta[1]', 'zeta[2]', 'mu[1]', 'mu[2]', 'sigma_m', 'gamma_m', 'xi_m', 'xi_M', 'sigma_M', 'gamma_M', paste0("log_alpha[", 1:length(cohort_years_model),"]"), 'sigma_i')
  res
}
)
par_90pc_tseries<-as.data.frame(do.call(rbind, par_90pc_tseries))

write.csv2(par_neff_tseries, paste0("parameters n_eff ", name_3comp, " ", years_3comp[1], "-", years_3comp[length(years_3comp)]," 3comp density-like smooth cohort effect all years v7 ", format(Sys.time(), format="%Y%m%d"),".csv"), row.names=FALSE)
write.csv2(par_median_tseries, paste0("parameters median ", name_3comp, " ", years_3comp[1], "-", years_3comp[length(years_3comp)]," 3comp density-like smooth cohort effect all years v7 ", format(Sys.time(), format="%Y%m%d"),".csv"), row.names=FALSE)
write.csv2(par_Rhat_tseries, paste0("parameters Rhat ", name_3comp, " ", years_3comp[1], "-", years_3comp[length(years_3comp)]," 3comp density-like smooth cohort effect all years v7 ", format(Sys.time(), format="%Y%m%d"),".csv"), row.names=FALSE)
write.csv2(par_10pc_tseries, paste0("parameters 10th percentile ", name_3comp, " ", years_3comp[1], "-", years_3comp[length(years_3comp)]," 3comp density-like smooth cohort effect all years v7 ", format(Sys.time(), format="%Y%m%d"),".csv"), row.names=FALSE)
write.csv2(par_90pc_tseries, paste0("parameters 90th percentile ", name_3comp, " ", years_3comp[1], "-", years_3comp[length(years_3comp)]," 3comp density-like smooth cohort effect all years v7 ", format(Sys.time(), format="%Y%m%d"),".csv"), row.names=FALSE)



#plots
for(i in 1:length(years_3comp)){
  trace_graph<-mcmc_trace(fit_3comp_sm_coh_ay, pars = c(paste0('zeta[', i, ',1]'), paste0('zeta[', i, ',2]'), paste0('mu[', i, ',1]'), paste0('xi_M[', i, ']'), paste0('lambda_M[', i, ']')))
  png(filename=paste0("trace_plot_", name_3comp, "_3c_coh_sm_ay_v7_",years_3comp[[i]],".png"), width=1680, height=1200)
  print(trace_graph)
  test<-dev.off()

  pairs_graph<-mcmc_pairs(fit_3comp_sm_coh_ay, pars = c(paste0('zeta[', i, ',1]'), paste0('zeta[', i, ',2]'), paste0('mu[', i, ',1]'), paste0('xi_M[', i, ']'), paste0('lambda_M[', i, ']'), paste0('sigma_i[', i, ']')))

  png(filename=paste0("pairs_plot_", name_3comp, "_3c_coh_sm_ay_v7_",years_3comp[[i]],".png"), width=1680, height=1200)
  print(pairs_graph)
  test<-dev.off()

  sims_3comp<-rstan::extract(fit_3comp_sm_coh_ay, pars=paste0('dx_rep[', i, ", ", 1:(max.age-min.age+1), ']'))
  dx_3comp<-as.data.frame(sims_3comp)
  names(dx_3comp)<-paste0("dx", (min.age:max.age))
  dx_3comp$iterations=as.integer(row.names(dx_3comp))
  dx_3comp<-pivot_longer(dx_3comp, cols=!iterations, names_to="age", names_prefix="dx", names_transform=as.integer, values_to="dx")
  p_dx_3comp<-ggplot() +
    geom_line(aes(x=data_3comp[[i]]$x1:data_3comp[[i]]$x2, y=data_3comp[[i]]$dx), size=2) +
    geom_line(data=filter(dx_3comp, iterations<=200), aes(x=age, y=dx, group=iterations), size=0.5, color="lightblue", alpha=0.2) +
    labs(x="Age", y="Number of deaths", title=paste0("Deaths by age, 3-component model with centered param. and smooth cohort effect, all years estimated together, year ", years_3comp[[i]])) +
    theme_minimal()

  png(filename=paste0("line_overlay_", name_3comp, "_3c_coh_sm_ay_v7_",years_3comp[[i]],".png"), width=1680, height=1200)
  print(p_dx_3comp)
  test<-dev.off()
}

#zip files with graphs
zip(zipfile=paste0("trace plots ", name_3comp, " ", years_3comp[1], "-", years_3comp[length(years_3comp)]," 3comp density-like smooth cohort effect all years v7 ", format(Sys.time(), format="%Y%m%d"),".zip"), files=list.files(pattern=paste0("trace_plot_", name_3comp, "_3c_coh_sm_ay_v7_")))
zip(zipfile=paste0("pairs plots ", name_3comp, " ", years_3comp[1], "-", years_3comp[length(years_3comp)]," 3comp density-like smooth cohort effect all years v7 ", format(Sys.time(), format="%Y%m%d"),".zip"), files=list.files(pattern=paste0("pairs_plot_", name_3comp, "_3c_coh_sm_ay_v7_")))
zip(zipfile=paste0("line overlays ", name_3comp, " ", years_3comp[1], "-", years_3comp[length(years_3comp)]," 3comp density-like smooth cohort effect all years v7 ", format(Sys.time(), format="%Y%m%d"),".zip"), files=list.files(pattern=paste0("line_overlay_", name_3comp, "_3c_coh_sm_ay_v7_")))

#residuals
res_3comp_ay<-mapply(function(i, d, m){
  #d=data, m=model
  #take the median of dx_rep from m
  d$dx - summary(m)$summary[paste0("dx_rep[",i, ",", min.age:max.age+1, "]"),"50%"]
}, 1:length(years_3comp), data_3comp, MoreArgs=list(m=fit_3comp_sm_coh_ay), SIMPLIFY=TRUE)

res_3comp_ay_df<-pivot_longer(data.frame(as.data.frame(res_3comp_ay), age=min.age:max.age), cols=!age, names_prefix="V", names_to="year", values_to="residual", names_transform=list(year=function(x) {as.integer(x)+m_years[1]-1 }))

p1001=ggplot(data=res_3comp_ay_df, mapping=aes(x=age, y=year, z=residual)) +
  geom_contour_filled() + theme_minimal() +
  labs(title="Raw residuals by age and year, 3 comp. model with cohort effects calculated on all years together", x="Age", y="Year", fill="Actual deaths")

#output residuals plot
png(paste0("res_abs_", name_3comp, "_3c_coh_sm_ay_v7.png"), width=1680, height=1200)
print(p1001)
test=dev.off()

#plot of medians, facet_wrap
#par_years_conv<-apply(par_neff_tseries, 1, function(x) min(x[c(2:5, 7:13, 15)]))
par_neff_df<-par_neff_tseries %>% select(1:11, ncol(par_neff_tseries)) %>%
  pivot_longer(cols=!year, names_to="parameter", values_to="neff")

par_median_df<-par_median_tseries %>% select(1:11, ncol(par_median_tseries)) %>%
  pivot_longer(cols=!year, names_to="parameter", values_to="value") %>%
  left_join(y=par_neff_df, by=c("year", "parameter"))
par_plot<-ggplot(par_median_df, aes(x=year, y=value)) +
  geom_line() +
  theme_minimal() +
  facet_wrap("parameter", scales="free_y")

par_plot2<-ggplot(filter(par_median_df, neff>100), aes(x=year, y=value)) +
  geom_line() +
  theme_minimal() +
  facet_wrap("parameter", scales="free_y")

png(paste0("par_plots_", name_3comp, "_3c_coh_sm_ay_v7.png"), width=1680, height=1200)
print(par_plot)
test=dev.off()

png(paste0("par_plots_", name_3comp, "_3c_coh_sm_ay_v7_filtered.png"), width=1680, height=1200)
print(par_plot2)
test=dev.off()