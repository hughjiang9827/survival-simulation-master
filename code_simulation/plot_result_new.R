load("Moss-simulation-master/code_simulation/df_metric.r")


t_range <- 2:11
df_metric[,se_est:=sqrt(ED2-ED^2)/sqrt(n)]
df_metric[,lower:=estimate-qnorm(0.975)*se_est]
df_metric[,upper:=estimate+qnorm(0.975)*se_est]
df_metric[,coverage:=between(true_surv,lower,upper)]
df_metric[,ci_length:=upper-lower]
perf <- df_metric[t%in%t_range,list(mean_est=mean(estimate),
                                    true_surv=mean(true_surv),
                                    mse = mean((estimate-true_surv)^2),
                                    mse_se = sd((estimate-true_surv)^2)/sqrt(.N),
                                    bias = mean(estimate-true_surv),
                                    bias2 = mean(estimate-true_surv)^2,
                                    var = (.N-1)*var(estimate)/.N,
                                    coverage = mean(coverage),
                                    ci_length = mean(ci_length),
                                    mean_n_t = mean(n_t),
                                    sim_count = .N,
                                    EED = mean(abs(ED)),
                                    EED2 = mean(ED2)),by=list(method,t)]


ggplot(perf,aes(x=t,y=abs(bias),color=method))+geom_line()
ggplot(perf,aes(x=t,y=EED,color=method))+geom_line()
ggplot(perf,aes(x=t,y=mse,color=method))+geom_line()
ggplot(perf,aes(x=t,y=coverage,color=method))+geom_line()
ggplot(perf,aes(x=t,y=ci_length,color=method))+geom_line()
ggplot(perf,aes(x=t,y=mean_est,color=method))+geom_line()+geom_line(aes(y=true_surv),color="black",linetype="dashed")
ggplot(perf,aes(x=t,y=mean_n_t))+geom_line()

long <- melt(perf,id=c("method","t","true_surv"),measure=c("mse","bias2","var"))
ggplot(long,aes(x=t,y=value,color=variable))+geom_line()+facet_wrap(~method)

# also present results in terms of haz
surv_to_haz <- function(x){
  c(-1*diff(x),0)/x
}
df_metric_haz <- copy(df_metric)
df_metric_haz[,estimate:=surv_to_haz(estimate),by=list(method,id_mcmc)]
df_metric_haz[,true_surv:=surv_to_haz(true_surv),by=list(method,id_mcmc)]
perf_haz <- df_metric_haz[t<=10,list(mean_est=mean(estimate),
                                        true_surv=mean(true_surv),
                                        mse = mean((estimate-true_surv)^2),
                                        mse_se = sd((estimate-true_surv)^2)/sqrt(.N),
                                        bias = mean(estimate-true_surv),
                                        bias2 = mean(estimate-true_surv)^2,
                                        var = (.N-1)*var(estimate)/.N),
                      by=list(method,t)]

ggplot(perf_haz,aes(x=t,y=bias,color=method))+geom_line()
ggplot(perf_haz,aes(x=t,y=mse,color=method))+geom_line()
long <- melt(perf_haz,id=c("method","t","true_surv"),measure=c("mse","bias2","var"))
ggplot(long,aes(x=t,y=value,color=variable))+geom_line()+facet_wrap(~method)

