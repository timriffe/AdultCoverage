##########################################################################################################################
ggb_bh<-function(x){           ##  banco
	# head(x)
	colnames(x) <- c("cod","sex","idade","pop1","ano1","pop2","ano2","morte")
tab<-data.frame(x)             ##
dif<-x$ano2[1]-x$ano1[1]
cod<-factor(x$cod)
tab1<-split(tab,cod)
for (i in seq(along = tab1)){
tab1[[i]]$pop1_mais<-(rev(cumsum(rev(tab1[[i]]$pop1))))
tab1[[i]]$pop2_mais<-(rev(cumsum(rev(tab1[[i]]$pop2))))
tab1[[i]]$morte_mais<-(rev(cumsum(rev(tab1[[i]]$morte))))}
for (i in seq(along = tab1)){tab1[[i]]$nascim<-0}
for (i in seq(along = tab1)){
for (j in 4:17) {tab1[[i]]$nascim[j]<-round(0.2*sqrt(tab1[[i]]$pop1[j-1]*tab1[[i]]$pop2[j]))}
        }
for (i in seq(along = tab1)){tab1[[i]]$Lx<-round(sqrt(tab1[[i]]$pop1_mais*tab1[[i]]$pop2_mais))}
for (i in seq(along = tab1)){tab1[[i]]$cresc_mais<-round(log(tab1[[i]]$pop2_mais/tab1[[i]]$pop1_mais)/dif, digits=5)}
for (i in seq(along = tab1)){tab1[[i]]$direito<-round(tab1[[i]]$morte_mais/tab1[[i]]$Lx, digits=5)}
for (i in seq(along = tab1)){tab1[[i]]$esquerdo<-0}
for (i in seq(along = tab1)){
for (j in 4:17) {tab1[[i]]$esquerdo[j]<-round((tab1[[i]]$nascim[j]/tab1[[i]]$Lx[j])-tab1[[i]]$cresc_mais[j],digits=5)}
        }
		
		
		
for (i in seq(along = tab1)){
tab1[[i]]$slope1<-round((sd(tab1[[i]]$esquerdo[4:15])/sd(tab1[[i]]$direito[4:15])), digits=2)                   ## para 10 a 65+
tab1[[i]]$slope2<-round((sd(tab1[[i]]$esquerdo[5:15])/sd(tab1[[i]]$direito[5:15])), digits=2)                   ## para 15 a 65+
tab1[[i]]$slope3<-round((sd(tab1[[i]]$esquerdo[6:15])/sd(tab1[[i]]$direito[6:15])), digits=2)                   ## para 20 a 65+
tab1[[i]]$slope4<-round((sd(tab1[[i]]$esquerdo[7:15])/sd(tab1[[i]]$direito[7:15])), digits=2)                   ## para 25 a 65+
tab1[[i]]$slope5<-round((sd(tab1[[i]]$esquerdo[8:15])/sd(tab1[[i]]$direito[8:15])), digits=2)                   ## para 30 a 65+
tab1[[i]]$intercep<-round(mean(tab1[[i]]$esquerdo[8:15])-tab1[[i]]$slope2*mean(tab1[[i]]$direito[8:15]),digits=4) ##ajuste pelo grupo de 30 a 65+
tab1[[i]]$relcomp<-round(exp(tab1[[i]]$intercep*dif),digits=2)
                          }
tab2<-unsplit(tab1,cod)
relacao<-list(relacao=tab2$relcomp)
relcensos<-unlist(relacao)
return(relcensos)}
### Fun��o para estimar B-H ajustado, segundo a rela��o entre os censos ########
bhadjtable<-function(x){                           ##  banco
tab<-data.frame(x)                                 ##
tab$pop1<-round(tab$pop1/ggb_bh(x),digits=1)
dif<-x$ano2[1]-x$ano1[1]
cod<-factor(x$cod)
tab1<-split(tab,cod)
for (i in seq(along = tab1)){tab1[[i]]$nascim<-0}
for (i in seq(along = tab1)){
for (j in 4:17) {tab1[[i]]$nascim[j]<-round(0.2*sqrt(tab1[[i]]$pop1[j-1]*tab1[[i]]$pop2[j]))}
        }
for (i in seq(along = tab1)){tab1[[i]]$cresc<-round(log(tab1[[i]]$pop2/tab1[[i]]$pop1)/dif, digits=5)}
for (i in seq(along = tab1)){tab1[[i]]$cum_cresc<-round(0, digits=5)
                             tab1[[i]]$cum_cresc[1]<-round(2.5*tab1[[i]]$cresc[1],digits=5)
                            }
for (i in seq(along = tab1)){
for (j in 2:18){tab1[[i]]$cum_cresc[j]<-round(2.5*tab1[[i]]$cresc[j]+5*sum(tab1[[i]]$cresc[(j-1):1]), digits=5)}
        }
for (i in seq(along = tab1)){tab1[[i]]$morte_tab<-round(tab1[[i]]$morte*exp(tab1[[i]]$cum_cresc), digits=0)}
for (i in seq(along = tab1)){tab1[[i]]$razao<-round(sum(tab1[[i]]$morte_tab[4:9])/sum(tab1[[i]]$morte_tab[10:13]), digits=2)}
for (i in seq(along = tab1)){tab1[[i]]$aberto<-round(4.88+((0.872-tab1[[i]]$razao[1])/(0.872-0.827))*(5-4.88), digits=2)} #valores sexo feminino
for (i in seq(along = tab1)){tab1[[i]]$pop_a<-0
tab1[[i]]$pop_a<-round(tab1[[i]]$morte[18]*(exp(tab1[[i]]$aberto[18]*tab1[[i]]$cresc[18])-((tab1[[i]]$aberto*tab1[[i]]$cresc[18])^2/6)),digits=0)
                            }
for (i in seq(along = tab1)){
for(j in 18:1){tab1[[i]]$pop_a[j-1]<-round(tab1[[i]]$pop_a[j]*exp(5*tab1[[i]]$cresc[j-1])+tab1[[i]]$morte[j-1]*exp(2.5*tab1[[i]]$cresc[j-1]),digits=0)}
        }
for (i in seq(along = tab1)){tab1[[i]]$Cx<-0}
for (i in seq(along = tab1)){
for (j in 4:17){tab1[[i]]$Cx[j]<-round(tab1[[i]]$pop_a[j]/tab1[[i]]$nascim[j],digits=2)}
        }
for (i in seq(along = tab1)){
tab1[[i]]$cob1<-round((sum(tab1[[i]]$Cx[4:15])/12), digits=2)                     ## para 10 a 65+
tab1[[i]]$cob2<-round((sum(tab1[[i]]$Cx[5:15])/11), digits=2)                     ## para 15 a 65+
tab1[[i]]$cob3<-round((sum(tab1[[i]]$Cx[6:15])/10), digits=2)                     ## para 20 a 65+
tab1[[i]]$cob4<-round((sum(tab1[[i]]$Cx[7:15])/9), digits=2)                      ## para 25 a 65+
tab1[[i]]$cob5<-round((sum(tab1[[i]]$Cx[8:15])/8), digits=2)                      ## para 30 a 65+
                            }
tab2<-unsplit(tab1,cod)
grau1<-tapply(tab2$cob1,tab2$cod,mean)
grau2<-tapply(tab2$cob2,tab2$cod,mean)
grau3<-tapply(tab2$cob3,tab2$cod,mean)
grau4<-tapply(tab2$cob4,tab2$cod,mean)
grau5<-tapply(tab2$cob5,tab2$cod,mean)
resultado<-cbind(grau1,grau2,grau3,grau4,grau5)
results<-list(resultado=resultado)
return(results)}
