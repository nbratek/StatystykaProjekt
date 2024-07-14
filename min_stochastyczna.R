library(smoof)
library(tidyverse) ##do przeksztalcania danych
library(ecr) ## GA 
library(foreach) ## do zrowwnoleglania obliczen
library(doParallel) ## do zrowwnoleglania obliczen



# PRS ---------------------------------------------------------------------


PRS <- function(l_wymiarow=2,funkcja="Rosenbrock",liczba_wywolan=1000) {

  fn <- if (funkcja=="Rosenbrock") {
    makeRosenbrockFunction(l_wymiarow)
  } else {makeRastriginFunction(l_wymiarow)}
  
  lower = getLowerBoxConstraints(fn)
  upper = getUpperBoxConstraints(fn)
  
  x <- rerun(liczba_wywolan,fn(runif(l_wymiarow,lower,upper))) 

 return(min(lapply(x,'[',1) %>%  unlist))
  
}


# GA ----------------------------------------------------------------------


alg_gen <-  function(l_wymiarow=2,funkcja="Ackley",liczba_wywolan=1000) {
  
  fn <- if (funkcja=="Ackley") {
    makeAckleyFunction(l_wymiarow)
  } else {makeRastriginFunction(l_wymiarow)}
  
  MU = 30L
  LAMBDA = 5L
  MAX.ITER = liczba_wywolan
  lower = getLowerBoxConstraints(fn)
  upper = getUpperBoxConstraints(fn)
  
  
  ecr(fitness.fun = fn, representation = "float",
      n.dim = getNumberOfParameters(fn), survival.strategy = "plus", minimize = TRUE,
      lower = lower, upper = upper,
      mu = MU, lambda = LAMBDA,
      mutator = setup(mutGauss, sdev = 2, lower = lower, upper = upper),
      terminators = list(stopOnIters(MAX.ITER))) $best.y
  
}


alg_gen()


# generowanie wynikow --------------------------------------------------------------


run_50_PRS <- function(wymiar,nazwa_funkcji) {
  rerun(50,PRS(wymiar,nazwa_funkcji)) %>% 
  unlist() %>%
  as.data.frame() %>% 
  mutate(algorytm="poszukiwanie_przypadkowe") ##dodawanie kolumny
}

run_50_alg_gen <- function(wymiar,nazwa_funkcji) {
  rerun(50,alg_gen(wymiar,nazwa_funkcji)) %>% 
  unlist() %>%
  as.data.frame() %>% 
  mutate(algorytm="algorytm_genetyczny")
}



grid <- expand.grid(wymiar <- c(2,10,20)
                    ,nazwa_funkcji <- c("Rosenbrock","Rastrigin")
                    )

myCluster = makeCluster(4) #liczba rdzeni do wykorzystania
registerDoParallel(myCluster)
clusterCall(myCluster, function() library(tidyverse))
clusterCall(myCluster, function() library(smoof))
clusterCall(myCluster, function() library(ecr))

clusterExport(myCluster,varlist = c("alg_gen"
                                    ,"PRS"
                                    ,"run_50_PRS"
                                    ,"run_50_alg_gen"
                                    ))

wynik_alg_gen <- plyr::mdply(.data=unname(grid)
             ,.fun=run_50_alg_gen
             ,.parallel = TRUE
)

wynik_PRS <- plyr::mdply(.data=unname(grid)
                       ,.fun=run_50_PRS
                       ,.parallel = TRUE
)


stopCluster(myCluster)
gc()


wynik_razem <- 
  bind_rows(wynik_alg_gen,wynik_PRS) %>% 
 `colnames<-`(c("liczba_wymiarow", "funkcja", "wartosc","algorytm"))



# wykresy -----------------------------------------------------------------



wynik_razem %>% 
  ggplot() +
  geom_boxplot(aes(x=wartosc,group=algorytm,fill=algorytm)) +
  facet_wrap(funkcja~liczba_wymiarow,scales = "free") +
  theme_bw()


wynik_razem %>% 
  ggplot() +
  geom_histogram(aes(x=wartosc,group=algorytm,fill=algorytm)) +
  facet_wrap(funkcja~liczba_wymiarow,scales = "free") +
  theme_bw()


# istotnosc_statystyczna --------------------------------------------------

options(scipen=100000)

porownanie_srednich <- 
  wynik_razem %>% 
  group_by(funkcja,algorytm,liczba_wymiarow) %>% 
  summarise(sr_wartosc=mean(wartosc)) 

wyniki_pivot <- 
wynik_razem %>% 
  pivot_wider(names_from = algorytm,values_from=wartosc)

rm(ttest_wyniki)
ttest_wyniki <- list()

for (i in seq(1,nrow(wyniki_pivot))) {
  
     tst <- t.test(unlist(wyniki_pivot$algorytm_genetyczny[i]),
                 unlist(wyniki_pivot$poszukiwanie_przypadkowe[i])
                 ,paired = F
                 )


    ttest_wyniki <- rbind(tst,ttest_wyniki)
  
}

wyniki_istotnosc_stat <- 
  bind_cols(wyniki_pivot,ttest_wyniki %>% as.data.frame()) %>% 
  unnest(p.value)


# MS ----------------------------------------------------------------------

MS <-  function(l_wymiarow=2,funkcja="Schwefel",liczba_wywolan=100) {
  
  fn <- if (funkcja=="Schwefel") {
    makeSchwefelFunction(l_wymiarow)
  } else {makeRastriginFunction(l_wymiarow)}
  
  lower = getLowerBoxConstraints(fn)
  upper = getUpperBoxConstraints(fn)
  
  x <- rerun(liczba_wywolan
             ,optim(par = runif(l_wymiarow,-500,500)
                    ,fn = fn
                    ,method = "L-BFGS-B"
             )
  ) 
  
  vals <- lapply(x,'[[',2) %>% unlist
  counts <- lapply(lapply(x,'[[',3),'[[',1) %>% unlist
  
  c(min(vals),mean(vals),mean(counts))
  
}

MS()

