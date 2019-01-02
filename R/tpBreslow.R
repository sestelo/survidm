#library(survival)
tpBreslow<-

  function (object, s, t, z.value, z.name, parte2=NULL, conf = FALSE, conf.type='log', conf.level = 0.95, n.boot = 199)

  {

    if (missing(object))
      stop("Argument 'object' is missing, with no default")
    # if (!inherits(object, "survIDM")) stop("'object' must be of class 'survIDM'")
    if (missing(s)) s <- 0

    obj <- object[[1]]

    #dim(obj) #911

    ptimes <- which(obj$event1 == 1 | obj$event == 1)

    if (missing(t)) t <- obj$Stime[ptimes]

    #t <- obj$Stime[ptimes]

    if (any(t <= 0)) stop("The values of 't' must be positive")
    if (s > max(obj$time1)) stop("The value of 's' is too large")
    if (s < 0) stop("'s' must be nonnegative")
    if (length(s) > 1) stop("Length of 's' must be 1")

    t <- t[t >= s]
    t <- sort(unique(t))

    if(length(t) == 0) stop("Invalid values for 't'.")

    p00 <- rep(NA, length(t))
    p01 <- rep(NA, length(t))
    p02 <- rep(NA, length(t))
    p11 <- rep(NA, length(t))
    p12 <- rep(NA, length(t))


    if(length(z.name)==1 & !('pspline' %in% substr(parte2,1,7))){

      if(substr(colnames(obj)[5],1,6)=='factor'){

        colnames(obj)[5]<-substr(colnames(obj)[5],8,nchar(colnames(obj)[5])-1)
      }

      p <- which(obj$time1>s)

      obj2 <- obj[p,]

      p <- which(obj$time1 <= s & obj$Stime > s)

      obj3<- obj[p,]

      new<-data.frame(z.value)

      colnames(new)<-colnames(obj)[5]
      #Breslow Pij(s,t|Age)

      if(conf.type=='log' | conf.type=='plain' | conf.type=='log-log' | conf.type=='logit'){


        formula<-as.formula(paste("Surv(time1, event1) ~ ",colnames(obj)[5], sep=''))

        p00 <- coxph(formula, data =  obj2)

        new<-data.frame(z.value)

        colnames(new)<-colnames(obj)[5]


        p00.s <- survfit(p00,newdata=new, conf.type=conf.type, conf.int=conf.level) #para todos os tipos de conf.type

        p00.s.t_sum<-summary(p00.s, time=t)

        times1<-p00.s.t_sum$time

        p00.s.t<-p00.s.t_sum$surv

        p00.s.t.l<-p00.s.t_sum$lower

        p00.s.t.u<-p00.s.t_sum$upper


        #

        formula2<-as.formula(paste("Surv(Stime, event) ~ ",colnames(obj)[5], sep=''))


        p01_1 <- coxph(formula2, data =  obj2)

        p01_1.s <- survfit(p01_1,newdata=new, conf.type=conf.type, conf.int=conf.level)

        p01_1.s.t<-summary(p01_1.s , time=t)

        p01.s.t<-p01_1.s.t$surv-p00.s.t

        p02.s.t<-1-p01_1.s.t$surv

        #p02.s.t + p01.s.t + p00.s.t


        if(conf==TRUE){

          #bootstrap p01:

          res.ci <- array(NA, dim=c(length(t), n.boot, 1))

          p01.ci <- matrix(NA, length(t), 2)

          n <- dim(object[[1]])[1]

          for (j in 1.:n.boot){

            #print(paste('--------------------------------------------------count: ', j))

            #j<-1
            xx <- sample.int(n, size = n, replace = TRUE)
            ndata <- object[[1]][xx,]

            if(substr(colnames(ndata )[5],1,6)=='factor'){

              colnames(ndata )[5]<-substr(colnames(ndata )[5],8,nchar(colnames(ndata )[5])-1)
            }

            p0 <- which(ndata$time1 > s)
            p1 <- which(ndata$time1 <= s & ndata$Stime > s)
            obj2b <- ndata [p0,]
            obj3b<- ndata [p1,]

            formula<-as.formula(paste("Surv(time1, event1) ~ ",colnames(obj)[5], sep=''))

            p00b <- coxph(formula, data =  obj2b)

            p00b.s <- survfit(p00b,conf.type='log',newdata=new, conf.int=conf.level)

            p00b.s.t_sum<-summary(p00b.s, time=t)

            times1b<-p00b.s.t_sum$time

            p00b.s.t<-p00b.s.t_sum$surv

            formula2<-as.formula(paste("Surv(Stime, event) ~ ",colnames(obj)[5], sep=''))


            p01b_1 <- coxph(formula2, data =  obj2b)

            p01b_1.s <- survfit(p01b_1, newdata=new,conf.type='log', conf.int=conf.level)

            p01b_1.s.t<-summary(p01b_1.s , time=t)

            p01b.s.t<-p01b_1.s.t$surv-p00b.s.t

            p02b.s.t<-1-p01b_1.s.t$surv

            #p00b.s.t+p01b.s.t+ p02b.s.t

            tabF<-as.data.frame(cbind(times1b,p00b.s.t,p01b.s.t,p02b.s.t))
            names(tabF)<-c('times','p00b.s.t','p01b.s.t','p02b.s.t')


            for (k in 1: length(t)) {

              res.ci[k, j, 1]<-tabF[k,3]

            }

          }

          for (k in 1: length(t)) {
            p01.ci[k,1] <- quantile(res.ci[k,,1], (1 - conf.level) / 2,na.rm=T)
            p01.ci[k,2] <- quantile(res.ci[k,,1], 1 - (1 - conf.level) / 2,na.rm=T)

          }


          p01.s.t.l<- p01.ci[,1]
          p01.s.t.u<-p01.ci[,2]

        }else{
          p01.s.t.l<-rep(NA,length(t))
          p01.s.t.u<-rep(NA,length(t))

        } #end conf=TRUE

        #-------------------
        #p02:

        p02.s.t.l<-1-p01_1.s.t$surv+qnorm((1-conf.level)/2)*p01_1.s.t$std.err
        p02.s.t.u<-1-p01_1.s.t$surv+qnorm(conf.level+(1-conf.level)/2)*p01_1.s.t$std.err


        #--------

        p11 <- coxph(formula2, data =  obj3)

        p11.s <- survfit(p11,newdata=new, conf.type=conf.type, conf.int=conf.level)

        p11.s.t_sum<-summary(p11.s, time=t)

        times2<-p11.s.t_sum$time

        p11.s.t<-p11.s.t_sum$surv

        p12.s.t<-1-p11.s.t

        #p11.s.t+p12.s.t

        p11.s.t.l<-p11.s.t_sum$lower

        p11.s.t.u<-p11.s.t_sum$upper

        p12.s.t.l<-1-p11.s.t_sum$surv+qnorm((1-conf.level)/2)*p11.s.t_sum$std.err

        p12.s.t.u<-1-p11.s.t_sum$surv+qnorm(conf.level+(1-conf.level)/2)*p11.s.t_sum$std.err


      }else{#bootstrap

        formula<-as.formula(paste("Surv(time1, event1) ~ ",colnames(obj)[5], sep=''))

        p00 <- coxph(formula, data =  obj2)


        p00.s <- survfit(p00, newdata=new,conf.type='log',  conf.int=conf.level)

        p00.s.t_sum<-summary(p00.s, time=t)

        times1<-p00.s.t_sum$time

        p00.s.t<-p00.s.t_sum$surv

        formula2<-as.formula(paste("Surv(Stime, event) ~ ",colnames(obj)[5], sep=''))

        p01_1 <- coxph(formula2, data =  obj2)

        p01_1.s <- survfit(p01_1, newdata=new,conf.type='log',conf.int=conf.level)

        p01_1.s.t<-summary(p01_1.s , time=t)

        p01.s.t<-p01_1.s.t$surv-p00.s.t

        p02.s.t<-1-p01_1.s.t$surv

        #p02.s.t + p01.s.t + p00.s.t

        #length(p02.s.t) #392

        p11 <- coxph(formula2, data =  obj3)

        p11.s <- survfit(p11, newdata=new, conf.type='log',conf.int=conf.level)

        p11.s.t_sum<-summary(p11.s, time=t)

        times2<-p11.s.t_sum$time

        p11.s.t<-p11.s.t_sum$surv

        p12.s.t<-1-p11.s.t

        #p11.s.t+p12.s.t

        length(p00.s.t) #383

        length(p11.s.t) #378

        which(!times1 %in% times2)

        p00.s.t<-p00.s.t[which(times1 %in% times2)]

        p01.s.t<-p01.s.t[which(times1 %in% times2)]

        p02.s.t<-p02.s.t[which(times1 %in% times2)]


        #to estimate conf intervals using bootstrap
        res.ci <- array(NA, dim=c(length(times2), n.boot, 5))

        p00.ci <- matrix(NA, length(times2), 2)
        p01.ci <- matrix(NA, length(times2), 2)
        p02.ci <- matrix(NA, length(times2), 2)
        p11.ci <- matrix(NA, length(times2), 2)
        p12.ci <- matrix(NA, length(times2), 2)

        n <- dim(object[[1]])[1]

        for (j in 1.:n.boot){

          #j<-1
          #print(paste('--------------------------------------------------count: ', j))


          xx <- sample.int(n, size = n, replace = TRUE)
          ndata <- object[[1]][xx,]

          if(substr(colnames(ndata )[5],1,6)=='factor'){

            colnames(ndata )[5]<-substr(colnames(ndata )[5],8,nchar(colnames(ndata )[5])-1)
          }
          p0 <- which(ndata$time1 > s)


          p1 <- which(ndata$time1 <= s & ndata$Stime > s)


          obj2b <- ndata [p0,]
          obj3b<- ndata [p1,]


          formula<-as.formula(paste("Surv(time1, event1) ~ ",colnames(obj)[5], sep=''))

          p00b <- coxph(formula, data =  obj2b)

          p00b.s <- survfit(p00b, newdata=new, conf.type='log', conf.int=conf.level)

          p00b.s.t_sum<-summary(p00b.s, time=t)

          times1b<-p00b.s.t_sum$time


          p00b.s.t<-p00b.s.t_sum$surv

          formula2<-as.formula(paste("Surv(Stime, event) ~ ",colnames(obj)[5], sep=''))

          p01b_1 <- coxph(formula2, data =  obj2b)

          p01b_1.s <- survfit(p01b_1, newdata=new,conf.type='log', conf.int=conf.level)

          p01b_1.s.t<-summary(p01b_1.s , time=t)

          p01b.s.t<-p01b_1.s.t$surv-p00b.s.t

          p02b.s.t<-1-p01b_1.s.t$surv

          #

          p11b <- coxph(formula2, data =  obj3b)

          p11b.s <- survfit(p11b, newdata=new,conf.type='log', conf.int=conf.level)

          p11b.s.t_sum<-summary(p11b.s, time=t)

          times2b<-p11b.s.t_sum$time

          p11b.s.t<-p11b.s.t_sum$surv

          p12b.s.t<-1-p11b.s.t

          tab1<-as.data.frame(cbind(times1b,p00b.s.t,p01b.s.t,p02b.s.t))
          names(tab1)<-c('times','p00b.s.t','p01b.s.t','p02b.s.t')

          #dim(tab1[tab1$times %in% times2b,])

          tab1<-tab1[tab1$times %in% times2b,]


          tab1<-tab1[tab1$times %in% times2b,]

          tab2<-as.data.frame(cbind(times2b,p11b.s.t,p12b.s.t))
          names(tab2)<-c('times','p11b.s.t','p12b.s.t')

          tabF<-merge(tab1,tab2,by='times',all = T)

          for (k in 1: length(times2)) {

            res.ci[k, j, 1]<-tabF[k,2]
            res.ci[k, j, 2]<-tabF[k,3]
            res.ci[k, j, 3]<-tabF[k,4]
            res.ci[k, j, 4]<-tabF[k,5]
            res.ci[k, j, 5]<-tabF[k,6]

          }
        }#end bootstrap


        for (k in 1: length(times2)) {

          p00.ci[k,1] <- quantile(res.ci[k,,1], (1 - conf.level) / 2,na.rm=T)
          p00.ci[k,2] <- quantile(res.ci[k,,1], 1 - (1 - conf.level) / 2,na.rm=T)

          p01.ci[k,1] <- quantile(res.ci[k,,2], (1 - conf.level) / 2,na.rm=T)
          p01.ci[k,2] <- quantile(res.ci[k,,2], 1 - (1 - conf.level) / 2,na.rm=T)

          p02.ci[k,1] <- quantile(res.ci[k,,3], (1 - conf.level) / 2,na.rm=T)
          p02.ci[k,2] <- quantile(res.ci[k,,3], 1 - (1 - conf.level) / 2,na.rm=T)


          p11.ci[k,1] <- quantile(res.ci[k,,4], (1 - conf.level) / 2,na.rm=T)
          p11.ci[k,2] <- quantile(res.ci[k,,4], 1 - (1 - conf.level) / 2,na.rm=T)

          p12.ci[k,1] <- quantile(res.ci[k,,5], (1 - conf.level) / 2,na.rm=T)
          p12.ci[k,2] <- quantile(res.ci[k,,5], 1 - (1 - conf.level) / 2,na.rm=T)

        }


        p00.s.t.l<- p00.ci[,1]
        p00.s.t.u<-p00.ci[,2]

        p01.s.t.l<- p01.ci[,1]
        p01.s.t.u<-p01.ci[,2]

        p02.s.t.l<- p02.ci[,1]
        p02.s.t.u<-p02.ci[,2]

        p11.s.t.l<- p11.ci[,1]
        p11.s.t.u<-p11.ci[,2]

        p12.s.t.l<- p11.ci[,1]
        p12.s.t.u<-p12.ci[,2]


      }#end bootstrap (else)

      index<-which(times1 %in% times2)


      if(conf.type!='bootstrap'){
        #table(times2 %in% times1) #all true

        #table(times1 %in% times2) #5 false

        #index<-which(times2 %in% times1) #previous version

        index<-which(times1 %in% times2)
        #which(!times1 %in% times2) #the last one

        p00 <- p00.s.t[index]
        p01 <- p01.s.t[index]
        p02 <- p02.s.t[index]
        p11 <- p11.s.t
        p12 <- p12.s.t
        #p11 <- p11.s.t[index]
        #p12 <- p12.s.t[index]

        #IC (lower)

        p00.l <-p00.s.t.l[index]
        p01.l <-p01.s.t.l[index]
        p02.l <-p02.s.t.l[index]
        p11.l <-p11.s.t.l
        p12.l <-p12.s.t.l
        #p11.l <-p11.s.t.l[index]
        #p12.l <-p12.s.t.l[index]
        #IC (upper)

        p00.u <-p00.s.t.u[index]
        p01.u <-p01.s.t.u[index]
        p02.u <-p02.s.t.u[index]
        #p11.u <-p11.s.t.u[index]
        #p12.u <-p12.s.t.u[index]
        p11.u <-p11.s.t.u
        p12.u <-p12.s.t.u


      }else{

        p00 <- p00.s.t
        p01 <- p01.s.t
        p02 <- p02.s.t
        p11 <- p11.s.t
        p12 <- p12.s.t

        #IC (lower)

        p00.l <-p00.s.t.l
        p01.l <-p01.s.t.l
        p02.l <-p02.s.t.l
        p11.l <-p11.s.t.l
        p12.l <-p12.s.t.l

        #IC (upper)

        p00.u <-p00.s.t.u
        p01.u <-p01.s.t.u
        p02.u <-p02.s.t.u
        p11.u <-p11.s.t.u
        p12.u <-p12.s.t.u

      }

      #------------------------------
      #ARRAY

      #est

      resS.all.probs_f<-array(NA, c(length(index),3,5))

      resS.all.probs_f[,1,1]<-p00

      resS.all.probs_f[,1,2]<-p01

      resS.all.probs_f[,1,3]<-p02

      #resS.all.probs_f[1:5,1,]

      resS.all.probs_f[,1,4]<-p11

      resS.all.probs_f[,1,5]<-p12

      #lower:

      resS.all.probs_f[,2,1]<-p00.l
      resS.all.probs_f[,2,2]<-p01.l
      resS.all.probs_f[,2,3]<-p02.l
      resS.all.probs_f[,2,4]<-p11.l
      resS.all.probs_f[,2,5]<-p12.l

      #upper:

      resS.all.probs_f[,3,1]<-p00.u
      resS.all.probs_f[,3,2]<-p01.u
      resS.all.probs_f[,3,3]<-p02.u
      resS.all.probs_f[,3,4]<-p11.u
      resS.all.probs_f[,3,5]<-p12.u


      if(conf==FALSE){

        resu <- data.frame(cbind(t[index], p00, p01, p02, p11, p12))

        names(resu) <- c("t", "p00", "p01", "p02", "p11", "p12")

        result <- list(est = resu,s = s, t = t, conf = conf,z.name = z.name,
                       z.value=z.value,conf.level = conf.level)


      }else{

        #est

        resu <- data.frame(cbind(t[index], p00, p01, p02, p11, p12))

        names(resu) <- c("t", "p00", "p01", "p02", "p11", "p12")

        #ci
        auxci <- resS.all.probs_f[,2:3,]

        suppressWarnings(auxci <- data.frame(matrix(auxci, ncol = 10, nrow =length(resu$t))))

        names(auxci) <- c("p00.li.ci", "p00.ls.ci", "p01.li.ci", "p01.ls.ci",
                          "p02.li.ci", "p02.ls.ci", "p11.li.ci", "p11.ls.ci",
                          "p12.li.ci", "p12.ls.ci")

        result <- list(est = resu,CI = auxci,s = s, t = t, conf = conf,z.name = z.name,
                       z.value=z.value,conf.level = conf.level)  #atual

      }


      class(result) <- c("Breslow", "tp")


      return(invisible(result))


    }else{

      p <- which(obj$time1>s)

      obj2 <- obj[p,]

      p <- which(obj$time1 <= s & obj$Stime > s)

      obj3<- obj[p,]


      for1<-NULL

      for (i in 1:length(parte2)){

        for1<-paste(for1,parte2[i],sep='+')
      }

      for1<-substr(for1,2,nchar(for1))


      if(!is.null(z.value)){

        ind<-as.numeric(which(sapply(obj[(ncol(obj)-length(z.value)+1):ncol(obj)], class)=='numeric'))

        new<-as.data.frame(as.matrix(t(z.value)))

        colnames(new)<-z.name

        for (j in 1:ncol(new)){


          if(j %in% ind){

            new[,j]<-as.numeric(as.character(new[,j]))

          }else{

            new[,j]<-as.character(new[,j])

          }

        }
      }

      if(conf.type=='log' | conf.type=='plain' | conf.type=='log-log' | conf.type=='logit'){

        formula<-as.formula(paste("Surv(time1, event1) ~ ", for1, sep=''))

        p00 <- coxph(formula, data =  obj2)

        if(!is.null(z.value)){

          p00.s <- survfit(p00, newdata=new, conf.type=conf.type, conf.int=conf.level)

        }else{

          p00.s <- survfit(p00,conf.type=conf.type, conf.int=conf.level)

        }

        p00.s.t_sum<-summary(p00.s, time=t)

        #summary(p00.s, time=365*2:5)

        times1<-p00.s.t_sum$time

        p00.s.t<-p00.s.t_sum$surv

        p00.s.t.l<-p00.s.t_sum$lower

        p00.s.t.u<-p00.s.t_sum$upper

        #

        formula2<-as.formula(paste("Surv(Stime, event) ~ ",for1, sep=''))


        p01_1 <- coxph(formula2, data =  obj2)

        if(!is.null(z.value)){

          p01_1.s <- survfit(p01_1, newdata=new,conf.type=conf.type, conf.int=conf.level)

        }else{

          p01_1.s <- survfit(p01_1,conf.type=conf.type, conf.int=conf.level)

        }

        p01_1.s.t<-summary(p01_1.s , time=t)

        p01.s.t<-p01_1.s.t$surv-p00.s.t

        p02.s.t<-1-p01_1.s.t$surv

        #p02.s.t + p01.s.t + p00.s.t


        if(conf==TRUE){

          #bootstrap  p01:
          res.ci <- array(NA, dim=c(length(t), n.boot, 1))

          p01.ci <- matrix(NA, length(t), 2)

          n <- dim(object[[1]])[1]

          for (j in 1.:n.boot){

            #print(paste('--------------------------------------------------count: ', j))

            xx <- sample.int(n, size = n, replace = TRUE)
            ndata <- object[[1]][xx,]

            p0 <- which(ndata$time1 > s)
            p1 <- which(ndata$time1 <= s & ndata$Stime > s)
            obj2b <- ndata [p0,]
            obj3b<- ndata [p1,]

            for1<-NULL

            for (i in 1:length(parte2)){

              for1<-paste(for1,parte2[i],sep='+')
            }

            for1<-NULL

            for (i in 1:length(parte2)){

              for1<-paste(for1,parte2[i],sep='+')
            }

            for1<-substr(for1,2,nchar(for1))


            formula<-as.formula(paste("Surv(time1, event1) ~ ", for1, sep=''))


            p00b <- coxph(formula, data =  obj2b)


            if(!is.null(z.value)){

              p00b.s <- survfit(p00b, newdata=new, conf.type='log',conf.int=conf.level)

            }else{

              p00b.s <- survfit(p00b,conf.type='log',conf.int=conf.level)

            }

            p00b.s.t_sum<-summary(p00b.s, time=t)

            times1b<-p00b.s.t_sum$time


            p00b.s.t<-p00b.s.t_sum$surv


            formula2<-as.formula(paste("Surv(Stime, event) ~ ",for1, sep=''))


            p01b_1 <- coxph(formula2, data =  obj2b)

            if(!is.null(z.value)){

              p01b_1.s <- survfit(p01b_1, newdata=new,conf.type='log',conf.int=conf.level)

            }else{

              p01b_1.s <- survfit(p01b_1,conf.type='log',conf.int=conf.level)

            }

            p01b_1.s.t<-summary(p01b_1.s , time=t)

            p01b.s.t<-p01b_1.s.t$surv-p00b.s.t

            p02b.s.t<-1-p01b_1.s.t$surv

            p00b.s.t+p01b.s.t+ p02b.s.t

            tabF<-as.data.frame(cbind(times1b,p00b.s.t,p01b.s.t,p02b.s.t))
            names(tabF)<-c('times','p00b.s.t','p01b.s.t','p02b.s.t')

            for (k in 1: length(t)) {

              res.ci[k, j, 1]<-tabF[k,3]

            }

          }#end bootstrap


          for (k in 1: length(t)) {
            p01.ci[k,1] <- quantile(res.ci[k,,1], (1 - conf.level) / 2,na.rm=T)
            p01.ci[k,2] <- quantile(res.ci[k,,1], 1 - (1 - conf.level) / 2,na.rm=T)

          }


          p01.s.t.l<- p01.ci[,1]
          p01.s.t.u<-p01.ci[,2]

        }else{

          p01.s.t.l<-rep(NA,length(t))
          p01.s.t.u<-rep(NA,length(t))

        }

        #-------------------
        #p02:

        p02.s.t.l<-1-p01_1.s.t$surv+qnorm((1-conf.level)/2)*p01_1.s.t$std.err
        p02.s.t.u<-1-p01_1.s.t$surv+qnorm(conf.level+(1-conf.level)/2)*p01_1.s.t$std.err


        #--------

        p11 <- coxph(formula2, data =  obj3)

        if(!is.null(z.value)){
          p11.s <- survfit(p11, newdata=new, conf.type=conf.type,conf.int=conf.level)

        }else{
          p11.s<- survfit(p11,conf.int=conf.level)

        }

        p11.s.t_sum<-summary(p11.s, time=t)

        #summary(p11.s, time=365*2:5)

        times2<-p11.s.t_sum$time

        p11.s.t<-p11.s.t_sum$surv

        p12.s.t<-1-p11.s.t

        #

        p11.s.t.l<-p11.s.t_sum$lower

        p11.s.t.u<-p11.s.t_sum$upper


        p12.s.t.l<-1-p11.s.t_sum$surv+qnorm((1-conf.level)/2)*p11.s.t_sum$std.err

        p12.s.t.u<-1-p11.s.t_sum$surv+qnorm(conf.level+(1-conf.level)/2)*p11.s.t_sum$std.err

      }else{#bootstrap


        formula<-as.formula(paste("Surv(time1, event1) ~ ", for1, sep=''))

        p00 <- coxph(formula, data =  obj2)

        if(!is.null(z.value)){

          p00.s <- survfit(p00, newdata=new, conf.type='log',conf.int=conf.level)

        }else{

          p00.s <- survfit(p00,conf.type='log',conf.int=conf.level)

        }

        p00.s.t_sum<-summary(p00.s, time=t)

        times1<-p00.s.t_sum$time

        p00.s.t<-p00.s.t_sum$surv

        formula2<-as.formula(paste("Surv(Stime, event) ~ ",for1, sep=''))

        p01_1 <- coxph(formula2, data =  obj2)

        if(!is.null(z.value)){

          p01_1.s <- survfit(p01_1, newdata=new,conf.type='log',conf.int=conf.level)

        }else{

          p01_1.s <- survfit(p01_1,conf.type='log',conf.int=conf.level)

        }

        p01_1.s.t<-summary(p01_1.s , time=t)

        p01.s.t<-p01_1.s.t$surv-p00.s.t

        p02.s.t<-1-p01_1.s.t$surv

        #p02.s.t + p01.s.t + p00.s.t

        p11 <- coxph(formula2, data =  obj3)

        if(!is.null(z.value)){

          p11.s <- survfit(p11, newdata=new, conf.type='log',conf.int=conf.level)

        }else{

          p11.s <- survfit(p11,conf.type='log',conf.int=conf.level)

        }

        p11.s.t_sum<-summary(p11.s, time=t)

        times2<-p11.s.t_sum$time

        p11.s.t<-p11.s.t_sum$surv

        p12.s.t<-1-p11.s.t

        #p11.s.t+p12.s.t

        length(p00.s.t) #383

        length(p11.s.t) #378

        which(!times1 %in% times2)

        p00.s.t<-p00.s.t[which(times1 %in% times2)]

        p01.s.t<-p01.s.t[which(times1 %in% times2)]

        p02.s.t<-p02.s.t[which(times1 %in% times2)]

        #conf intervals using bootstrap

        res.ci <- array(NA, dim=c(length(times2), n.boot, 5))

        p00.ci <- matrix(NA, length(times2), 2)
        p01.ci <- matrix(NA, length(times2), 2)
        p02.ci <- matrix(NA, length(times2), 2)
        p11.ci <- matrix(NA, length(times2), 2)
        p12.ci <- matrix(NA, length(times2), 2)

        n <- dim(object[[1]])[1]

        for (j in 1.:n.boot){
          #j<-1
          #print(paste('--------------------------------------------------count: ', j))

          xx <- sample.int(n, size = n, replace = TRUE)
          ndata <- object[[1]][xx,]

          p0 <- which(ndata$time1 > s)
          p1 <- which(ndata$time1 <= s & ndata$Stime > s)

          obj2b <- ndata [p0,]
          obj3b<- ndata [p1,]


          for1<-NULL

          for (i in 1:length(parte2)){

            for1<-paste(for1,parte2[i],sep='+')
          }

          for1<-substr(for1,2,nchar(for1))


          formula<-as.formula(paste("Surv(time1, event1) ~ ", for1, sep=''))


          p00b <- coxph(formula, data =  obj2b)


          if(!is.null(z.value)){

            p00b.s <- survfit(p00b, newdata=new, conf.type='log',conf.int=conf.level)

          }else{

            p00b.s <- survfit(p00b,conf.type='log',conf.int=conf.level)

          }

          p00b.s.t_sum<-summary(p00b.s, time=t)

          times1b<-p00b.s.t_sum$time

          p00b.s.t<-p00b.s.t_sum$surv


          formula2<-as.formula(paste("Surv(Stime, event) ~ ",for1, sep=''))


          p01b_1 <- coxph(formula2, data =  obj2b)

          if(!is.null(z.value)){

            p01b_1.s <- survfit(p01b_1, newdata=new,conf.type='log',conf.int=conf.level)

          }else{

            p01b_1.s <- survfit(p01b_1,conf.type='log',conf.int=conf.level)

          }

          p01b_1.s.t<-summary(p01b_1.s , time=t)

          p01b.s.t<-p01b_1.s.t$surv-p00b.s.t

          p02b.s.t<-1-p01b_1.s.t$surv

          #

          p11b <- coxph(formula2, data =  obj3b)

          if(!is.null(z.value)){
            p11b.s <- survfit(p11b, newdata=new,conf.int=conf.level)

          }else{
            p11b.s<- survfit(p11b,conf.int=conf.level)

          }

          p11b.s.t_sum<-summary(p11b.s, time=t)

          times2b<-p11b.s.t_sum$time

          p11b.s.t<-p11b.s.t_sum$surv

          p12b.s.t<-1-p11b.s.t

          tab1<-as.data.frame(cbind(times1b,p00b.s.t,p01b.s.t,p02b.s.t))
          names(tab1)<-c('times','p00b.s.t','p01b.s.t','p02b.s.t')

          #dim(tab1[tab1$times %in% times2b,])

          tab1<-tab1[tab1$times %in% times2b,]

          tab2<-as.data.frame(cbind(times2b,p11b.s.t,p12b.s.t))
          names(tab2)<-c('times','p11b.s.t','p12b.s.t')

          tabF<-merge(tab1,tab2,by='times',all = T)

          for (k in 1: length(times2)) {

            res.ci[k, j, 1]<-tabF[k,2]
            res.ci[k, j, 2]<-tabF[k,3]
            res.ci[k, j, 3]<-tabF[k,4]
            res.ci[k, j, 4]<-tabF[k,5]
            res.ci[k, j, 5]<-tabF[k,6]

          }

        }#end bootstrap

        for (k in 1: length(times2)) {

          p00.ci[k,1] <- quantile(res.ci[k,,1], (1 - conf.level) / 2,na.rm=T)
          p00.ci[k,2] <- quantile(res.ci[k,,1], 1 - (1 - conf.level) / 2,na.rm=T)

          p01.ci[k,1] <- quantile(res.ci[k,,2], (1 - conf.level) / 2,na.rm=T)
          p01.ci[k,2] <- quantile(res.ci[k,,2], 1 - (1 - conf.level) / 2,na.rm=T)

          p02.ci[k,1] <- quantile(res.ci[k,,3], (1 - conf.level) / 2,na.rm=T)
          p02.ci[k,2] <- quantile(res.ci[k,,3], 1 - (1 - conf.level) / 2,na.rm=T)


          p11.ci[k,1] <- quantile(res.ci[k,,4], (1 - conf.level) / 2,na.rm=T)
          p11.ci[k,2] <- quantile(res.ci[k,,4], 1 - (1 - conf.level) / 2,na.rm=T)

          p12.ci[k,1] <- quantile(res.ci[k,,5], (1 - conf.level) / 2,na.rm=T)
          p12.ci[k,2] <- quantile(res.ci[k,,5], 1 - (1 - conf.level) / 2,na.rm=T)

        }


        p00.s.t.l<- p00.ci[,1]
        p00.s.t.u<-p00.ci[,2]

        p01.s.t.l<- p01.ci[,1]
        p01.s.t.u<-p01.ci[,2]

        p02.s.t.l<- p02.ci[,1]
        p02.s.t.u<-p02.ci[,2]

        p11.s.t.l<- p11.ci[,1]
        p11.s.t.u<-p11.ci[,2]

        p12.s.t.l<- p11.ci[,1]
        p12.s.t.u<-p12.ci[,2]

      }#end bootstrap (else)

      index<-which(times1 %in% times2)


      if(conf.type!='bootstrap'){
        #table(times2 %in% times1) #all true

        #table(times1 %in% times2) #5 false

        #index<-which(times2 %in% times1) #previous version

        index<-which(times1 %in% times2)
        #which(!times1 %in% times2) #the last one

        p00 <- p00.s.t[index]
        p01 <- p01.s.t[index]
        p02 <- p02.s.t[index]
        p11 <- p11.s.t
        p12 <- p12.s.t
        #p11 <- p11.s.t[index]
        #p12 <- p12.s.t[index]

        #IC (lower)

        p00.l <-p00.s.t.l[index]
        p01.l <-p01.s.t.l[index]
        p02.l <-p02.s.t.l[index]
        p11.l <-p11.s.t.l
        p12.l <-p12.s.t.l
        #p11.l <-p11.s.t.l[index]
        #p12.l <-p12.s.t.l[index]
        #IC (upper)

        p00.u <-p00.s.t.u[index]
        p01.u <-p01.s.t.u[index]
        p02.u <-p02.s.t.u[index]
        #p11.u <-p11.s.t.u[index]
        #p12.u <-p12.s.t.u[index]
        p11.u <-p11.s.t.u
        p12.u <-p12.s.t.u


      }else{

        p00 <- p00.s.t
        p01 <- p01.s.t
        p02 <- p02.s.t
        p11 <- p11.s.t
        p12 <- p12.s.t

        #IC (lower)

        p00.l <-p00.s.t.l
        p01.l <-p01.s.t.l
        p02.l <-p02.s.t.l
        p11.l <-p11.s.t.l
        p12.l <-p12.s.t.l

        #IC (upper)

        p00.u <-p00.s.t.u
        p01.u <-p01.s.t.u
        p02.u <-p02.s.t.u
        p11.u <-p11.s.t.u
        p12.u <-p12.s.t.u

      }

      #------------------------------
      #ARRAY

      #est

      resS.all.probs_f<-array(NA, c(length(index),3,5))

      resS.all.probs_f[,1,1]<-p00

      resS.all.probs_f[,1,2]<-p01

      resS.all.probs_f[,1,3]<-p02

      resS.all.probs_f[,1,4]<-p11

      resS.all.probs_f[,1,5]<-p12


      #lower:

      resS.all.probs_f[,2,1]<-p00.l
      resS.all.probs_f[,2,2]<-p01.l
      resS.all.probs_f[,2,3]<-p02.l
      resS.all.probs_f[,2,4]<-p11.l
      resS.all.probs_f[,2,5]<-p12.l

      #upper:

      resS.all.probs_f[,3,1]<-p00.u
      resS.all.probs_f[,3,2]<-p01.u
      resS.all.probs_f[,3,3]<-p02.u
      resS.all.probs_f[,3,4]<-p11.u
      resS.all.probs_f[,3,5]<-p12.u


      if(conf==FALSE){

        #p00+p01+p02

        #p11+p12

        resu <- data.frame(cbind(t[index], p00, p01, p02, p11, p12))

        names(resu) <- c("t", "p00", "p01", "p02", "p11", "p12")

        head(resu)

        result <- list(est = resu,s = s, t = t, conf = conf,z.name = z.name,
                       z.value=z.value,conf.level = conf.level)



      }else{

        #est

        resu <- data.frame(cbind(t[index], p00, p01, p02, p11, p12))

        names(resu) <- c("t", "p00", "p01", "p02", "p11", "p12")

        #head(resu)

        #ci
        auxci <- resS.all.probs_f[,2:3,]

        suppressWarnings(auxci <- data.frame(matrix(auxci, ncol = 10, nrow =length(resu$t))))

        names(auxci) <- c("p00.li.ci", "p00.ls.ci", "p01.li.ci", "p01.ls.ci",
                          "p02.li.ci", "p02.ls.ci", "p11.li.ci", "p11.ls.ci",
                          "p12.li.ci", "p12.ls.ci")


        result <- list(est = resu,CI = auxci,s = s, t = t, conf = conf,z.name = z.name,
                       z.value=z.value,conf.level = conf.level)


      }


      class(result) <- c("Breslow", "tp")

      return(invisible(result))

    }

  }
