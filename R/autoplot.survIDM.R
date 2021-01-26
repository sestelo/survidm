ggplot2::autoplot

autoplot.survIDM <- function(x = object, y = NULL, trans = "all", func = "distribution",
                             conf = NULL, type = NULL,conftype = NULL, col = 1:6,
                             confcol = 1:6, lty = 1, conflty = 2, xlab = "Time (years)",
                             ylab = NULL, ylim = NULL, xlim = NULL, interactive = FALSE,...) {


  object <- x


  if(inherits(x, "survIDM") & class(x)[1] =='markov'){



    est <- object$TPestimates
    if(missing(xlab))  xlab <- paste(object$nm.method, "estimator")
    if(missing(ylab))  ylab <- paste("AJ estimator")



    #0->1


    p1<-ggplot(est, aes(est$aj01,est$nm01))+theme_bw()+labs(x = xlab) +labs(y = ylab)+labs(title='p01')+
      geom_point(aes(est$aj01, est$nm01))+geom_abline(intercept = 0, slope = 1,color='red',size=0.8)

    #0->2


    p2<-ggplot(est, aes(est$aj02,est$nm02))+theme_bw()+labs(x = xlab) +labs(y = ylab)+labs(title='p02')+
      geom_point(aes(est$aj02, est$nm02))+geom_abline(intercept = 0, slope = 1,color='red',size=0.8)

    #1->2


    p3<-ggplot(est, aes(est$aj12,est$nm12))+theme_bw()+labs(x = xlab) +labs(y = ylab)+labs(title='p12')+
      geom_point(aes(est$aj12, est$nm12))+geom_abline(intercept = 0, slope = 1,color='red',size=0.8)

    #1->2 NM with CI and AJ
    y <- cbind(est$nm12, est$nm12LCI, est$nm12UCI, est$aj12)
    x <- est$times


    p4<-ggplot(est, aes(est$times,est$nm12))+theme_bw()+labs(x = xlab) +labs(y = ylab)+
      geom_ribbon(aes(ymin=est$nm12LCI,ymax=est$nm12UCI),fill='gray60',alpha=0.7)+
      geom_line(aes(est$times,est$nm12),color='black',size=1)+
      geom_line(aes(est$times,est$aj12),color='red',size=1)

    grid.arrange(p1, p2, p3,p4,  layout_matrix = rbind(c(1, 2, 3), c(4, 4, 4)))

  }else{

    if (inherits(x, "survIDM")) {


      if (class(x)[1] == "data.frame") {

        plot(x)
      }


      if (!class(object)[1] %in% c("AJ", "LIDA", "LM", "PLM", "tpIPCW", "CIF",
                                   "cifIPCW", "soj", "sojIPCW", "LMAJ",
                                   "PLMAJ", "PAJ",'tpBreslow','markov')) {
        stop("The argumment 'Object' must be of one of the following classes
              'AJ', 'LIDA', 'LM', 'PLM', 'LMAJ', 'PLMAJ', 'PAJ', 'tpIPCW',
              'CIF', 'cifIPCW', 'soj', 'sojIPCW','tpBreslow' or 'markov'")
      }


      # for all
      #-----------------------------------------------

      object <- x

      object$Nlevels
      if (object$Nlevels != length(col))
        col <- rep(col, times = object$Nlevels)
      if (object$Nlevels != length(confcol))
        confcol <- rep(confcol, times = object$Nlevels)


      if (is.null(type))
        type <- "s"
      if (is.null(conftype))
        conftype <- "s"


      if (is.null(conf)) {
        ci <- object$conf
      } else {
        if (conf == TRUE & object$conf == FALSE) {
          stop("The surv object does not contain confidence intervals")
        }
        if (conf == TRUE & object$conf == TRUE)
          ci <- TRUE
        if (conf == FALSE)
          ci <- FALSE
      }



      if (is.null(ylim)) ylim <- c(0, 1)

      ob <- object$est

      obCI <- object$CI

      if (is.null(xlim) & object$Nlevels > 1) {
        xlim <- c(min(sapply(ob, function(x) min(x[, 1]), simplify = TRUE)),
                  max(sapply(ob, function(x) max(x[, 1]), simplify = TRUE)))
      }

      #---------------------------


      #tp
      #----------------------

      if (class(object)[1] %in% c("AJ", "LIDA", "LM", "PLM", "tpIPCW",
                                  "LMAJ", "PLMAJ", "PAJ",'tpBreslow')) {

        if(is.null(ylab) & class(object)[1] != "tpIPCW")
          ylab <- paste("p[ij](", (x$s), ",t)") #ylab <- bquote(paste('p[ij]', "(", .(x$s), ",t)"))

        if(is.null(ylab) & class(object)[1] == "tpIPCW")
          ylab <- paste("p[ij]", (x$s), ",t|", (x$z.name),")")#ylab <- bquote(paste(p[ij], "(", .(x$s), ",t|", .(x$z.name),")"))

        #-----------------------
        trans2 = trans
        tp <- c("00", "01", "02", "11", "12")
        if(trans == "all") {trans2 = tp}
        ii <- trans2 == tp
        itp <- 2:6
        itpCI <- c(1, 3, 5, 7, 9)
        #-----------------------

        if (object$Nlevels == 1) {

          if (ci == TRUE) {


            ob2<-as.data.frame(rbind(cbind(ob[,1],ob[,2],rep('00', length(ob[,1]))),
                                     cbind(ob[,1],ob[,3],rep('01', length(ob[,3]))),
                                     cbind(ob[,1],ob[,4],rep('02', length(ob[,4]))),
                                     cbind(ob[,1],ob[,5],rep('11', length(ob[,5]))),
                                     cbind(ob[,1],ob[,6],rep('12', length(ob[,6])))))



            names(ob2)<-c('time','TP','type')


            ob3<-obCI[itpCI[ii]]



            ob3.2<-as.data.frame(rbind(cbind(ob3[,1],rep('00', length(ob3[,1]))),
                                       cbind(ob3[,2],rep('01', length(ob3[,2]))),
                                       cbind(ob3[,3],rep('02', length(ob3[,3]))),
                                       cbind(ob3[,4],rep('11', length(ob3[,4]))),
                                       cbind(ob3[,5],rep('12', length(ob3[,5])))))


            names(ob3.2)<-c('tp_min','type')


            ob4<-obCI[itpCI[ii]+1]

            ob4.2<-as.data.frame(rbind(cbind(ob4[,1],rep('00', length(ob4[,1]))),
                                       cbind(ob4[,2],rep('01', length(ob4[,2]))),
                                       cbind(ob4[,3],rep('02', length(ob4[,3]))),
                                       cbind(ob4[,4],rep('11', length(ob4[,4]))),
                                       cbind(ob4[,5],rep('12', length(ob4[,5])))))

            names(ob4.2)<-c('tp_max','type')


            p1<-ggplot(ob2, aes(x=as.numeric(ob2$time), y=as.numeric(ob2$TP), group=factor(ob2$type),
                                fill=factor(ob2$type)))+theme_bw()+labs(x = xlab)+labs(y = ylab)


            p2<-p1+geom_ribbon(aes(ymin=as.numeric(ob3.2$tp_min),ymax=as.numeric(ob4.2$tp_max)),alpha=.4)


            p3<-p2+geom_line(aes(x=as.numeric(ob2$time), y=as.numeric(ob2$TP), color=ob2$type), size=1.1)+
              theme(legend.title=element_blank())


            if(isTRUE(interactive)){
              if (requireNamespace("plotly", quietly=TRUE)) {return(plotly::ggplotly(p3))}
            }else{
              return(p3)
            }



          }else{


            ob2<-as.data.frame(rbind(
              cbind(ob[,1],ob[,2],rep('00', length(ob[,2]))),
              cbind(ob[,1],ob[,3],rep('01', length(ob[,3]))),
              cbind(ob[,1],ob[,4],rep('02', length(ob[,4]))),
              cbind(ob[,1],ob[,5],rep('11', length(ob[,5]))),
              cbind(ob[,1],ob[,6],rep('12', length(ob[,6])))))


            names(ob2)<-c('time','TP','type')


            p3<-ggplot(ob2, aes(x=as.numeric(ob2$time), y=as.numeric(ob2$TP),group=factor(ob2$type),
                                color=factor(ob2$type)))



            p4<-p3+theme_bw()+labs(x = xlab)+labs(y = ylab)+geom_line(size=1)+ theme(legend.title=element_blank())

            if(isTRUE(interactive)){
              if (requireNamespace("plotly", quietly=TRUE)) {return(plotly::ggplotly(p4))}
            }else{
              return(p4)
            }


          }

          #if(trans == "all") legend("topright", c("00", "01", "02", "11", "12")[itp - 1], col = col, lty = lty)

        } else { # more than 1 level

          if (trans == "all") {
            stop(paste("The argumment 'trans' can't be 'all' if the factor", attr(terms(object$formula),"term.labels"),
                       "is included in the formula, you must select one of the transition probabilities."))
          }



          ob2<-NULL

          ob3.2<-NULL

          ob4.2<-NULL

          for(i in 1:object$Nlevels){

            #i<-1

            ob2<-as.data.frame(rbind(ob2,rbind(
              cbind(ob[[i]][,1],ob[[i]][,2],rep('00', length(ob[[i]][,2])),rep(object$levels[i], length(ob[[i]][,2]))),
              cbind(ob[[i]][,1],ob[[i]][,3],rep('01', length(ob[[i]][,3])),rep(object$levels[i], length(ob[[i]][,3]))),
              cbind(ob[[i]][,1],ob[[i]][,4],rep('02', length(ob[[i]][,4])),rep(object$levels[i], length(ob[[i]][,4]))),
              cbind(ob[[i]][,1],ob[[i]][,5],rep('11', length(ob[[i]][,5])),rep(object$levels[i], length(ob[[i]][,5]))),
              cbind(ob[[i]][,1],ob[[i]][,6],rep('12', length(ob[[i]][,6])),rep(object$levels[i], length(ob[[i]][,6]))))))



            if (ci == TRUE) {

              ob3<-obCI[[i]][c(1,3,5,7,9)]#[itpCI[ii]]

              ob3.2<-as.data.frame(rbind(ob3.2,rbind(cbind(ob3[,1],rep('00', length(ob3[,1])),rep(object$levels[i], length(ob3[,1]))),
                                                     cbind(ob3[,2],rep('01', length(ob3[,2])),rep(object$levels[i], length(ob3[,2]))),
                                                     cbind(ob3[,3],rep('02', length(ob3[,3])),rep(object$levels[i], length(ob3[,3]))),
                                                     cbind(ob3[,4],rep('11', length(ob3[,4])),rep(object$levels[i], length(ob3[,4]))),
                                                     cbind(ob3[,5],rep('12', length(ob3[,5])),rep(object$levels[i], length(ob3[,5]))))))



              ob4<-obCI[[i]][c(1,3,5,7,9)+1]#[itpCI[ii]+1]


              ob4.2<-as.data.frame(rbind(ob4.2,rbind(cbind(ob4[,1],rep('00', length(ob4[,1])),rep(object$levels[i], length(ob4[,1]))),
                                                     cbind(ob4[,2],rep('01', length(ob4[,2])),rep(object$levels[i], length(ob4[,2]))),
                                                     cbind(ob4[,3],rep('02', length(ob4[,3])),rep(object$levels[i], length(ob4[,3]))),
                                                     cbind(ob4[,4],rep('11', length(ob4[,4])),rep(object$levels[i], length(ob4[,4]))),
                                                     cbind(ob4[,5],rep('12', length(ob4[,5])),rep(object$levels[i], length(ob4[,5]))))))


            }


          }

          names(ob2)<-c('time','TP','trans','type')


          if (ci == TRUE) {

            names(ob3.2)<-c('tp_min','trans','type')


            names(ob4.2)<-c('tp_max','trans','type')


          }


          ob2f<-ob2[ob2$trans==trans,]



          if (ci == TRUE) {

            ob3.2f<-ob3.2[ob3.2$trans==trans,]


            ob4.2f<-ob4.2[ob4.2$trans==trans,]


          }

          if (ci == TRUE) {


            p1<-ggplot(ob2f, aes(x=as.numeric(ob2f$time), y=as.numeric(ob2f$TP), group=factor(ob2f$type),
                                 fill=factor(ob2f$type)))+theme_bw()+labs(x = xlab)+labs(y = ylab)



            p2<-p1+geom_ribbon(aes(ymin=as.numeric(ob3.2f$tp_min),ymax=as.numeric(ob4.2f$tp_max)),alpha=.3)

            p3<-p2+geom_line(aes(x=as.numeric(ob2f$time), y=as.numeric(ob2f$TP), color=ob2f$type), size=1.1)+
              theme(legend.title=element_blank())

            if(isTRUE(interactive)){
              if (requireNamespace("plotly", quietly=TRUE)) {return(plotly::ggplotly(p3))}
            }else{
              return(p3)
            }



          }else{

            p3<-ggplot(ob2f, aes(x=as.numeric(ob2f$time), y=as.numeric(ob2f$TP),group=factor(ob2f$type),
                                 color=factor(ob2f$type)))



            p4<-p3+theme_bw()+labs(x = xlab)+labs(y = ylab)+geom_line(size=1)+ theme(legend.title=element_blank())


            if(isTRUE(interactive)){
              if (requireNamespace("plotly", quietly=TRUE)) {return(plotly::ggplotly(p4))}
            }else{
              return(p4)
            }



          }


        }
      } #end plot for tp



      #cif
      #---------------------

      if (class(object)[1] %in%  c("CIF", "cifIPCW")) {

        if (is.null(ylab) &  class(object)[1] == "CIF") ylab <- "CIF(t)"

        if (is.null(ylab) &  class(object)[1] == "cifIPCW")
          ylab <- bquote(paste("CIF(t|", .(x$z.name), ")"))


        if(class(object)[1] == "cifIPCW") object$s <- 0


        if (object$Nlevels == 1) {

          if(object$s != 0){
            ob <- ob[, -2]
            #obCI <- obCI[, -c(1:2)]
          }


          if(class(object)[1] == "cifIPCW" & ci == TRUE) {
            obCI <- ob[, 3:4] # in order to corerct the out of cifIPCW
          }


          if (ci == TRUE) {

            if(object$s != 0){
              #ob <- ob[, -2]
              obCI <- obCI[, -c(1:2)]
            }



            p1<-ggplot(ob,aes(ob[,1],ob[,2]))+theme_bw()+labs(x=xlab,y = ylab)

            p2<-p1+geom_ribbon(aes(ymin=as.numeric(obCI[, 1]), ymax=as.numeric(obCI[, 2]),alpha=0.7),fill='gray')

            p3<-p2+geom_line(aes(ob[,1],ob[,2]))+geom_line(color='black',size=1)+#theme(legend.position="none")
              theme(legend.title=element_blank())

            if(isTRUE(interactive)){
              if (requireNamespace("plotly", quietly=TRUE)) {return(plotly::ggplotly(p3))}
            }else{
              return(p3)
            }

          }else{



            p1<-ggplot(ob,aes(ob[,1],ob[,2]))+theme_bw()+labs(x=xlab)+labs(y = ylab)

            p2<-p1+geom_line(aes(ob[,1],ob[,2]))+geom_line(color='red',size=1)



            if(isTRUE(interactive)){
              if (requireNamespace("plotly", quietly=TRUE)) {return(plotly::ggplotly(p2))}
            }else{
              p2
            }



          }

        }else{ # more than 1 level



          ob2<-NULL

          ob3.2<-NULL

          ob4.2<-NULL


          for(i in 1:object$Nlevels){

            #i<-1

            ob2<-as.data.frame(rbind(ob2,rbind(cbind(ob[[i]][,1],ob[[i]][,2],rep(object$levels[i], length(ob[[i]][,2]))))))


            if (ci == TRUE) {

              ob3<-obCI[[i]][1]


              ob3.2<-as.data.frame(rbind(ob3.2,rbind(cbind(ob3[,1],rep(object$levels[i], length(ob3[,1]))))))


              ob4<-obCI[[i]][2]

              ob4.2<-as.data.frame(rbind(ob4.2,rbind(cbind(ob4[,1],rep(object$levels[i], length(ob4[,1]))))))


            }

          }


          names(ob2)<-c('time','cif','type')

          head(ob2)

          if (ci == TRUE) {

            names(ob3.2)<-c('cif_min','type')


            names(ob4.2)<-c('cif_max','type')


            p1<-ggplot(ob2, aes(x=as.numeric(time), y=as.numeric(cif), group=factor(type),
                                fill=factor(type)))+theme_bw()+labs(x = xlab,y = ylab)

            p2<-p1+geom_ribbon(aes(ymin=as.numeric(ob3.2$cif_min),ymax=as.numeric(ob4.2$cif_max)),alpha=.3)


            p3<-p2+geom_line(aes(x=as.numeric(time), y=as.numeric(cif), color=ob2$type), size=1.1)+
              theme(legend.title=element_blank())


            if(isTRUE(interactive)){
              if (requireNamespace("plotly", quietly=TRUE)) {return(plotly::ggplotly(p3))}
            }else{
              return(p3)
            }



          }else{

            p3<-ggplot(ob2, aes(x=as.numeric(time), y=as.numeric(cif),group=factor(type),
                                color=factor(type)))


            p4<-p3+theme_bw()+labs(x = xlab,y = ylab)+geom_line(size=1)+ theme(legend.title=element_blank())

            if(isTRUE(interactive)){
              if (requireNamespace("plotly", quietly=TRUE)) {return(plotly::ggplotly(p4))}
            }else{
              return(p4)
            }



          }

        }#end more



      } #ends for CIF



      # soj
      #--------------------------------

      if (class(object)[1] %in%  c("soj", "sojIPCW")) {
        if (is.null(ylab) &  class(object)[1] == "soj") ylab <- "Sojourn(t)"

        if (is.null(ylab) &  class(object)[1] == "sojIPCW")
          ylab <- bquote(paste("Sojourn(t|", .(x$z.name), ")"))


        if (object$Nlevels == 1) {


          #   if(class(object)[1] == "cifIPCW") obCI <- ob[, 3:4] # in order to corerct the out of cifIPCW
          if(func == "survival"){ob[, 2] <- 1 - ob[, 2]}
          if(func == "survival" & ci == TRUE){obCI <- 1 - obCI}


          if (ci == TRUE) {


            p1<-ggplot(ob,aes(ob[,1],ob[,2]))+theme_bw()+labs(x=xlab,y = ylab)


            p2<-p1+geom_ribbon(aes(ymin=as.numeric(obCI[, 1]), ymax=as.numeric(obCI[, 2]),alpha=0.8),fill='gray')

            p3<-p2+geom_line(aes(ob[,1],ob[,2]))+geom_line(color='black',size=1)+theme(legend.position="none")


            if(isTRUE(interactive)){
              if (requireNamespace("plotly", quietly=TRUE)) {return(plotly::ggplotly(p3))}
            }else{
              return(p3)
            }



          }else{


            p1<-ggplot(ob,aes(ob[,1],ob[,2]))+theme_bw()+labs(x=xlab)+labs(y = ylab)

            p2<-p1+geom_line(aes(ob[,1],ob[,2]))+geom_line(color='red',size=1)

            if(isTRUE(interactive)){
              if (requireNamespace("plotly", quietly=TRUE)) {return(plotly::ggplotly(p2))}
            }else{
              return(p2)
            }


          }


        }else{ # more than 1 level


          ob2<-NULL

          ob3.2<-NULL

          ob4.2<-NULL


          for(i in 1:object$Nlevels){

            #i<-1

            ob2<-as.data.frame(rbind(ob2,rbind(cbind(ob[[i]][,1],ob[[i]][,2],rep(object$levels[i], length(ob[[i]][,2]))))))


            if (ci == TRUE) {

              ob3<-obCI[[i]][1]

              ob3.2<-as.data.frame(rbind(ob3.2,rbind(cbind(ob3[,1],rep(object$levels[i], length(ob3[,1]))))))

              ob4<-obCI[[i]][2]

              ob4.2<-as.data.frame(rbind(ob4.2,rbind(cbind(ob4[,1],rep(object$levels[i], length(ob4[,1]))))))


            }

          }


          names(ob2)<-c('time','sojourn','type')


          if (ci == TRUE) {

            names(ob3.2)<-c('sojourn_min','type')

            names(ob4.2)<-c('sojourn_max','type')

            p1<-ggplot(ob2, aes(x=as.numeric(time), y=as.numeric(sojourn), group=factor(type),
                                fill=factor(type)))+theme_bw()+labs(x = xlab,y = ylab)

            p2<-p1+geom_ribbon(aes(ymin=as.numeric(ob3.2$sojourn_min),ymax=as.numeric(ob4.2$sojourn_max)),alpha=.3)


            p3<-p2+geom_line(aes(x=as.numeric(time), y=as.numeric(sojourn), color=type), size=1.1)+
              theme(legend.title=element_blank())


            if(isTRUE(interactive)){
              if (requireNamespace("plotly", quietly=TRUE)) {return(plotly::ggplotly(p3))}
            }else{
              return(p3)
            }



          }else{

            p3<-ggplot(ob2, aes(x=as.numeric(time), y=as.numeric(sojourn),group=factor(type),
                                color=factor(type)))


            p4<-p3+theme_bw()+labs(x = xlab,y = ylab)+geom_line(size=1)+ theme(legend.title=element_blank())

            if(isTRUE(interactive)){
              if (requireNamespace("plotly", quietly=TRUE)) {return(plotly::ggplotly(p4))}
            }else{
              return(p4)
            }




          }

        }#end more



      } #ends for sojourn



    }else{
      stop("Argument x must be either survIDM object.")
    }

  }

}


