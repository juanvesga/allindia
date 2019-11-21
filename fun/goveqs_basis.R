
goveqs_basis <- function (t, state, parameters) {
  
  
  with(as.list(c(state,parameters)),             
       {
         
         
         invec <- as.numeric(state[1:i_all$nstates])
         mat<-M$lambda
         mat[,s$lo] <- mat[,s$lo]/sum(invec[s$lo]); 
         mat[,s$hi] <- mat[,s$hi]/sum(invec[s$hi]); 
         n_aux=i$nx-i_all$nstates
         
         dx<-compute_dx(invec,
                        mat,
                        s$lo-1,
                        s$hi-1,
                        M$lin,
                        M$nlin$lo$ds,
                        M$nlin$lo$mdr,
                        M$nlin$hi$ds,
                        M$nlin$hi$mdr,
                        as.numeric(M$mortvec),
                        M$xp,
                        M$sm,
                        M$xr,
                        M$mac,
                        p$growth, 
                        c(i$U$lo,i$U$hi,i_all$S$lo,i_all$S$hi)-1,
                        s$nstates -1,
                        s$respSymp -1,
                        p$resp_symptomatic,
                        p$hi,
                        agg$inc,
                        agg$notif,
                        agg$mort,
                        agg$pt,
                        agg$txfp,
                        agg$dx,
                        agg$acf,
                        agg$pmo,
                        sel$inc,
                        sel$acqu,
                        sel$remo,
                        sel$notif,
                        sel$pt,
                        sel$txfp,
                        sel$dx,
                        sel$acf,
                        i$aux$inc-1,
                        i$aux$notif-1,
                        i$aux$mort-1,
                        i$aux$remo-1,
                        i$aux$pt-1,
                        i$aux$txfp-1,
                        i$aux$dx-1,
                        i$aux$pmo-1,
                        i$aux$acf-1,
                        n_aux,
                        t)
         list(dx)
         
         
       }
  )
}