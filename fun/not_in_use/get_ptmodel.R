get_ptmodel<-function(location,x,ptpars,pt_in, prm, ref, sel, agg, gps, hivpoints, sfin){

id=prm$id

objitv<- function(pt_itv) run_model(ptpars,x,pt_itv, prm, ref, sel, agg, gps, hivpoints, sfin)
objbase<-function(pt_itv) run_model(c(x[id], 0),x,pt_itv, prm, ref, sel, agg, gps, hivpoints, sfin)



#Baseline
pt_itv0<-c(9,0,0.05,0.69)
base<-objbase(pt_itv0)

#Intervention
intv<-objitv(pt_in)

return(intv)

}