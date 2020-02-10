Es<-function(r,qI,qO,qC,N){
a<-qC/(1-qC/N)
z<-qI*(1-r)*(N-1)
x<-(a+qO*r+qI*(1-r))
y<-(qO*r+qI*(1-r)-(r*a*qO)/(a+N*qI*(1-r)))
ans<-y/x
return(ans)
}

P0<-function(r,qI,qO,qC,N){
A<-r*qO
B<-(1-r)*qI
ans<-qC/(qC+(N-qC)*B)
return(ans)
}

r_Pareto<-function(qI,qO,qC,N){
a<-qC/(1-qC/N)
qIc<-((N-1)*qO-a)/N
if(qI<qIc){
X<-sqrt((a+qO)*(qO-qI)*(N-1))
Y<-sqrt(N*qO*(a+N*qI))
ans<-(1-(-a*X+(qO+a)*Y)/(N*qI*X+(qO-qI)*Y))
}
else{
ans<-0
}
return(ans)
}

Es_Pareto<-function(qI,qO,qC,N){
r<-r_Pareto(qI,qO,qC,N)
ans<-Es(r,qI,qO,qC,N)
return(ans)
}

r_Nash<-function(qI,qO,qC,N){
a<-qC/(1-qC/N)
D2<-((a+qO)*(qO-qI*N)+(qO-qI)*(a*N+qO))**2-4*(qO-qI)*(qO-qI*N*N)*(a+qO)*(a+qO)
qIc<-((N-1)*qO-a)/N
if(qI<qIc){
ans<-1-2*(a+qO)*(a+qO)/((a+qO)*(qO-qI*N)+(qO-qI)*(a*N+qO)+sqrt(D2))
}
else{
ans<-0
}
return(ans)
}

Es_Nash<-function(qI,qO,qC,N){
r<-r_Nash(qI,qO,qC,N)
ans<-Es(r,qI,qO,qC,N)
return(ans)
}

Es_I<-function(qI,qO,qC,N){
a<-qC/(1-qC/N)
ans<-qI/(a+qI)
return(ans)
}

Es_S<-function(qI,qO,qC,N,NI){
a<-qC/(1-qC/N)
ans<-NI*qI*qO/((a+qO)*(a+NI*qI))
return(ans)
}

r_optim<-function(ro,qI,qO,qC,N){
r<-qI*(N-1)*(1-ro)
a<-qC/(1-qC/N)
x<-(-qI*(a+qI+r)+sqrt(qI*qO*r)*sqrt(a+qI+r))
y<-qI*(qO-qI)
ans<-x/y
if(ans>1){ans<-1}
if(ans<0){ans<-0}
return(ans)
}

wn<-function(r,ro,qI,qO,qC,N){
a<-qC/(1-qC/N)
x<-(a+qO*r+qI*(1-r))
A<-(a*qO*r)/(a+N*qI-qI*((N-1)*ro+r))
y<-(qO*r+qI*(1-r)-A)
ans<-y/x
return(ans)
}

dwn<-function(r,ro,qI,qO,qC,N){
a<-qC/(1-qC/N)
a2<-a*(qO-qI)
a3<-a*qO*(a+qI)
b<-(a+qO*r+qI*(1-r))
d<-(a+N*qI-qI*((N-1)*ro+r))
ans<-a2/(b*b)-a3/(b*b*d)-(a*qO*r*qI)/(b*d*d)
if(r==0 && ans< 0) ans<-0
if(r==1 && ans>0)  ans<-0
return(ans)
}


W<-function(r,qI,qO,qC){
N<-length(r)
ans<-0
for(i in 1:N){
ro<-(sum(r)-r[i])/(N-1)
ans<-ans+wn(r[i],ro,qI,qO,qC,N)
}
return(ans)
}


dWn<-function(n,r,qI,qO,qC){
N<-length(r)
ro<-(sum(r)-r[n])/(N-1)
ans<-dwn(r[n],ro,qI,qO,qC,N)
a<-qC/(1-qC/N)
for(i in 1:N){
if(i!=n){
b<-(a+qO*r[i]+qI*(1-r[i]))
d<-(a+N*qI-qI*sum(r))
ans<-ans-a*qI*qO*r[i]/(b*d*d)
}
}
if(r[n]==0 && ans< 0) ans<- 0
if(r[n]>=0.99999 && ans>0)  ans<- 0
return(ans)
}

