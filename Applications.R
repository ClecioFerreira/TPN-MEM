
source('functions.R')
source('mem_normal.r')
require(mnormt)

# Barnett data set (1969)
z=as.matrix(read.table('barnett.txt',header=F),72,4)
z=z/100

# TPN-MEM model
theta=EM.NTP_MEM(z)
MI=MIapprox2(theta,z)
ep=sqrt(diag(solve(MI)))
logntp<-logNTP_MEM(theta,z)
d=3
n=72
p1=3*d+4
AIC=-2*logntp+2*p1
BIC=-2*logntp+p1*log(n)

# Normal-MEM model
thetan=mem.em.normal(z,1e-5)
logn=veron(thetan,z)
MI=MIapproxn(thetan,z)
ep=sqrt(diag(solve(MI)))
p0=3*d+3
AIC=-2*logn+2*p0
BIC=-2*logn+p0*log(n)


###################################################################
#  systolic blood pressure data set

z=as.matrix(read.table('testiculo.txt',header=F),42,5)
z=z[,2:6]

# TPN-MEM model
theta=EM.NTP_MEM(z)
MI=MIapprox2(theta,z)
ep=sqrt(diag(solve(MI)))
logntp<-logNTP_MEM(theta,z)
d=4
n=42
p1=3*d+4
AIC=-2*logntp+2*p1
BIC=-2*logntp+p1*log(n)


# Normal-MEM Model
thetan=mem.em.normal(z,1e-5)
logn=veron(thetan,z)
MI=MIapproxn(thetan,z)
ep=sqrt(diag(solve(MI)))
p0=3*d+3
AIC=-2*logn+2*p0
BIC=-2*logn+p0*log(n)

