def AttOrtho(eta=50,n=30,rho=10):
  p=[random_prime(2^eta,False,15*2^(eta-4)) for i in range(n)]
  x0=prod(p)
  r=[[ZZ.random_element(-2^rho+1,2^rho) for i in range(n)] 
     for j in range(n)]
  x=[crt(rj,p) for rj in r]
  y=crt([ZZ.random_element(-2^rho+1,2^rho) for i in range(n)],p)

  x=x+[mod(xi*y,x0).lift() for xi in x]

  tau=2*n
  M=matrix(ZZ,tau,tau)
  for i in range(tau-1):
    M[i,i]=1
    M[i,tau-1]=mod(-x[i]*inverse_mod(x[tau-1],x0),x0)
  M[tau-1,tau-1]=x0

  ML=M.LLL()
  V=ML[:tau-n].right_kernel().matrix()
  W0,W1=V[:,:n].T,V[:,n:tau].T
  v=(W1*W0^-1).eigenvalues()
  rprimes=Set([gcd(y-vi,x0) for vi in v])
  print "Number of primes recovered:",len(Set(p).intersection(rprimes)),
  print "out of",n
  
def smallMat(nrows,ncols,rho,p):
  n=len(p)
  r=[Matrix([[ZZ.random_element(-2^rho+1,2^rho) for k in range(ncols)] 
             for j in range(nrows)]) for i in range(n)]
  return Matrix(ZZ,[[crt([r[i][j,k] for i in range(n)],p) 
                   for k in range(ncols)] 
                    for j in range(nrows)])
      
def AttOrthoVec(eta=60,n=5,rho=10,m=2):
  p=[random_prime(2^eta,False,2^(eta-1)) for i in range(n)]
  x0=prod(p)
  V=smallMat(n*m,m,rho,p)
  K=random_matrix(Integers(x0),m)
  VT=V*K
  Kp=random_matrix(Integers(x0),m)
  A0,A1=smallMat(m,m,rho,p),smallMat(m,m,rho,p)
  C0,C1=K^-1*A0*Kp,K^-1*A1*Kp
  Vp0,Vp1=VT*C0,VT*C1

  Vp=Matrix(ZZ,2*n*m,m)
  Vp[:n*m,:],Vp[n*m:,:]=Vp0,Vp1

  tau=2*n*m
  M=matrix(ZZ,tau,tau)
  for i in range(tau-m):
    M[i,i]=1
  M[:,tau-m:]=-Vp*Matrix(Integers(x0),Vp[tau-m:,:]).inverse()
  M[tau-m:,tau-m:]=x0*matrix.identity(m)

  ML=M.LLL()

  VO=ML[:m*n].right_kernel().matrix()
  W0,W1=VO[:,:n*m].T,VO[:,n*m:].T
  f=(W0*W1^-1).charpoly()
  B=matrix(Integers(x0),C0)*matrix(Integers(x0),C1)^-1
  rprimes=Set([gcd(f1[0].change_ring(Integers(x0))(B)[0,0],x0) 
               for f1 in factor(f)])
  print "Number of primes recovered:",
  print len(Set(p).intersection(rprimes)),"out of",n
