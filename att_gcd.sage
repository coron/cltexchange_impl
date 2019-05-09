import sage.rings.polynomial.polynomial_ring as pring

def genXi(rho,eta,gam,p):
  return p*ZZ.random_element(2^(gam-eta))+ZZ.random_element(2^rho)

def multiPointEval(f,li,R):
  n,x=len(li),R.0
  if n==1: return [f(x=li[0])]
  li1,li2=li[:n//2],li[n//2:]
  f1=f.quo_rem(prod((x-xi for xi in li1)))[1]
  f2=f.quo_rem(prod((x-xi for xi in li2)))[1]
  return multiPointEval(f1,li1,R)+multiPointEval(f2,li2,R)

def attackGCDMultiPoint(rho=12,eta=1000,gam=40000,verbose=True):
  p=random_prime(2^eta,lbound=2^(eta-1),proof=False)
  if verbose: print p
  x0=p*prod([random_prime(2^eta,lbound=2^(eta-1),proof=False) 
             for i in range(gam/eta-1)])
  R=pring.PolynomialRing_dense_mod_n(Integers(x0),'x')
  x=R.0
  c=genXi(rho,eta,gam,p)

  t=cputime(subprocesses=True)
  f=prod((x+i-c for i in range(2^(rho/2))))
  ev=multiPointEval(f,range(0,2^rho,2^(rho/2)),R)
  pp=gcd(x0,prod(ev))
  if verbose: print pp
  return cputime(t)

def attackGCDMaskedMultiPoint(rho=20,eta=100,gam=1000,verbose=True):
  p=random_prime(2^eta,lbound=2^(eta-1),proof=False)
  print "p=",p
  x0=p*prod([random_prime(2^eta,lbound=2^(eta-1),proof=False) 
             for i in range(gam/eta-1)])
  R=pring.PolynomialRing_dense_mod_n(Integers(x0),'x')
  x=R.0
  z=Integers(x0).random_element()

  L1=[genXi(rho,eta,gam,p)*z for i in range(2^(rho//2+1))]
  L2=[genXi(rho,eta,gam,p)*z for i in range(2^(rho//2+1))]

  t=cputime(subprocesses=True)
  f=prod((x-c for c in L1))
  ev=multiPointEval(f,L2,R)
  pp=gcd(x0,prod(ev)).lift()
  print pp
  return cputime(t)

def genVecXi(rho,eta,gam,m,p,K,x0):
  v=vector([genXi(rho,eta,gam,p) for i in range(m)])
  return K*v

def attackGCDVectorMultiPoint(rho=10,eta=100,gam=1000,m=2):
  p=random_prime(2^eta,lbound=2^(eta-1),proof=False)
  print "p=",p
  x0=p*prod([random_prime(2^eta,lbound=2^(eta-1),proof=False) 
  	     for i in range(gam/eta-1)])

  R=pring.PolynomialRing_dense_mod_n(Integers(x0),'x')	     
  x=R.0
  
  K=matrix(Integers(x0),m,m)
  for i in range(m):
    for j in range(m):
      K[i,j]=ZZ.random_element(x0)
      
  L1=[genVecXi(rho,eta,gam,m,p,K,x0)[0] for i in range(2^(m*rho//2+1))]
  L2=[genVecXi(rho,eta,gam,m,p,K,x0)[0] for i in range(2^(m*rho//2+1))]

  t=cputime(subprocesses=True)
  f=prod((x-c for c in L1))
  ev=multiPointEval(f,L2,R)
  pp=gcd(x0,prod(ev)).lift()
  print pp
  return cputime(t)

def genMatXi(rho,eta,gam,m,p,K,K2,x0):
  M=Matrix([[genXi(rho,eta,gam,p) for i in range(m)] for j in range(m)])
  return K*M*K2

def attackGCDMatrixMultiPoint(rho=4,eta=100,gam=200,m=2,verbose=True):
  p=random_prime(2^eta,lbound=2^(eta-1),proof=False)
  print "p=",p
  x0=p*prod([random_prime(2^eta,lbound=2^(eta-1),proof=False) for i in range(gam/eta-1)])

  R=sage.rings.polynomial.polynomial_ring.PolynomialRing_dense_mod_n(Integers(x0),'x')
  x=R.0
  
  K=matrix(Integers(x0),m,m)
  K2=matrix(Integers(x0),m,m)
  for i in range(m):
    for j in range(m):
      K[i,j]=ZZ.random_element(x0)
      K2[i,j]=ZZ.random_element(x0)
      
  L1=[genMatXi(rho,eta,gam,m,p,K,K2,x0)[0,0] for i in range(2^(m*m*rho//2+1))]
  L2=[genMatXi(rho,eta,gam,m,p,K,K2,x0)[0,0] for i in range(2^(m*m*rho//2+1))]

  t=cputime(subprocesses=True)
  f=prod((x-c for c in L1))
  ev=multiPointEval(f,L2,R)
  pp=gcd(x0,prod(ev)).lift()
  print pp
  return cputime(t)