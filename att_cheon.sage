def encodeRand(rho,p):
  return crt([ZZ.random_element(2^rho) for i in range(len(p))],p)

def AtkCheon():
  eta=100
  n=10
  rho=20
  p=[random_prime(2^eta,False,2^(eta-1)) for i in range(n)]
  x0=prod(p)
  pzt=crt([ZZ.random_element(2^rho)*x0/pi for pi in p],p)
  c=[encodeRand(rho,p) for i in range(n)]
  d=[encodeRand(rho,p) for i in range(n)]
  b0=encodeRand(rho,p)
  b1=encodeRand(rho,p)

  W0=matrix(ZZ,[[ci*dj*b0*pzt % x0 for dj in d] for ci in c])
  W1=matrix(ZZ,[[ci*dj*b1*pzt % x0 for dj in d] for ci in c])
  
  print "Basic attack",
  W=W0*W1^-1
  rec_primes=sorted([gcd(x0,a.denominator()*b0-a.numerator()*b1) 
                     for a in W.eigenvalues()])
  assert rec_primes==sorted(p)
  print "OK"

  print "Attack modulo q",
  q=random_prime(2^eta,False,2^(eta-1))
  W0q=matrix(Integers(q),W0)
  W1q=matrix(Integers(q),W1)
    
  Wq=W0q*W1q^-1
  eigen=[a.rational_reconstruction() for a in Wq.eigenvalues(extend=False)]
  rec_primes=sorted([gcd(x0,a.denominator()*b0-a.numerator()*b1) for a in eigen])
  assert rec_primes==sorted(p)
  print "OK"

def encodeVec(m,rho,p):
  return vector([encodeRand(rho,p) for i in range(m)])

def AtkMatrixCheon(m=2):
  eta=100
  n=5
  rho=20
  p=[random_prime(2^eta,False,2^(eta-1)) for i in range(n)]
  x0=prod(p)
  pzt=crt([ZZ.random_element(2^rho)*x0/pi for pi in p],p)
  d=m*n
  c=[encodeVec(m,rho,p) for i in range(d)]
  d=[encodeVec(m,rho,p) for i in range(d)]
  b0=matrix([encodeVec(m,rho,p) for i in range(m)])
  b1=matrix([encodeVec(m,rho,p) for i in range(m)])

  W0=matrix(ZZ,[[ci*b0*dj*pzt % x0 for dj in d] for ci in c])
  W1=matrix(ZZ,[[ci*b1*dj*pzt % x0 for dj in d] for ci in c])
  
  print "Basic attack",
  W=W0*W1^-1
  f=W.charpoly()
  B=matrix(Integers(x0),b0)*matrix(Integers(x0),b1)^-1
  rec_primes=sorted([gcd(f1[0].change_ring(Integers(x0))(B)[0,0],x0) 
                     for f1 in factor(f)])
  assert rec_primes==sorted(p)
  print "OK"

  print "Variant attack",
  W=W0*W1^-1
  f=W.charpoly()
  rec_primes = sorted([gcd(matrix(Integers(x0),f1[0](W)*W0)[0,0],x0) 
                       for f1 in factor(f)])
  assert rec_primes==sorted(p)
  print "OK"

def encodeMat(m,rho,p):
  return matrix(ZZ,[[encodeRand(rho,p) for i in range(m)] for j in range(m)])

def AtkTensoringCheon(m=2):
  eta=180
  n=3
  rho=20
  p=[random_prime(2^eta,False,2^(eta-1)) for i in range(n)]
  x0=prod(p)
  pzt=crt([ZZ.random_element(2^rho)*x0/pi for pi in p],p)
  d=m^3*n
  a1=[encodeVec(m,rho,p) for i in range(d)]
  B1=[encodeMat(m,rho,p) for i in range(2)]
  C1=[encodeMat(m,rho,p) for i in range(d)]
  A2=[encodeMat(m,rho,p) for i in range(d)]
  B2=[encodeMat(m,rho,p) for i in range(2)]
  c2=[encodeVec(m,rho,p) for i in range(d)]
  
  W=[matrix(ZZ,[[a1[i]*B1[k]*C1[j]*A2[i]*B2[k]*c2[j]*pzt % x0 
                 for j in range(d)] for i in range(d)]) 
     for k in range(2)]

  WW=W[0]*W[1]^-1
  f=WW.charpoly()

  print "Basic attack",
  BB=[matrix(Integers(x0),
      B2[k].tensor_product(matrix.identity(m)).tensor_product(B1[k])) 
      for k in range(2)]
  BinvB=BB[0]*BB[1]^-1
  rec_primes=sorted([gcd(f1[0].change_ring(Integers(x0))(BinvB)[0,0],x0) 
                     for f1 in factor(f)])
  assert rec_primes==sorted(p)
  print "OK"

  print "Variant attack",
  rec_primes = sorted([gcd(matrix(Integers(x0),f1[0](WW)*W[0])[0,0],x0) 
                       for f1 in factor(f)])
  assert rec_primes==sorted(p)
  print "OK"
