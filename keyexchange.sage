import sys

##################################################
####    Defining parameters
##################################################

PARAM_SET = "small"
#PARAM_SET = "medium"
#PARAM_SET = "large"

# eta  : bit-length of each prime p_i
# m    : matrix size
# n    : number of primes 
# mu   : number of steps per each user in each repetition
# alpha: maximum bit-length of the g_i
# k    : number of repetitions
# rhom : noise of encoding in the matrices
# nu   : number of bits to extract
# Np   : number of parties
# tau  : number of matrices per step is 2^tau
# beta : bit length of the randomness in the zero-testing parameter

Np=4
rhom=2
k=2
tau=3

if "small" == PARAM_SET:
    sec_param=52 ; eta = 1759 ; m=6 ; n=160 ; mu = 15 ; alpha=11 ; nu = 62 
elif "medium" == PARAM_SET:
    sec_param = 62 ; eta=2602 ; m=6 ; n=294 ; mu=21; alpha=12 ; nu=72
elif "large" == PARAM_SET:
    sec_param = 72 ; eta=3761 ; m=6 ; n=1349 ; mu=27 ; alpha=14 ; nu=82

beta=alpha
m0=m/3 ; m1=m/3     # length of submatrices
    
from time import gmtime, strftime

def print_time(msg):
    print strftime(msg + " (%a, %d %b %Y %H:%M:%S)", gmtime())
    sys.stdout.flush()

def symmetric_mod(x, p):
    x=ZZ(x) % p
    if x>p//2:
        return x-p
    return x

# Returns an invertible matrix modulo gi
def sample_matrix_slot(_m,gi):
    M = Matrix(ZZ, _m, _m) 
    while True:
        for i in range(_m):
            for j in range(_m):
                M[i,j] = ZZ.random_element(gi)
        if det(M) !=0:
            return M

# Returns an invertible matrix modulo the product of the gi's
def sample_matrix(_m,g):
    M=Matrix(ZZ,_m,_m)
    Mg=[sample_matrix_slot(_m,gi) for gi in g]
    for i in range(_m):
        for j in range(_m):
            M[i,j]=crt([Mgi[i,j] for Mgi in Mg],g)
    return M
                
def sample_matrix_A_with_block_B(B,g):
    return block_matrix(ZZ, [[sample_matrix(m0,g), 0, 0], [0, sample_matrix(m0,g), 0], [0, 0, B]])

# we return an encoding element modulo x0
def encode(value,g,p,x0,rho):
    values = [ZZ.random_element(-2^rho,2^rho)*g[i] + (value % g[i]) for i in range(len(p))]
    return mod(crt(values, p),x0)

# We encode a matrix, with noise rhom
def encode_mat(A,g,p,x0):
    C = Matrix(Integers(x0), A.nrows(), A.ncols())
    values=[[ZZ.random_element(-2^rhom,2^rhom)*gk + (A[i,j] % gk) for i in range(A.nrows()) for j in range(A.ncols())]
            for gk in g]
    enc_values=CRT_vectors(values,p)
    for i in xrange(A.nrows()):
        for j in xrange(A.ncols()):
            C[i,j] = enc_values[i*A.ncols()+j]
    return C


def sample_K_and_invK(m, x0):
    MK=random_matrix(Integers(x0),m)
    return MK,MK^-1

def sample_bookend_vectors(s_asterisk,t_asterisk,g,p,x0):
    pg=prod(g)
    rhob=alpha
    tmp = [0 for i in xrange(m0)] + [ZZ.random_element(pg) for i in xrange(m0)] + s_asterisk
    s=vector(Integers(x0),[encode(tmpi,g,p,x0,rhob) for tmpi in tmp])
    tmp = [ZZ.random_element(pg) for i in xrange(m0)] + [0 for i in xrange(m0)] +  t_asterisk
    t=vector(Integers(x0),[encode(tmpi,g,p,x0,rhob) for tmpi in tmp])
    return s,t

# the bundling scalars are all +-1 modulo each gi
def sample_pm1_g(g):
    return crt([(-1)^ZZ.random_element(2) for i in range(len(g))],g)

# sample the bundling scalars
def sample_alphas(g):
    kappa = mu * Np * k # degree of multilinearity
    pg=prod(g)
    alphas=[[[sample_pm1_g(g) for b in range(2^tau)] for i in range(kappa)] for u in range(Np)]
    for u in range(1,Np):
        for b in range(2^tau):
            for i in range(mu*Np):
                alphas[u][i][b]=(mod(prod(alphas[0][Np*mu*j+i][b] for j in range(k)),pg)/mod(prod(alphas[u][Np*mu*j+i][b] for j in range(1,k)),pg)).lift()
    return alphas

def genpseudoprime(eta,etamin=211):  # etamin=211 for lam=52, etamin=268 for lam=62
  if eta<=etamin:
    raise IOError
  if eta<=(2*etamin):
    return random_prime(2^eta,False,2^(eta-1))
  else:
    return random_prime(2^etamin,False,2^(etamin-1))*genpseudoprime(eta-etamin)

def runtimep(Np,kappa,tau,n,eta):
    return N(Np*kappa*2^tau*n^1.58*eta)

def test():
    tim=cputime()
    print_time("Generating parameters with n=%d, mu=%d, tau=%d ..." % (n,mu,tau))
    base_memory = get_memory_usage()
    print "       Memory usage (MB): 0"

    p = [genpseudoprime(eta) for i in range(n)] # CLT primes

    # the primes gi must be distinct
    P=Primes()
    g = [P.unrank(2^alpha//alpha+i) for i in range(n)] # g1, ..., gn that define the CLT plaintext space
    inv_g = [inverse_mod(g[i], p[i]) for i in xrange(n)]

    x0 = prod(p)

    # we don't publish pzt
    h = [ZZ.random_element(2^beta) for i in xrange(n)]
    pzt=crt([h[i]*inv_g[i]*x0/p[i] for i in range(n)],p)

    pg=prod(g)
    
    kappa = mu * Np * k # degree of multilinearity
    
    print "Predicted noise:",int(4*alpha+kappa*(alpha+rhom+.3*log(m,2))+beta-alpha)
    print "Predicted runtime parameter generation:",runtimep(Np,kappa,tau,n,eta)/runtimep(4,16,3,40,1759)*473
    
    # these are the same for all users
    s_asterisk = [ZZ.random_element(pg) for i in xrange(m1)]
    t_asterisk = [ZZ.random_element(pg) for i in xrange(m1)]
    B_matrices = [[sample_matrix(m1,g) for b in range(2^tau)] for i in range(kappa)]

    alphas=sample_alphas(g)  # these are the bundling scalars
    
    # matrices[u][i][b] = row u, i-th matrix, corresponding to bit b
    matrices = [[[ 0 for b in xrange(2^tau) ] for i in xrange(kappa)] for u in xrange(Np)]
    bookends_left = []
    bookends_right = []

    for u in xrange(Np):
        print_time("   starting generation of matrices of row " + str(u))

        # for each party
        Ki, invKi = sample_K_and_invK(m, x0)
        s, t = sample_bookend_vectors(s_asterisk,t_asterisk,g,p,x0)
        s = s*invKi
        for i in xrange(kappa):
            Ki_plus_1, invKi_plus_i = sample_K_and_invK(m, x0)
            for b in range(2^tau):
                B = alphas[u][i][b]*B_matrices[i][b]
                A = sample_matrix_A_with_block_B(B,g)
                matrices[u][i][b] = Ki * encode_mat(A,g,p,x0) * invKi_plus_i
            Ki = Ki_plus_1

        t = pzt * Ki * t
        bookends_left.append(s)
        bookends_right.append(t)
        print "       Memory usage (MB):", get_memory_usage() - base_memory
    print "Time parameter generation:", cputime(tim)
    print
        
    tim=cputime()
    ####### Publishing phase: each party u publishes the shared information based in his secret

    # users' secret bit strings
    secrets = [[ZZ.random_element(2^tau) for i in range(mu)] for u in range(Np)]

    # public_products[u][v][j] = j-th product from user u to user v (using matrices from row v and secret of u, j-th repetition)
    public_products = [[[0 for j in range(k)] for v in xrange(Np)] for u in xrange(Np)]
    
    for u in range(Np):
        print_time("Party " + str(u) + " publishing...")
        for v in range(Np):
            if u != v:
                for j in range(k):
                    public_products[u][v][j] = prod(matrices[v][mu*Np*j+mu*u+i][secrets[u][i]] for i in range(mu))
                
    print "Memory usage (MB):", get_memory_usage() - base_memory
    print "Time publishing phase generation:", cputime(tim)
    print

    tim=cputime()
    
    ####### Key derivation phase
    shared_keys = []
    prolist=[]
    for u in range(Np):
        print_time("Starting to derive key of party " + str(u))
        v = u
        # party u computes his products based on his secret
        for j in xrange(k):
             public_products[u][v][j] = prod(matrices[v][mu*Np*j+mu*u+i][secrets[u][i]] for i in range(mu))
    
             # party u compute his whole product based in all the public products
        pro = bookends_left[u]*prod(public_products[v][u][j] for j in range(k) for v in range(Np))*bookends_right[u]
        prolist.append(pro)

        # party u extracts the shared key
        shared_keys.append(pro.lift() >> (pro.lift().nbits() - nu))

        print_time("Finished deriving key of party " + str(u))

    print "shared_keys =", shared_keys
    print "Time key derivation:", cputime(tim)

    print "Noise:",symmetric_mod(prolist[0]-prolist[1],x0).nbits()+eta-x0.nbits()
    

    print "kappa=",kappa
    print

def runtest():
    global n,mu,tau,eta
    #n=20
    #mu=2
    test()
