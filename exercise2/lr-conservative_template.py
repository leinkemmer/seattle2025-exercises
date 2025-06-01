from pylab import *
from numpy.linalg import solve
from numpy.linalg import qr

problem = 'ts'
nx = 128
nv = 128
r = 20

t_final = 40.0
deltat  = 0.01

if problem == 'll':
    L = 4*pi # linear Landau damping
else:
    L = 10*pi


xs = linspace(0, L, nx, endpoint=False)
hx = xs[1]-xs[0]

vs = linspace(-6.0, 6.0, nv, endpoint=False)
hv = vs[1]-vs[0]

#
# Differentiation matrices
#
A_cd_x = diag([1.0/(2.0*hx)]*(nx-1), 1) + diag([-1.0/(2.0*hx)]*(nx-1), -1)
A_cd_x[0, -1] = -1.0/(2.0*hx)
A_cd_x[-1, 0] =  1.0/(2.0*hx)

A_cd_v = diag([1.0/(2.0*hv)]*(nv-1), 1) + diag([-1.0/(2.0*hv)]*(nv-1), -1)
A_cd_v[0, -1] = -1.0/(2.0*hv)
A_cd_v[-1, 0] =  1.0/(2.0*hv)


A_v = diag(vs)

#
# Compute the coefficients
#
def compute_c1(V):
    # TODO

def compute_c2(V):
    # TODO

def compute_d1(X, E):
    # TODO

def compute_d2(X):
    # TODO

#
# RHS of the equations for K, S, L
#

# rhs of equation (2.13)
def rhs_K(K, E, c1, c2):
    # TODO

def rhs_K_withE(K, V, c1, c2):
    # TODO

# rhs of equation (2.14)
def rhs_S(S, c1, c2, d1, d2):
    # TODO

# rhs of equation (2.15)
def rhs_L(L, E, d1, d2):
    # TODO

#
# Computation of the electric field
#
def compute_rho(K, V):
    # TODO

def compute_mass(X, S, V):
    K = X @ S
    return hx*ones(nx).transpose() @ compute_rho(K, V)

def compute_E(K, V):
    rhohat = fft(compute_rho(K, V))
    Ehat = 1j*zeros(len(rhohat))
    Ehat[1:] = 1.0/(1j*2*pi/L*fftfreq(len(rhohat), 1)[1:]*len(rhohat))*rhohat[1:]
    Ehat[0] = 0.0
    return real(ifft(Ehat))

def electric_energy(E):
    return 0.5*sum(E**2)*hx

#
# Projector splitting integrator
#
def rk4(deltat, U, rhs):
    k1 = rhs(U)
    k2 = rhs(U + 0.5*deltat*k1)
    k3 = rhs(U + 0.5*deltat*k2)
    k4 = rhs(U + deltat*k3)
    return U + 1.0/6.0*deltat*(k1 + 2.0*k2 + 2.0*k3 + k4)


def time_step_augBUG(X, S, V):
    # K step
    # TODO
    
    # L step
    # TODO

    # setup the augmented basis
    Xa = zeros((nx,2*r))
    Xa[:,0:r] = X
    Xa[:,r:] = K
    Xa, _ = qr(Xa, mode='reduced')
    Xa *= 1.0/sqrt(hx)
    Xa[:,0:r] = X
    
    Va = zeros((nv,2*r))
    Va[:,0:r] = V
    Va[:,r:] = L
    Va, _ = qr(Va, mode='reduced')
    Va *= 1.0/sqrt(hv)
    Va[:,0:r] = V

    Sa = zeros((2*r,2*r))
    Sa[0:r,0:r] = S

    # S step
    c1 = compute_c1(Va)
    c2 = compute_c2(Va)
    d1 = compute_d1(Xa,E)
    d2 = compute_d2(Xa)
    # TODO

    # truncate
    # TODO: use an SVD to truncate back to rank r

    return X1, S1, V1


def conservative_truncate(X, S, V):
    cons_trunc = 1
    r_augm = S.shape[1]

    # Step 5
    K_tilde = X @ S
    K_constt = K_tilde[:nx, :cons_trunc]
    K_remt = K_tilde[:nx, cons_trunc:r_augm]

    # Step 6
    Q_cons, R_cons = qr(K_constt, mode='reduced')
    X_cons = Q_cons / np.sqrt(hx)
    S_cons = R_cons * np.sqrt(hx)

    # Step 7
    Q_remt, R_remt = qr(K_remt, mode='reduced')
    X_remt = Q_remt / np.sqrt(hx)
    S_remt = R_remt * np.sqrt(hx)

    # Step 8
    U, s, Vh = svd(S_remt, full_matrices=False)
    S_rem = np.diag(s[:r - cons_trunc])
    U_hat = U[:, :r - cons_trunc]
    W_hat = Vh[:r - cons_trunc, :].T

    X_rem = X_remt @ U_hat
    W_np1 = V[:, cons_trunc:r_augm] @ W_hat

    # Step 9
    Vout = np.zeros((nv, r))
    Vout[:, :cons_trunc] = V[:, :cons_trunc]
    Vout[:, cons_trunc:] = W_np1

    # Step 10
    X_hat = np.zeros((nx, r))
    X_hat[:, :cons_trunc] = X_cons
    X_hat[:, cons_trunc:] = X_rem

    Qx, Rx = qr(X_hat, mode='reduced')
    Xout = Qx / np.sqrt(hx)
    R = Rx * np.sqrt(hx)

    # Step 11
    S_tmp = np.zeros((r, r))
    S_tmp[:cons_trunc, :cons_trunc] = S_cons
    S_tmp[cons_trunc:, cons_trunc:] = S_rem
    Sout = R @ S_tmp

    return Xout, Sout, Vout



def time_step_consBUG(X, S, V):
    # K step
    # TODO
    
    # L step
    # TODO

    # setup the augmented basis
    Xa = zeros((nx,2*r))
    Xa[:,0:r] = X
    Xa[:,r:] = K
    Xa, _ = qr(Xa, mode='reduced')
    Xa *= 1.0/sqrt(hx)
    Xa[:,0:r] = X
    
    Va = zeros((nv,2*r))
    Va[:,0:r] = V
    Va[:,r:] = L
    Va, _ = qr(Va, mode='reduced')
    Va *= 1.0/sqrt(hv)
    Va[:,0:r] = V

    Sa = zeros((2*r,2*r))
    Sa[0:r,0:r] = S

    # Compute the augmentation necessary to get conservation with RK4
    c1 = compute_c1(Va)
    c2 = compute_c2(Va)
    
    # stage 1
    FU_s1 = rhs_K_withE(Xa@Sa, Va, c1, c2)
    Shat_s1 = hx*Xa.transpose()@FU_s1

    # stage 2
    FU_s2 = rhs_K_withE(Xa@(Sa + 0.5*deltat*Shat_s1), Va, c1, c2)
    Shat_s2 = hx*Xa.transpose()@FU_s2
    
    # stage 3
    FU_s3 = rhs_K_withE(Xa@(Sa + 0.5*deltat*Shat_s2), Va, c1, c2)
    Shat_s3 = hx*Xa.transpose()@FU_s3
    
    # stage 4
    FU_s4 = rhs_K_withE(Xa@(Sa + deltat*Shat_s3), Va, c1, c2)

    FU_augm = (FU_s1 + 2.0*FU_s2 + 2.0*FU_s3 + FU_s4)[:,0]

    # augmented with the result of RK4
    # TODO: add FU_augm to Xa to obtain a new basis Xb
    Xbar, _ = qr(Xb, mode='reduced')
    Xbar *= 1.0/sqrt(hx)
    Xbar[0:nx,0:r] = X
    
    # S step
    Sbar_s1 = Xbar.transpose()@FU_s1*hx
    Sbar_s2 = Xbar.transpose()@FU_s2*hx
    Sbar_s3 = Xbar.transpose()@FU_s3*hx
    Sbar_s4 = Xbar.transpose()@FU_s4*hx

    Sbar0 = zeros((2*r+1, 2*r))
    Sbar0[0:r,0:r] = S
    Sout = Sbar0 + deltat/6.0*(Sbar_s1 + 2.0*Sbar_s2 + 2.0*Sbar_s3 + Sbar_s4)
    
    return conservative_truncate(Xbar, Sout, Va)




def plot_all(X, S, V, ee=[], mass_err=[]):
    figure(figsize=(15,5))
    subplot(1,2,1)
    imshow(X.dot(S.dot(V.transpose())).transpose(), extent=[xs[0], xs[-1], vs[0], vs[-1]])
    colorbar()
    xlabel('x')
    ylabel('v')

    subplot(1,2,2)
    plot(xs, compute_E(X.dot(S), V))
    xlabel('x')
    ylabel('E')
    xlim([0,L])
    
    savefig('final.pdf', bbox_inches='tight', pad_inches=0)

    if len(ee) != 0:
        ts = linspace(0,t_final, len(ee))
        figure(figsize=(10,5))
        semilogy(ts, ee)
        xlim([0,t_final])
        savefig('evolution-ee.pdf', bbox_inches='tight', pad_inches=0)
        
        figure(figsize=(10,5))
        semilogy(ts, mass_err)
        xlim([0,t_final])
        savefig('mass_error.pdf', bbox_inches='tight', pad_inches=0)


#
# set the initial value
# (we need to supply the algorithm with an orthonormalized set of basis functions)
#
X = identity(nx)[:,0:r]
V = zeros((nv, r))
S = zeros((r, r))

# linear Landau damping
if problem == 'll':
    V[:,0] = 1.0  # basis function needed for conservation
    X[:,0] = 1.0 + 1e-2*cos(0.5*xs)
    V[:,1] = exp(-0.5*vs**2)/sqrt(2*pi)
    S[0,1] = 1.0
else:
    V[:,0] = 1.0  # basis function needed for conservation
    X[:,0] = 1.0+0.001*cos(0.2*xs)
    V[:,1] = 0.5*(exp(-0.5*(vs-2.4)**2) + exp(-0.5*(vs+2.4)**2))/sqrt(2*pi)
    S[0,1] = 1.0


X, S1 = qr(X, mode='reduced')
V, S2 = qr(V, mode='reduced')
X *= 1.0/sqrt(hx)
V *= 1.0/sqrt(hv)
S = S1.dot(S.dot(S2.transpose()))*sqrt(hx)*sqrt(hv)


#
# run the simulation
#
fs = open('evolution.data', 'w')
fs.write('# t electric_energy\n')

t = 0.0
num_steps = int(ceil(t_final/deltat))
ees = []; mass_errs = [];
for i in range(num_steps):
    if t_final - t < deltat:
        deltat = t_final - t

    E = compute_E(X.dot(S), V)
    ee = electric_energy(E)
    ees.append(ee)

    mass = compute_mass(X, S, V)
    if i==0:
        mass0 = mass
    mass_errs.append(abs(mass-mass0)/mass0)

    fs.write('{0: <30} {1: <30} {1: <30}\n'.format(t, ee, mass_errs))
    fs.flush()

    X, S, V = time_step_augBUG(X, S, V)
    #X, S, V = time_step_consBUG(X, S, V)

    t += deltat
    print('\r', end='')
    print('t={}'.format(t), end='')

fs.close()

plot_all(X, S, V, ees, mass_errs)
