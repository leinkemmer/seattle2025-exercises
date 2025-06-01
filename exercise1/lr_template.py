from pylab import *
from numpy.linalg import solve
from numpy.linalg import qr

nx = 128
nv = 128
r = 5

t_final = 40.0
deltat  = 0.01

L = 4*pi # linear Landau damping


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
    return hv*V.transpose().dot(A_v.dot(V))

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
    return -A_cd_x.dot(K.dot(c1.transpose())) + diag(E).dot(K.dot(c2.transpose()))

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

def time_step(X, S, V):
    # K step
    c1 = compute_c1(V)
    c2 = compute_c2(V)
    K = X.dot(S)
    E = compute_E(K, V)
    K = rk4(deltat, K, lambda K: rhs_K(K, E, c1, c2))
    
    X, S = qr(K, mode='reduced')
    X *= 1.0/sqrt(hx)
    S *= sqrt(hx)

    # S step
    # TODO

    # L step
    # TODO (output is V and Str)

    return X, Str.transpose(), V


def plot_all(X, S, V, ee=[]):
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


#
# set the initial value
# (we need to supply the algorithm with an orthonormalized set of basis functions)
#
X = identity(nx)[:,0:r]
V = zeros((nv, r))
S = zeros((r, r))

# linear Landau damping
S[0,0] = 1.0
X[:,0] = 1.0 + 1e-2*cos(0.5*xs)
V[:,0] = exp(-0.5*vs**2)/sqrt(2*pi)


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
ees = []
for i in range(num_steps):
    if t_final - t < deltat:
        deltat = t_final - t

    E = compute_E(X.dot(S), V)
    ee = electric_energy(E)
    ees.append(ees)

    fs.write('{0: <30} {1: <30}\n'.format(t, ee))
    fs.flush()

    X, S, V = time_step(X, S, V)

    t += deltat
    print('\r', end='')
    print('t={}'.format(t), end='')

fs.close()

plot_all(X, S, V, ees)
