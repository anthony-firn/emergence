# Inviscid Rayleigh eigenproblem for the diffuse slab profile used by capillary_probe rt:
#   d/dz(rho w') - k^2 rho w = -(g k^2/n^2) rho' w,  periodic z.
import numpy as np, scipy.linalg as sl
def growth(N, ratio=857.0, k=2.0, g=1.0, sigma=0.02, M=1024):
    eta = 1.2*2*np.pi/N; w = np.sqrt(2)*eta
    z = np.arange(M)*2*np.pi/M; h = z[1]-z[0]
    zlo, zhi = 0.5*np.pi, 1.5*np.pi
    def rho_of(z):
        # periodic images of the slab
        phi = np.tanh((z-zlo)/w) - np.tanh((z-zhi)/w) - 1.0
        c = np.clip(phi, -1, 1)
        r1, r2 = 1.0, 1.0/ratio
        return 0.5*(r1+r2) + 0.5*(r1-r2)*c
    rho = rho_of(z); rhoh = rho_of(z+0.5*h)  # rho at i+1/2
    # (-L) w : -(rho_{i+1/2}(w_{i+1}-w_i) - rho_{i-1/2}(w_i-w_{i-1}))/h^2 + k^2 rho_i w_i
    A = np.zeros((M,M)); 
    for i in range(M):
        ip, im = (i+1)%M, (i-1)%M
        rp, rm = rhoh[i], rhoh[im]
        A[i,i] = (rp+rm)/h**2 + k*k*rho[i]
        A[i,ip] -= rp/h**2; A[i,im] -= rm/h**2
    drho = (rhoh - np.roll(rhoh,1))/h
    B = np.diag(drho)
    mu = sl.eigh(B, A, eigvals_only=True)   # B w = mu A w
    n2 = g*k*k*mu.max() - sigma*k**3/(1+1/ratio)
    return np.sqrt(n2)
sharp = np.sqrt((1-1/857)/(1+1/857)*2 - 0.02*8/(1+1/857))
for N in (32, 64):
    n = growth(N); print(f"N={N}: diffuse-profile theory n={n:.4f} ({100*(n/sharp-1):+.2f}% vs sharp {sharp:.4f})")
