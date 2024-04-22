import numpy as np
import math
import matplotlib.pyplot as plt

def m_local(l,m):
    m_ = np.array([[156, 22 * l, 54, -13 * l], [22 * l, 4 * l * l, 13 * l, -3 * l * l], [54, 13 * l, 156, -22 * l],
                   [-13 * l, - 3 * l * l, - 22 * l, 4 * l * l]])
    m_ = (m * l / 420) * m_
    return (m_)

def k_local(EI,m,l,T,order,i):
    k = np.array([[12, 6 * l, - 12, 6 * l], [6 * l ,4 * l * l, -6 * l, 2 * l * l], [-12, - 6 * l ,12, - 6 * l],
                  [6 * l ,2 * l * l, - 6 * l ,4 * l * l]])
    k = (EI / (l ** 3)) * k
    t=T-(i*l*m)
    kg = np.array(
        [[6 / 5, l / 10, -6 / 5, l / 10], [l / 10, 2 * l * l / 15, -l / 10, -l * l / 30], [-6 / 5, -l / 10, 6 / 5, -l / 10],
         [l / 10, -l * l / 30, -l / 10, 2 * l * l / 15]])
    kg = (t / l) * kg
    k_local = k + kg
    return(k_local)

def local_force(q,l):
    f= np.array([6, l, 6, -l]).reshape(-1, 1)
    f = (q * l / 12) * f
    return(f)

def newmark(M, K, C, x0, v0, F0, dt):
     a = 0.5
     b = 0.25
     a0 = np.zeros((1, 3))
     a0 = ((F0) - np.dot(C, v0) - np.dot(K, x0))
     a0 = np.dot(np.linalg.inv(M), a0)
     # a0=((F0)-(C*v0)-(K*x0))/M
     k_cap = K + ((a / (b * dt)) * C) + ((1 / (b * dt * dt)) * M)
     A = ((1 / (b * dt)) * M) + ((a / b) * C)
     B = ((1 / (2 * b)) * M) + (((a / (2 * b)) - 1) * C)
     p_cap = F0 + np.dot(A, v0) + np.dot(B, a0)
     dx = np.dot(np.linalg.inv(k_cap), p_cap)
     dv = ((a / (b * dt)) * dx) - ((a / b) * v0) + (((a / (2 * b)) - 1) * dt * a0)
     da = ((1 / (b * dt * dt)) * dx) - ((1 / (b * dt)) * v0) - ((1 / (2 * b)) * a0)
     return (x0 + dx, v0 + dv, a0 + da)

if __name__ == '__main__':
    element_no=100
    node_no=element_no+1
    d=node_no
    order=2*node_no
    L=2262
    l=L/element_no
    dt=0.001
    Do = 533.4 / 1000
    t = 15.875 / 1000
    Db = 1.3716
    D_r = 7850
    E = 210*(10**9)
    D_s = 1025
    D_rf = 1200
    Kbt = 8800
    Kbb = 127400
## Envirnmental data
    wh = 6
    wl = 10.2
    wp = 11.2
    Vw = 15.43
    Vd = 0.93
    Cd = 1.2
    Cm = 2.0
    g = 9.81
    L0 = 1.56 * wp**2
    k = 2 *np.pi / L0
    Di = Do - 2 * t
    Ae = np.pi *( Do**2) / 4
    Ai = np.pi * (Di**2) / 4
    An = np.pi*((Do**2)-(Di**2))/4
    I = np.pi*((Do**4) -(Di**4)) / 64
    EI = E * I
    FD = 100
    mf = D_rf * Ai
    mr = D_r * An
    m = (mr + mf)/10
    T=(D_r*An+Ai*D_rf-Ae*D_s)*g*2262
    Vcm = 0.03 * Vw
    # Initialize arrays
    Uw = np.zeros((d, 1))
    Uc = np.zeros((d, 1))
    q_f = np.zeros((d, 1))
    F0 = np.zeros((2 * d, 1))

    # Load calculation loop
    i = 0
    for y in range(1, d + 1):
         Uw[i, 0] = np.pi * wh * np.cosh(k * (d - y)) / (wp * np.sinh(k * d))  # Wave velocity in m/s
         if y >= d - FD:
             Uc[i, 0] = Vcm * ((y - d + FD) / FD) + Vd * (y / d) ** (1 / 7)
         else:
             Uc[i, 0] = Vd * (y / d) ** (1 / 7)
         q_f[i, 0] = 0.5 * Cd * D_s * Do * np.abs(Uw[i, 0]) * Uw[i, 0] + 0.5 * Cd * Do * D_s * Uc[i, 0] * Uc[i, 0]
         i += 1

    M=np.zeros((order,order))

    K = np.zeros((order, order))
    for i in range(1,element_no+1):
       s = slice(2 * i - 2, 2 * i + 2)
       K[s,s]=K[s,s]+k_local(EI,m,l,T,order,i)
       M[s,s]=M[s,s]+m_local(l,m)
       q = q_f[i - 1]
       F0[2*i-2:2*i+2] = F0[2*i-2:2*i+2]+local_force(q,l)


    #Natural frequency
    eigenvalues_squared, eigenvectors = np.linalg.eig(np.dot(np.linalg.inv(M),K))
    natural_frequencies = np.sqrt(eigenvalues_squared)
    #print(natural_frequencies)
    natural_frequencies = sorted(natural_frequencies)
    w1=natural_frequencies[0]
    w2=natural_frequencies[1]
    w1=np.pi*2/w1
    w2=np.pi*2/w2
    #w1=2.68
    #w2=5.46

    # damping matrix
    zeta=0.03
    pi=math.pi
    alpha=(2*zeta*w1*w2*2*pi)/(w1+w2)
    beta=zeta/((w1+w2)*pi)
    print(alpha, beta)
    C=(alpha*M)+(beta*K)

    #Newmarks method
    x0= np.zeros((order, 1))
    v0= np.zeros((order, 1))
    x1, v1, a1=newmark(M, K, C, x0, v0, F0, dt)
    x1[0]=0
    x1[-2]=0
    K= np.array_str(K)
    print("K=",K)
    displ=x1[::2]
    #print("dimension1 =", len(x1))
    #print(displ)
    #print("v=", x1)
    #v1 = np.array_str(v1, precision=2, suppress_small=True)
    #print("v1=",v1)
    #a1 = np.array_str(a1, precision=2, suppress_small=True)
    #print("a1=",a1)
    column_vector = np.arange(0,d,1).reshape(-1,1)
    plt.plot(displ,column_vector, label='Data Line')
    plt.show()
