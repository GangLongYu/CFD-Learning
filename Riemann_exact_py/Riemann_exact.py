import numpy as np
import matplotlib.pyplot as plt

def f_Riemann(p, pi, ri, gamma):
    ci = np.sqrt(gamma * pi / ri) # known sound speed
    if p < 0:
        print('Error, pressure < 0 !!!\n')
    elif p < 1e-4: # p = 0
        fp = -2 * ci / (gamma - 1)
    elif p > pi: # rarefraction wave
        fp = (p-pi)/(ri*ci)/np.sqrt((gamma+1)/(2*gamma)*(p/pi)+(gamma-1)/(2*gamma))
    else:
        fp = 2*ci/(gamma-1) * ((p/pi)**((gamma-1)/2/gamma)-1)
    return fp


def Fp_Riemann(p, p_L, r_L, p_R, r_R, gamma):
    Fp = f_Riemann(p, p_L, r_L, gamma) + f_Riemann(p, p_R, r_R, gamma)
    return Fp


# pressure in the central region by using Newton iteration
def pstar_Newton(p_L, u_L, r_L, p_R, u_R, r_R, gamma): 
    error_level = 1e-15
    error = 1e4
    delta_p = 1e-4 # get derivative dF/dp
    du = u_L - u_R

    pstar = np.min([p_L, p_R]) # initial value of pstar
    while error > error_level:
        dF_dp = (Fp_Riemann(pstar+delta_p,p_L,r_L,p_R,r_R,gamma)-Fp_Riemann(pstar,p_L,r_L,p_R,r_R,gamma))/delta_p # F'(pstar)
        pstar_new = pstar - (Fp_Riemann(pstar,p_L,r_L,p_R,r_R,gamma)-du)/dF_dp
        if pstar_new < 0:
            pstar_new = 0 # p<0 non-physical
        error = abs(pstar_new - pstar)
        pstar = pstar_new
        print('"pstar = {}, ", "abs(pstar_new - pstar) = {}"'.format(pstar, error))
    return pstar


def Riemann_exact_solver(p_L, u_L, r_L, p_R, u_R, r_R, gamma, nx, t):
    r = np.zeros((nx,))
    u = r.copy()
    p = r.copy()
    x = np.array([-1+2*i/(nx-1) for i in range(nx)]) # x in [-1, 1]

    c_L = np.sqrt(gamma*p_L/r_L) # sound speed
    c_R = np.sqrt(gamma*p_R/r_R)
    pstar = pstar_Newton(p_L, u_L, r_L, p_R, u_R, r_R, gamma)
    ustar = (u_L-f_Riemann(pstar,p_L,r_L,gamma) + u_R+f_Riemann(pstar,p_R,r_R,gamma))/2
    den_L = r_L*c_L*np.sqrt((gamma+1)/(2*gamma)*(pstar/p_L)+(gamma-1)/(2*gamma)) # f(p)激波表达式的分母
    den_R = r_R*c_R*np.sqrt((gamma+1)/(2*gamma)*(pstar/p_R)+(gamma-1)/(2*gamma)) 

    if pstar > p_L: # shock wave (left)
        z_L = u_L - den_L / r_L
        rstar_L = r_L*den_L/(den_L+r_L*(ustar-u_L))

        index = x < z_L * t                                      # 左激波左侧
        p[index] = p_L; u[index] = u_L; r[index] = r_L; 
        a = x >= z_L * t; b = x < ustar * t; index = a & b       # 左激波到接触间断
        p[index] = pstar; u[index] = ustar; r[index] = rstar_L; 
        if pstar > p_R: # shock wave (right)
            z_R = u_R + den_R / r_R
            rstar_R = r_R*den_R/(den_R-r_R*(ustar-u_R))

            a = x >= ustar * t; b = x < z_R * t; index = a & b   # 接触间断到右激波左侧
            p[index] = pstar; u[index] = ustar; r[index] = rstar_R; 
            index = x >= z_R * t                                 # 右激波右侧
            p[index] = p_R; u[index] = u_R; r[index] = r_R; 
        else:           # rarefraction wave (right)
            rstar_R = r_R * (pstar/p_R)**(1/gamma)
            cstar_R = np.sqrt(gamma*pstar/rstar_R)
            z_R_head = u_R + c_R
            z_R_tail = ustar + cstar_R

            a = x >= ustar * t; b = x < z_R_tail * t; index = a & b    # 接触间断到右稀疏波左侧
            p[index] = pstar; u[index] = ustar; r[index] = rstar_R; 
            a = x >= z_R_tail * t; b = x < z_R_head * t; index = a & b # 右稀疏波内部
            ca = (gamma-1)/(gamma+1)*(x[index]/t-u_R) + 2/(gamma+1)*c_R
            u[index] = x[index]/t - ca; p[index] = p_R*(ca/c_R)**(2*gamma/(gamma-1)); r[index] = gamma*p[index]/(ca*ca)
            index = x >= z_R_head * t                     # 右稀疏波右侧
            p[index] = p_R; u[index] = u_R; r[index] = r_R; 
    else:           # rarefraction wave (left)
        rstar_L = r_L * (pstar/p_L)**(1/gamma)
        z_L_head = u_L - c_L
        if rstar_L < 1e-4: # 两稀疏波中间为真空，rho=0，声速无定义
            z_L_tail = u_L + 2*c_L/(gamma-1)
        else:
            cstar_L = np.sqrt(gamma*pstar/rstar_L)
            z_L_tail = ustar - cstar_L

        index = x < z_L_head * t                                   # 左稀疏波左侧
        p[index] = p_L; u[index] = u_L; r[index] = r_L; 
        a = x >= z_L_head * t; b = x < z_L_tail * t; index = a & b # 左稀疏波内部
        ca = (gamma-1)/(gamma+1)*(u_L-x[index]/t)+2*c_L/(gamma+1)
        u[index] = x[index]/t + ca; p[index] = p_L*(ca/c_L)**(2*gamma/(gamma-1)); r[index] = gamma*p[index]/(ca*ca)
        if pstar > p_R: # shock wave (right)
            z_R = u_R + den_R / r_R
            rstar_R = r_R*den_R/(den_R-r_R*(ustar-u_R))

            a = x >= z_L_tail * t; b =  x < ustar * t; index = a & b
            p[index] = pstar; u[index] = ustar; r[index] = rstar_L; 
            a = x >= ustar * t; b = x < z_R * t; index = a & b   # 接触间断到右激波左侧
            p[index] = pstar; u[index] = ustar; r[index] = rstar_R; 
            index = x >= z_R * t                                 # 右激波右侧
            p[index] = p_R; u[index] = u_R; r[index] = r_R; 
        else:           # rarefraction wave (right)
            rstar_R = r_R * (pstar/p_R)**(1/gamma)
            z_R_head = u_R + c_R
            if rstar_R < 1e-4: # 中间真空
                z_R_tail = u_R - 2*c_R/(gamma-1)
                a = x >= z_L_tail * t; b =  x < z_R_tail * t; index = a & b
                p[index] = 0; u[index] = 0; r[index] = 0; 
            else:              # 中间非真空，存在接触间断
                cstar_R = np.sqrt(gamma*pstar/rstar_R)
                z_R_tail = ustar + cstar_R
                a = x >= z_L_tail * t; b =  x < ustar * t; index = a & b
                p[index] = pstar; u[index] = ustar; r[index] = rstar_L; 
                a = x >= ustar * t; b = x < z_R_tail * t; index = a & b    # 接触间断到右稀疏波左侧
                p[index] = pstar; u[index] = ustar; r[index] = rstar_R; 

            a = x >= z_R_tail * t; b = x < z_R_head * t; index = a & b # 右稀疏波内部
            ca = (gamma-1)/(gamma+1)*(x[index]/t-u_R) + 2/(gamma+1)*c_R
            u[index] = x[index]/t - ca; p[index] = p_R*(ca/c_R)**(2*gamma/(gamma-1)); r[index] = gamma*p[index]/(ca*ca)
            index = x >= z_R_head * t                                  # 右稀疏波右侧
            p[index] = p_R; u[index] = u_R; r[index] = r_R; 

    return r, u, p, x


if __name__ == '__main__':
    #p_L, u_L, r_L, p_R, u_R, r_R = [1, 0, 1, 0.1, 0, 0.125] # case 1
    #p_L, u_L, r_L, p_R, u_R, r_R = [1, 3, 1, 0.2, -1, 0.1]  # case 2
    #p_L, u_L, r_L, p_R, u_R, r_R = [1, -2, 1, 0.2, 2, 0.4]  # case 3
    #p_L, u_L, r_L, p_R, u_R, r_R = [0.1, 0, 0.125, 1, 0, 1]  # case 4
    p_L, u_L, r_L, p_R, u_R, r_R = [0.2, -8, 1, 0.3, 6, 0.1]  # case 5
    
    gamma = 1.4
    
    # check FP_Riemann
    '''
    p = np.linspace(0, 2, 50)
    Fp = []
    for i in p:
        temp = Fp_Riemann(i, p_L, r_L, p_R, r_R, gamma)
        Fp.append(temp)
    plt.plot(p, Fp, c='r', marker='>')
    plt.xlabel('p')
    plt.ylabel('F(p)')
    '''
    
    # main function
    print('===========Exact Riemann Solve, Python Version===========\n')
    nx = int(input('please input nx (grid number) for plot, e.g. nx = 201\n'))
    t = float(input('please input time t, e.g. t = 0.14\n'))
    
    r, u, p, x = Riemann_exact_solver(p_L, u_L, r_L, p_R, u_R, r_R, gamma, nx, t)
    plt.figure(1); plt.plot(x, r); plt.xlabel('x'); plt.ylabel('density'); plt.show()
    plt.figure(2); plt.plot(x, u); plt.xlabel('x'); plt.ylabel('velocity'); plt.show()
    plt.figure(3); plt.plot(x, p); plt.xlabel('x'); plt.ylabel('pressure'); plt.show()

    with open('Riemann.dat', 'w') as f:
        print('variables=x,r,u,p')
        for i in range(nx):
            print('%20.12e ' % x[i], '%20.12e' % r[i], '%20.12e' % u[i], '%20.12e' % p[i], file=f)

    print("OK, the flow data are writen to 'Riemann.dat' as a tecplot file\n")