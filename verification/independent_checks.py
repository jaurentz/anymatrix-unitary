#!/usr/bin/env python3
"""Independent NumPy/SciPy verification of the mathematics in MATLAB prototype.

This does NOT execute MATLAB or verify Anymatrix integration. It reconstructs
all algebraic objects independently for mathematical cross-checks.
"""
import json
import numpy as np
from scipy.special import roots_jacobi
from scipy.optimize import linear_sum_assignment


def alpha_formula(k,a,b):
    return (b-a)/(k+a+b+2) if k % 2 == 0 else -(a+b+1)/(k+a+b+2)


def independent_jacobi_recurrence(a,b,m):
    # Compute Jacobi recurrence using Gauss-Jacobi quadrature and a 3-term
    # orthonormal Stieltjes procedure (not the closed formulas in MATLAB).
    x,w = roots_jacobi(max(80,2*m+5),a,b)
    x=2*x;w=w/w.sum()
    pprev=np.zeros_like(x)
    p=np.ones_like(x)
    a_prev=0
    d=[];off=[]
    for k in range(m):
        dk=np.dot(w,x*p*p)
        q=(x-dk)*p-a_prev*pprev
        ak=np.sqrt(np.dot(w,q*q))
        d.append(dk)
        if k<m-1: off.append(ak)
        pprev,p=p,q/ak
        a_prev=ak
    return np.array(d),np.array(off)


def inverse_geronimus(d,off,count):
    alpha=np.empty(count)
    u=2-d[0];v=2+d[0]
    for j in range(count//2+1):
        if j>0:
            u=2-d[j]-off[j-1]**2/u
            v=2+d[j]-off[j-1]**2/v
            if 2*j<=count:alpha[2*j-1]=1-(u+v)/2
        if 2*j<count:alpha[2*j]=(v-u)/(v+u)
    return alpha


def core(U,k,a):
    # k zero indexed; apply core on right.
    rho=np.sqrt((1-abs(a))*(1+abs(a)))
    U[:,k:k+2] = U[:,k:k+2]@np.array([[np.conj(a),rho],[rho,-a]])
    return U


def matrix(n,a,b,beta,kind,twists=None):
    alpha=[alpha_formula(j,a,b) for j in range(n-1)]
    D=np.eye(n,dtype=complex);D[-1,-1]=np.conj(beta)
    if kind=='cmv':
        L=np.eye(n,dtype=complex);M=np.eye(n,dtype=complex)
        for j in range(0,n-1,2):core(L,j,alpha[j])
        for j in range(1,n-1,2):core(M,j,alpha[j])
        if n%2: L[-1,-1]=np.conj(beta)
        else: M[-1,-1]=np.conj(beta)
        return L@M
    if kind=='hessenberg':
        U=np.eye(n,dtype=complex)
        for j in range(n-1):core(U,j,alpha[j])
        U[:,-1]*=np.conj(beta)
        return U
    if kind=='inverse_hessenberg':
        U=D.copy()
        for j in range(n-2,-1,-1):core(U,j,alpha[j])
        return U
    if kind=='twisted':
        if twists is None:twists=[j%2==0 for j in range(n-2)]
        order=[0] if n>1 else []
        for j in range(1,n-1):
            if twists[j-1]:order.insert(0,j)
            else:order.append(j)
        U=np.eye(n,dtype=complex)
        for j in order:core(U,j,alpha[j])
        U[:,-1]*=np.conj(beta)
        return U
    raise ValueError(kind)


def test():
    max_unitarity=0.0;max_eig=0.0;max_moment=0.0;max_alpha=0.0
    n_tests=0
    for a,b in [(0,0),(-.5,-.5),(.5,.5),(.3,1.6),(-.6,.4),(2.4,-.7)]:
        # Independently reconstruct alpha via Gauss-Jacobi quadrature.
        d,off=independent_jacobi_recurrence(a,b,13)
        al=inverse_geronimus(d,off,24)
        err=np.max(np.abs(al-[alpha_formula(k,a,b) for k in range(24)]))
        max_alpha=max(max_alpha,err)
        assert err<1e-8, (a,b,err)
        xx,ww=roots_jacobi(100,a,b);ww=ww/ww.sum()
        for n in [1,2,3,5,8,12]:
            for beta in [1,np.exp(.37j)]:
                H=matrix(n,a,b,beta,'hessenberg')
                heig=np.linalg.eigvals(H)
                for kind in ['hessenberg','inverse_hessenberg','cmv','twisted']:
                    U=matrix(n,a,b,beta,kind)
                    res=np.linalg.norm(U.conj().T@U-np.eye(n),'fro')
                    max_unitarity=max(max_unitarity,res)
                    assert res<1e-10,(a,b,n,kind,res)
                    eig=np.linalg.eigvals(U)
                    ri,ci=linear_sum_assignment(np.abs(heig[:,None]-eig[None,:]))
                    err=np.max(np.abs(heig[ri]-eig[ci]))
                    max_eig=max(max_eig,err)
                    assert err<1e-9,(a,b,n,kind,err)
                    kval=min(n-1,6)
                    V=np.linalg.eig(U)[1]
                    weights=np.abs(V[0,:])**2
                    for k in range(1,kval+1):
                        mtrue=np.dot(ww,np.cos(k*np.arccos(xx)))
                        m=np.dot(weights,eig**k)
                        err=abs(m-mtrue)
                        max_moment=max(max_moment,err)
                        assert err<1e-9,(a,b,n,kind,k,err)
                    n_tests+=1
    report={'method':'Independent SciPy/NumPy reconstruction, NOT MATLAB execution',
        'matrix_instances':n_tests, 'max_frobenius_unitarity_residual':max_unitarity,
        'max_eigenvalue_matching_error':max_eig,
        'max_quadrature_moment_error':max_moment,
        'max_inverse_geronimus_coefficient_error':max_alpha}
    return report


if __name__=='__main__':
    result=test()
    with open(__file__.replace('.py','_results.json'),'w') as f:
        json.dump(result,f,indent=2)
    print(json.dumps(result,indent=2))
