#!/usr/bin/env python3
"""Regenerate the unitary spectral atlas examples.

Documentation helper; requires numpy and matplotlib (not MATLAB).
Reproduces the analytic examples and circular style in the atlas, not the
SVG assets byte-for-byte.

Run: python docs/figures/plot_spectral_atlas.py --out atlas_plots
"""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

R = 0.94
TWO_PI = 2 * np.pi

def wrap(x):
    return (np.asarray(x) + np.pi) % TWO_PI - np.pi

def quadrature(fn, n=2800):
    th = -np.pi + TWO_PI * (np.arange(n) + .5) / n
    w = np.maximum(0., fn(th))
    return th, w / w.sum()

def quantiles(th, weights, n=140):
    order = np.argsort(th)
    th, w = th[order], weights[order]
    return np.interp((np.arange(n)+.5)/n, np.cumsum(w), th)

def smoothing(grid, phases, weights):
    assert abs(np.sum(weights)-1.) < 1e-10
    result = np.empty(len(grid))
    for k,t in enumerate(grid):
        result[k] = np.sum(weights * (1-R*R) /
                           (TWO_PI * (1-2*R*np.cos(t-phases)+R*R)))
    return result

def examples():
    # Exact finite counting measures: DFT, Grover, cyclic shift.
    A = (np.array([0.,np.pi/2,-np.pi,-np.pi/2]),
         np.array([17,15,16,16],dtype=float)/64)
    alpha = 2*np.arcsin(np.sqrt(2/128))
    B = (np.array([-np.pi,0.,alpha,-alpha]),
         np.array([125,1,1,1],dtype=float)/128)
    C = (wrap(TWO_PI*np.arange(160)/160),np.ones(160)/160)
    # Underlying Bernstein-Szego measure, not CMV trace law.
    D = quadrature(lambda t:(1-.6**2)/
                   (TWO_PI*(1+.6**2-2*.6*np.cos(t))))
    # Exact finite cycle propagator and Dirac phase laws.
    E = (wrap(-1.6*np.cos(TWO_PI*np.arange(256)/256)),np.ones(256)/256)
    nu = np.sqrt(1-.55**2)
    omega = np.arccos(nu*np.cos(TWO_PI*np.arange(256)/256))
    F = (np.r_[omega,-omega],np.ones(512)/512)
    # Spherical rotation in scaling delta=c/L, exact multiplicity law.
    L,c=70,1.4
    m=np.arange(-L+1,L)
    G=(m*c/L,(L-np.abs(m))/L**2)
    # Underlying generalized circular Jacobi measure.
    H=quadrature(lambda t:(4*np.sin(t/2)**2)**(-.35),n=3400)
    # Underlying Geronimus e_1 measure with an isolated atom.
    a=.3+.2j
    edge=2*np.arcsin(abs(a))
    b=2*np.angle(1+np.conj(a))
    mass=2*(a.real+abs(a)**2)/abs(1+a)**2
    def continuous(t):
        t=wrap(t)
        ans=np.zeros_like(t)
        ok=np.abs(t)>edge
        ans[ok]=np.sqrt(np.maximum(0,np.sin(t[ok]/2)**2-abs(a)**2))/(
            TWO_PI*abs(1+a)*np.abs(np.sin((t[ok]-b)/2)))
        return ans
    iph,iw=quadrature(continuous,n=5000)
    I=(np.r_[iph,b],np.r_[iw*(1-mass),mass])
    return [
        ('A. Fourier atoms',A,A[0],True,None),
        ('B. Grover atom',B,B[0],True,None),
        ('C. Uniform limit',C,C[0],False,None),
        ('D. Smooth OPUC',D,quantiles(*D),False,'quantiles'),
        ('E. One singular band',E,E[0],False,None),
        ('F. Two singular bands',F,F[0],False,None),
        ('G. Soft band edges',G,G[0],False,None),
        ('H. Jacobi singularity',H,quantiles(*H),False,'quantiles'),
        ('I. Geronimus mixed',I,quantiles(iph,iw,105),True,b),
    ]

def draw(ax,entry,grid):
    title,(phases,weights),hairlines,show_atoms,marker=entry
    assert np.isfinite(weights).all() and (weights>=0).all()
    assert abs(weights.sum()-1)<1e-10
    rho=smoothing(grid,phases,weights)
    ax.set_theta_zero_location('E')
    ax.set_theta_direction(1)
    ax.set_rlim(0,1.7)
    ax.set_axis_off()
    for r in (.40,.65,.85,1.,1.18,1.4,1.6):
        ax.plot(grid,np.full_like(grid,r),color='#d5d8dc',
                lw=.65 if r!=1 else 1.0,zorder=0)
    for ang in np.arange(12)*np.pi/6:
        ax.plot([ang,ang],[0,1.69],color='#e5e7e9',lw=.65,zorder=0)
    for ang in hairlines:
        ax.plot([ang,ang],[.945,1.055],color='#424a51',
                lw=.48,alpha=.7,solid_capstyle='round',zorder=2)
    if show_atoms:
        if marker is not None and not isinstance(marker,str):
            ax.plot([marker,marker],[.84,1.16],color='#13191e',lw=2.8,zorder=3)
        else:
            for ang,w in zip(phases,weights):
                ax.plot([ang,ang],[.88,1.12],color='#15191e',
                        lw=1+3*np.sqrt(w/weights.max()),
                        solid_capstyle='round',zorder=3)
    ax.plot(grid,1+.54*rho/max(rho),color='#161b20',lw=1.55,zorder=4)
    ax.set_title(title,fontsize=11,pad=16)

def save_gallery(out,examples,indices,filename,rows,cols):
    fig=plt.figure(figsize=(4.3*cols,4.4*rows),dpi=180)
    grid=np.linspace(-np.pi,np.pi,600)
    for i,k in enumerate(indices):
        ax=fig.add_subplot(rows,cols,i+1,projection='polar')
        draw(ax,examples[k],grid)
    fig.suptitle('Unitary spectral-distribution atlas (Poisson r = 0.94)',
                 fontsize=17,y=.995)
    fig.tight_layout(rect=[.01,.01,.99,.97])
    for ext in ('png','svg'):
        fig.savefig(out/f'{filename}.{ext}',dpi=220,bbox_inches='tight')
    plt.close(fig)

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,default=Path('atlas_plots'))
    args=parser.parse_args()
    args.out.mkdir(exist_ok=True,parents=True)
    ex=examples()
    save_gallery(args.out,ex,range(9),
                 'spectral_distribution_classification_plate',3,3)
    save_gallery(args.out,ex,[4,5,6,8],
                 'spectral_distribution_edge_gallery',1,4)
    print('Wrote atlas images to',args.out)

if __name__ == '__main__':
    main()
