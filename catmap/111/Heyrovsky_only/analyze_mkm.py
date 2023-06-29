#!/usr/bin/env python
import sys, os
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib import cm
import numpy as np
import pickle as pkl
from scipy.optimize import curve_fit
from general_tools import lin_fun

plt.rcParams["figure.figsize"] = (10,5)
#plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.size"] = 16
plt.rc('axes', labelsize=28)    # fontsize of the x and y labels
plt.rcParams['xtick.labelsize'] = 18
plt.rcParams['ytick.labelsize'] = 18
plt.rcParams['figure.figsize'] = (12,7)
markersize=10

home=os.getcwd()
path_info=home.split('/')
add_total_rate_as_panel=False #Not working yet

nsteps=2

fig,ax=plt.subplots(nsteps)

def main():
    from catmap import ReactionModel
    model = ReactionModel(setup_file = 'reactions.mkm')
    model.output_variables+=['production_rate', 'free_energy','coverage']
    model.run()
    data,phs,pots=extract_data_from_mkm(model)
    plot_heatplot(data,phs,pots,nsteps)
    ax[1].annotate('*H',(-1.6,13),ha='left',va='top',fontsize=30)
    ax[0].annotate('*CO',(-1.6,13),ha='left',va='top',fontsize=30)
    plt.show()

def plot_2D_plot(data,phs,pots,nsteps,phout,title):
    plt.rcParams['figure.figsize'] = (10,5)
    plt.rcParams["font.family"] = "Times New Roman"
    plt.rcParams["font.size"] = 16
    plt.rc('axes', labelsize=28)    # fontsize of the x and y labels
    plt.rcParams['xtick.labelsize'] = 18
    plt.rcParams['ytick.labelsize'] = 18
    markersize=10
    fig2d,ax2d = plt.subplots(1,2,figsize=[15,5],sharex=True)


    colors=['b','peru','r','g']
    lines=['-','--']
    labels=['A','B']
    X=np.array(sorted(pots))
    Y=np.array(sorted(phs))

    iphs = np.where([i in phout for i in Y])[0]
    rate=np.ones((len(X)))*0.5
    dr=np.ones((len(X)-1))
    sel=np.ones((len(X)))*0.5
    for iph,ph in enumerate(phout):
      ax2d[0].plot(np.nan,np.nan,'k',label=f'pH{ph}',linestyle=lines[iph])
    #  ax2d[1].plot(np.nan,np.nan,'k',label=f'pH{ph}',linestyle=lines[iph])
      for istep in range(nsteps):
             for ipot,pot in enumerate(X):
                rate[ipot]=data[pot][ph][istep]
                if ipot:
                    dr[ipot-1]=(np.log10(rate[ipot])-np.log10(rate[ipot-1]))/(X[ipot]-X[ipot-1])*0.059

             ax2d[0].plot(X,rate,color=colors[istep],linestyle=lines[iph],linewidth=2)
             ax2d[1].plot(X[1:],-dr,color=colors[istep],linestyle=lines[iph],linewidth=2)
             ax2d[0].annotate(labels[istep],(min(X)+0.1,max(rate)),color=colors[istep],fontsize=30).draggable()
    ax2d[0].set_ylabel('TOF / s$^{-1}$')
    ax2d[1].set_ylabel('Transfer coefficient')
    ax2d[0].set_xlabel('U$_{\mathrm{SHE}}$ / V')
    ax2d[1].set_xlabel('U$_{\mathrm{SHE}}$ / V')
    ax2d[0].set_yscale('log')
#    ax2d[0].set_ylim([1e-4,10000])
    ax2d[1].set_ylim([0.0,0.4])
    ax2d[0].set_xlim([-1.6,-1.0])
#    ax2d[0].legend()
#    ax2d[1].legend()
    ax2d[0].set_title(title)#' '*120+r'$\alpha_{\mathrm{RLS}}$=0.2, $\alpha_{\mathrm{A}}$=0.5, $\alpha_{\mathrm{B}}$=0.5, $\Delta$G$_\mathrm{A}$=1.05eV, $\Delta$G$_\mathrm{B}$=1.0eV')
    fig2d.tight_layout()

    plt.savefig('figure.pdf')
    plt.show()



def plot_heatplot(data,phs,pots,nsteps):
    X=np.array(sorted(pots))
    Y=np.array(sorted(phs))

    for col in range(1):
     for istep in range(nsteps):
        R,S = get_rate_and_selectivity(col,istep,data,nsteps,X,Y)
        plot_it(R,S,fig,ax,col,istep,X,Y)
        if 0:
         if istep == nsteps-1:
            for thisax in ax[istep]:
                thisax.set_xlabel('U$_{SHE}$ [V]')
         else:
            for thisax in ax[istep]:
                thisax.set_xticks([])

        ax[istep].set_ylabel('pH')
        #ax[istep][1].set_yticks([])


def get_rate_and_selectivity(col,istep,data,steps,X,Y):

    Selectivity=np.ones((len(X),len(Y)))*0.5
    rate=np.ones((len(X),len(Y)))*0.5

    for ix,x in enumerate(X):
       for iy,y in enumerate(Y):
        try:
            if col == 1:
                Selectivity[ix][iy]=data[x][y][istep]/np.sum(data[x][y][:nsteps])
            else:
                rate[ix][iy]=data[x][y][istep]#/np.sum(data[x][y][:nsteps])
        except:
            Selectivity[ix][iy]=np.nan#data[x][y][istep]#/np.sum(data[x][y][:3])
            rate[ix][iy]=1e-20#np.nan#data[x][y][istep]#/np.sum(data[x][y][:nsteps])
    return rate, Selectivity

def plot_it(R,S,fig,ax,col,istep,X,Y):
        if col == 0:
         b = ax[istep].imshow(R.T,
                interpolation='bicubic',
                cmap=cm.jet,
                   origin='lower', extent=[X.min(), X.max(), Y.min(), Y.max()],norm=LogNorm(),#,
                    vmin=1e-3,
                    vmax=1e0,#)
                    aspect='auto')#, vmin=-abs(alldata[:,2]).max())

        else:
         a = ax[istep][1].imshow(S.T,
                interpolation='bicubic',
                cmap=cm.RdYlGn,
                   origin='lower', extent=[X.min(), X.max(), Y.min(), Y.max()],#norm=LogNorm(),#,
                    vmin=0,
                    vmax=1,
                    aspect='auto')#, vmin=-abs(alldata[:,2]).max())
        if istep == 1:
            fig.colorbar(b,ax=ax,shrink=1,label='Coverage')#,orientation='horizontal') #,location='top',orientation='horizontal')
#        plt.colorbar(b)


def extract_data_from_mkm(model):
    data={}
    pots,phs=[],[]
    for dat in model.coverage_map:
        pot,ph=np.around(dat[0][0],3),np.around(dat[0][1],3)
        if pot not in data:
            data[pot] = {}
        data[pot][ph] = dat[1]
        if pot not in pots:
            pots.append(pot)
        if ph not in phs:
            phs.append(ph)
    return data,phs,pots

def read_data(infile='mkm.pkl'):
    data_in = pkl.load(open(infile,'rb'),encoding='latin1')
    data={}
    pots,phs=[],[]
    for dat in data_in['coverage_map']:
        pot,ph=np.around(dat[0][0],3),np.around(dat[0][1],3)
        if pot not in data:
            data[pot] = {}
        data[pot][ph] = dat[1]
        if pot not in pots:
            pots.append(pot)
        if ph not in phs:
            phs.append(ph)
    return data,phs,pots
if __name__ == "__main__":
    main()
