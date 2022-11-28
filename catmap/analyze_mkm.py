#!/usr/bin/env python
import sys, os
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib import cm
import numpy as np
import pickle as pkl
from scipy.optimize import curve_fit
from general_tools import lin_fun
from catmap import analyze

plt.rcParams["figure.figsize"] = (10,5)
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.size"] = 16
plt.rc('axes', labelsize=28)    # fontsize of the x and y labels
plt.rcParams['xtick.labelsize'] = 18
plt.rcParams['ytick.labelsize'] = 18
plt.rcParams['figure.figsize'] = (12,7)
markersize=10

home=os.getcwd()
path_info=home.split('/')
if len(sys.argv) < 2:
    print('Facet needs to be given in the command line!')
    exit()
facet=sys.argv[1]

cov_and_rate={'coverage':['H_'+facet],'production_rate': ['H2_g']}

def main():
    model=run_catmap(facet)

    dattype='coverage'
    steps=[]
    for i,label in enumerate(model.output_labels[dattype]):
        if label in cov_and_rate[dattype]:
            steps.append(i)
    nsteps=len(steps)

    fig,ax=plt.subplots(nsteps)

    data,phs,pots=extract_data_from_mkm(model)
    plot_heatplot(data,phs,pots,steps,fig,ax)

    if nsteps == 1:
#        ax.annotate('*H',(-1.7,13.5),ha='left',va='top',fontsize=40,fontweight='bold')
#        ax.set_title(f'Cu({facet})',fontsize=40)
        ax.set_xlabel('U vs SHE / V')
        if facet == '100':
            ax.annotate('(a) Cu(100)',(-1.7,13.5),ha='left',va='top',fontsize=40)
        else:
            ax.annotate('(b) Cu(211)',(-1.7,13.5),ha='left',va='top',fontsize=40)
    else:
        ax[1].annotate('*H',(-1.7,13.5),ha='left',va='top',fontsize=40,fontweight='bold')
        ax[0].annotate('*CO',(-1.7,13.5),ha='left',va='top',fontsize=40,fontweight='bold')
        ax[1].set_xlabel('U vs SHE / V')

        ax[0].set_title(f'Cu({facet})',fontsize=40)
    plt.savefig(f'../results/Coverages_{facet}.pdf')
    plt.show()

def run_catmap(facet,runbasedir=home):
    os.chdir(runbasedir+'/'+facet)
    from catmap import ReactionModel
    model = ReactionModel(setup_file = 'reactions.mkm')
    model.output_variables+=['production_rate', 'free_energy','coverage']
    model.run()
    os.chdir(home)
    return model

def plot_heatplot(data,phs,pots,steps,fig,ax,dattype='coverage'):
    X=np.array(sorted(pots))
    Y=np.array(sorted(phs))
    nsteps=len(steps)

    for col in range(1):
     for istep in steps:
        R,S = get_rate_and_selectivity(col,istep,data,nsteps,X,Y)
        plot_it(R,S,fig,ax,col,istep,X,Y,nsteps)
        if 0:
         if istep == nsteps-1:
            for thisax in ax[istep]:
                thisax.set_xlabel('U$_{SHE}$ [V]')
         else:
            for thisax in ax[istep]:
                thisax.set_xticks([])
        if nsteps == 1:
            ax.set_ylabel('pH')
        else:
            ax[istep].set_ylabel('pH')



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
                print(rate,data[x][y][istep])
        except:
            Selectivity[ix][iy]=np.nan#data[x][y][istep]#/np.sum(data[x][y][:3])
            rate[ix][iy]=1e-20#np.nan#data[x][y][istep]#/np.sum(data[x][y][:nsteps])
    return rate, Selectivity

def plot_it(R,S,fig,ax,col,istep,X,Y,nsteps):
        if facet == '211': vmin=1e-3
        else: vmin=1e-6

        if nsteps == 1:
            thisax=ax
        else:
            thisax=ax[istep]
        if col == 0:
         b = thisax.imshow(R.T,
                interpolation='bicubic',
                cmap=cm.jet,
                   origin='lower', extent=[X.min(), X.max(), Y.min(), Y.max()],norm=LogNorm(),#,
                    vmin=vmin,
                    vmax=1e0,#)
                    aspect='auto')#, vmin=-abs(alldata[:,2]).max())

        else:
         a = thisax[1].imshow(S.T,
                interpolation='bicubic',
                cmap=cm.RdYlGn,
                   origin='lower', extent=[X.min(), X.max(), Y.min(), Y.max()],#norm=LogNorm(),#,
                    vmin=0,
                    vmax=1,
                    aspect='auto')#, vmin=-abs(alldata[:,2]).max())
        if istep == 1:
            fig.colorbar(b,ax=ax,shrink=1,label='*H Coverage')#,orientation='horizontal') #,location='top',orientation='horizontal')


def extract_data_from_mkm(model,dattype='coverage'):
    data={}
    pots,phs=[],[]

    if dattype=='coverage':
        datin=model.coverage_map
    elif dattype=='rate':
        datin=model.production_rate_map

    for dat in datin:
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

def run_catmaps_own_analysis(model):
        if not os.path.exists('output'):
            os.mkdir('output')

        vm = analyze.VectorMap(model)
        vm.plot_variable = 'production_rate'
        vm.log_scale = True
        vm.colorbar = True
        vm.min = 1e-5
        vm.max = 1e+2
        fig = vm.plot(save=False)
        fig.savefig('output/production_rate.pdf')

if __name__ == "__main__":
    main()
