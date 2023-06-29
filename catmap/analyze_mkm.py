#! /usr/bin/env python3
import sys, os
from catmap import analyze
sys.path.append('..')
from tools.plot_env import *
from catmaphelpers.catmap_plotter import _make_plots, _plot_model, _plot_CV, post_process_CV, _compute_tafel_prate, plot_tafel_screening, _get_unit_cell_area, _plot_pH, _plot_map, _plot_pH_bar_coverage


home=os.getcwd()

if len(sys.argv) < 2:
    print('Facet needs to be given in the command line!')
    exit()
facet=sys.argv[1]
AperSite = _get_unit_cell_area(a_input=3.6149)
ec = 1.602176e-19; A2cm = 1e-8
tof2cur = (1e3*ec/(A2cm*A2cm))/AperSite[facet]
tof2cur = {'C2+':8*tof2cur, 'Ac':4*tof2cur, 'tot':1.0,'EtOH':8*tof2cur,'C1':6*tof2cur,'HER':2*tof2cur,'C2H4':8*tof2cur}

#cov_and_rate={'coverage':['H_'+facet,'CO_'+facet],'production_rate': ['H2_g']}
cov_and_rate={'coverage':['H_'+facet],'production_rate': ['H2_g']}
facet_ratio={'111':0.5,'100':0.2,'211':0.05}
dattype='production_rate'
rateunit='TOF'
dattype='coverage'
def main():
    model=run_catmap(facet)

    steps=[i for i,label in enumerate(model.output_labels[dattype])
           if label in cov_and_rate[dattype]]
    data,phs,pots=extract_data_from_mkm(model,steps,dattype)

    fig,ax=plt.subplots(len(steps))
    plot_heatplot(data,phs,pots,steps,fig,ax,dattype=dattype)
    add_annotations(ax,len(steps))
    plt.savefig(f'../results/{dattype}_{facet}.pdf')
    plt.show()

def add_annotations(ax,nsteps):
    if nsteps == 1:
        topax,botax=ax,ax
    else:
        topax,botax=ax

    botax.set_xlabel('U vs SHE / V')
    if facet == '100':
            topax.annotate('(a) Cu(100)',(-1.7,13.5),ha='left',va='top',color='w',fontsize=40)
    elif facet == '211':
            topax.annotate('(b) Cu(211)',(-1.7,13.5),ha='left',va='top',fontsize=40)
    elif facet == '111':
            topax.annotate('(c) Cu(111)',(-1.7,13.5),ha='left',va='top',fontsize=40)
    H_eq_line=[]
    dg0={'100':0.295,'211':0.048,'111':0.197099}
    for pHs in range(4,15):
            H_eq_line.append([-dg0[facet]-0.059*pHs,pHs])
    H_eq_line=np.array(H_eq_line)
    botax.plot(H_eq_line[:6,0],H_eq_line[:6,1],'--k')
    botax.plot(H_eq_line[7:,0],H_eq_line[7:,1],'--k')
    botax.set_xlim(-1.7,-0.4)
    botax.annotate('U$_{RHE}$=-$\Delta$G$^0_\mathrm{H}$',H_eq_line[6],rotation=300,ha='center',va='center')

    if nsteps > 1:
        ax[1].annotate('*H',(-1.7,13.5),ha='left',va='top',fontsize=40,fontweight='bold')
        ax[0].annotate('*CO',(-1.7,4.5),ha='left',va='bottom',fontsize=40,fontweight='bold')

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
        if dattype == 'production_rate':
            if rateunit != 'TOF':
                R*=facet_ratio[facet]*tof2cur['HER']
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
                rate[ix][iy]=data[x][y][istep]
        except:
            Selectivity[ix][iy]=np.nan
            rate[ix][iy]=1e-20
    return rate, Selectivity

def plot_it(R,S,fig,ax,col,istep,X,Y,nsteps):
        if dattype=='coverage':
            if facet == '211': vmin=1e-3
            else: vmin=1e-8
            vmax=1
        elif dattype=='production_rate':
            vmin=1e-5
            vmax=1e3

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
                    vmax=vmax,#)
                    aspect='auto')

        else:
         a = thisax[1].imshow(S.T,
                interpolation='bicubic',
                cmap=cm.RdYlGn,
                   origin='lower', extent=[X.min(), X.max(), Y.min(), Y.max()],#norm=LogNorm(),#,
                    vmin=0,
                    vmax=1,
                    aspect='auto')

        if istep == 1:
            if dattype == 'coverage':
                label='Coverage'
            elif dattype == 'production_rate':
                if rateunit == 'TOF':
                    label = 'TOF / s$^{-1}$'
                else:
                    label='j$_{\mathrm{HER}}$ / (mA/cm$^2$)'
        if nsteps == 1:
            if dattype == 'coverage':
                label='*H Coverage'
            elif dattype == 'production_rate':
                if rateunit == 'TOF':
                    label = 'TOF / s$^{-1}$'
                else:
                    label='j$_{\mathrm{HER}}$ / (mA/cm$^2$)'

        if istep == 1 or nsteps == 1:
            fig.colorbar(b,ax=ax,shrink=1,label=label)


def extract_data_from_mkm(model,steps,dattype='coverage'):
    data={}
    pots,phs=[],[]

    if dattype=='coverage':
        datin=model.coverage_map
    elif dattype=='production_rate':
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
