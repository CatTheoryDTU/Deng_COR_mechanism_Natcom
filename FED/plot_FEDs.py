import pickle as pkl
import sys,os
sys.path.append('..')
import tools.plot_env
from tools.FED_tools import plot_FED_with_barrier,read_calculated_data
from tools.intermediates_dict import ads_and_electron
import numpy as np

facets=['100','211']
SHE_potentials=[-1.2]
entype='G'
pH=13
facet='100'
SHE_abs=4.4
labels=['(a)','(b)','(c)']

params={'pH':[pH],
        'potentials':[i+SHE_abs for i in SHE_potentials],
        'energy_type': entype,
        'annotate_intermediates':False,
        'annotate_reaction_conditions': False,
        'add_potential_response_in_x':False,
        'outdir': 'results/',
        'return_plt':True,
        'title':'',
        'ylim':[-1.2,1.1],
        'linestyles':['-','--',':'],
        'view': True
        }

for ifac,facet in enumerate(facets):

    if facet in ['211','111','100']:
        params['emphasize_barriers']=[['clean','H'],['H','H2_g']]
#    elif facet == '100':
#        params['emphasize_barriers']=[['clean','H'],['H','HH'],['HH','H2_g']]

    read_calculated_data(0,start_from_pkl=True,pklfile='../data/parsed_data.pckl',indict=ads_and_electron)

    # Add the gasphase boundaries - numbers on the right represent SHE potential response and energy at a workfunction of 0
    ads_and_electron['clean']['nHe']=0
    ads_and_electron['clean'][f'G_vs_pot_{facet}'] = np.array([0,0])
    ads_and_electron['clean'][f'E_vs_pot_{facet}'] = np.array([0,0])
    ads_and_electron['H2_g'][f'G_vs_pot_{facet}'] = np.array([0,0])
    ads_and_electron['H2_g'][f'E_vs_pot_{facet}'] = np.array([0,0])

    # Hack to include 2*H as an additional Volmer step from *H
    ads_and_electron['HH']={'nHe':2}
    ads_and_electron['HH'][f'G_vs_pot_{facet}'] = ads_and_electron['H'][f'G_vs_pot_{facet}']*2
    ads_and_electron['H'][f'G_ddag_vs_pot_{facet}']['HH']={'base':ads_and_electron['clean'][f'G_ddag_vs_pot_{facet}']['H']['base'].copy()+ads_and_electron['H'][f'G_vs_pot_{facet}']}
    ads_and_electron['HH'][f'G_ddag_vs_pot_{facet}']={'H2_g':{'base':ads_and_electron['H'][f'G_ddag_vs_pot_{facet}']['H2']['chemical'].copy()}}

    # Just a renaming from the parsing to connect H to the gasphase H2
    ads_and_electron['H'][f'G_ddag_vs_pot_{facet}']['H2_g']=ads_and_electron['H'][f'G_ddag_vs_pot_{facet}']['H2'].copy()
    ads_and_electron['H'][f'E_ddag_vs_pot_{facet}']['H2_g']=ads_and_electron['H'][f'E_ddag_vs_pot_{facet}']['H2'].copy()

    plt = plot_FED_with_barrier(ads_and_electron,
               facets=[facet],
               included_steps=[['clean','H','','H2_g'],
                   ['clean','H','HH','H2_g'],
                   ['clean','H']
                   ],
               proton_donor='base',
               colors=['k','b','g'],
               **params)

    Hey_y=0.85
    T_y=-1.7
    pot=SHE_potentials[0]+SHE_abs
    RHE_pot=SHE_potentials[0]+0.059*pH

    plt.annotate(labels[ifac],(-0.1,Hey_y),fontsize=40,ha='center').draggable()

    eh=ads_and_electron['H'][f'G_vs_pot_{facet}']
    plt.annotate('*H',(1,(eh[0]*pot+eh[1])+RHE_pot),fontsize=40,ha='center',va='bottom').draggable()
    eh=ads_and_electron['HH'][f'G_vs_pot_{facet}']
    plt.annotate('2*H',(2,eh[0]*pot+eh[1]+2*RHE_pot),fontsize=40,ha='center',va='bottom').draggable()
    plt.annotate('H$_{2(g)}$',(3,2*RHE_pot),fontsize=40,ha='center',va='bottom').draggable()
    eh=ads_and_electron['clean'][f'G_ddag_vs_pot_{facet}']['H']['base']
    plt.annotate('RDS',(0.5,(eh[0]*pot+eh[1])),fontsize=40,ha='center',va='bottom').draggable()

    SHE=r'V$_{\mathrm{SHE}}$'
    color='k'
    plt.annotate(f'Cu({facet}), {SHE_potentials[0]}{SHE}, pH {pH}',(-0.2,-1.1),fontsize=40,ha='left',color=color).draggable()
    plt.savefig(f'../results/FED_{facet}.pdf')
    plt.show()
    plt.close()

