import pickle as pkl
import sys,os
sys.path.append('tools_for_analysis')
from FED_tools import plot_FED_with_barrier,read_calculated_data
from scripts.intermediates_dict import ads_and_electron
from general_tools import get_reference_vibrational_contribution
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
plt.rcParams["figure.figsize"] = (10,5)
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["font.size"] = 16
plt.rc('axes', labelsize=28)    # fontsize of the x and y labels
plt.rcParams['xtick.labelsize'] = 18
plt.rcParams['ytick.labelsize'] = 18
plt.rcParams['figure.figsize'] = (10,5)
markersize=10

facets=['100','211']
SHE_potentials=[-1.2]
entype='G'
pH=13
facet='100'
SHE_abs=4.4
labels=['(a)','(b)']

data_basedir='results/parsed_data/'
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
#    fig,ax=plt.subplots(2,figsize=(8,5))
    if facet == '211':
        params['emphasize_barriers']=[['clean','H'],['H','H2_g']]
    elif facet == '100':
        params['emphasize_barriers']=[['clean','H'],['H','HH'],['HH','H2_g']]

    read_calculated_data('parsed_data.pckl',start_from_pkl=True,pklfile='parsed_data.pckl',indict=ads_and_electron)

    ads_and_electron['clean']['nHe']=0
    ads_and_electron['clean'][f'G_vs_pot_{facet}'] = np.array([0,0])
    ads_and_electron['clean'][f'E_vs_pot_{facet}'] = np.array([0,0])
    ads_and_electron['H2_g'][f'G_vs_pot_{facet}'] = np.array([0,0])
    ads_and_electron['H2_g'][f'E_vs_pot_{facet}'] = np.array([0,0])

    ads_and_electron['HH']={'nHe':2}
    ads_and_electron['HH'][f'G_vs_pot_{facet}'] = ads_and_electron['H'][f'G_vs_pot_{facet}']*2
    ads_and_electron['H'][f'G_ddag_vs_pot_{facet}']['HH']={'base':ads_and_electron['clean'][f'G_ddag_vs_pot_{facet}']['H']['base'].copy()+ads_and_electron['H'][f'G_vs_pot_{facet}']}
    ads_and_electron['HH'][f'G_ddag_vs_pot_{facet}']={'H2_g':{'base':ads_and_electron['H'][f'G_ddag_vs_pot_{facet}']['H2']['chemical'].copy()}}
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
          #     plotter=plt,
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
    plt.annotate(f'Cu({facet}), {SHE_potentials[0]}{SHE}, pH {pH}',(-0.2,-1.1),fontsize=40,ha='left').draggable()
    plt.savefig(f'results/FED_{facet}.pdf')
    plt.show()
    plt.close()

print('For catmap input:')
print(f"H 100: {ads_and_electron['H']['G_vs_pot_100'][0]*4.4+ads_and_electron['H']['G_vs_pot_100'][1]},beta={ads_and_electron['H']['G_vs_pot_100'][0]}")
print(f"H 211: {ads_and_electron['H']['G_vs_pot_211'][0]*4.4+ads_and_electron['H']['G_vs_pot_211'][1]}")
E=ads_and_electron['clean']['G_ddag_vs_pot_100']['H']['base']
print(f"H2O-ele 100: {E[0]*4.4+E[1]}, beta:{E[0]}")
E=ads_and_electron['clean']['G_ddag_vs_pot_211']['H']['base']
print(f"H2O-ele 211: {E[0]*4.4+E[1]}, beta:{E[0]}")
E=ads_and_electron['H']['G_ddag_vs_pot_100']['H2']['base']
print(f"H-H2O-ele 100: {E[0]*4.4+E[1]}, beta:{E[0]}")
E=ads_and_electron['H']['G_ddag_vs_pot_211']['H2']['base']
print(f"H-H2O-ele 211: {E[0]*4.4+E[1]}, beta:{E[0]}")
