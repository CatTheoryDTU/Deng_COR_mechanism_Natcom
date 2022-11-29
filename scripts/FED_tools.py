import pickle as pckl
import numpy as np
def plot_FED_with_barrier(alldata,facets,included_steps,potentials=[],pH=[13],
        ylim=[-1.2,0.6],view=False, annotate_intermediates=True,
        energy_type='G',proton_donor='base',V0_SHE=4.4,normalize_to_IS=False,
        figsize=(15,5),title=None,annotate_reaction_conditions=True,
        check_barriers_with_line=False,outdir='results',outformat='png',
        labelsize=None,colors=['k','b','r','g','y'],emphasize_barriers=[],
        linestyles=['-','--','-.',':'],return_plt=False,xticks=None,
        add_potential_response_in_x=False,xlabel='N$_\mathrm{H}$ + $\gamma$',
        plotter=None):

    if plotter is not None:
        plt=plotter
    else:
        from matplotlib import pyplot as plt
    plt.rcParams["figure.figsize"] = figsize


    if not potentials:
        print('No potentials where given for FED with barriers')
        return

    if isinstance(included_steps[0],str):
        included_steps=[included_steps]

    if isinstance(facets,str):
        facets=[facets]

    if len(facets) > 1:
        # If more than one facet is included the facets have different colors
        # and the mechanism different linestyles
         colors=list(reversed(colors[:len(facets)]))
         linestyles=list(reversed(linestyles[:len(included_steps)]))
         print('Several facets have been given, they will have varying '
               'colors and the mechanism will change linestyle')

    else:
        # If only one facet is included the mechanism have different colors
        # and the potentials different linestyles
         colors=list(reversed(colors[:len(included_steps)]))
         #linestyles=list(linestyles[0])*len(included_steps)
         linestyles=list(linestyles)*len(included_steps)
         if len(list(linestyles)) >= len(potentials):
             linestyles=list(linestyles)*len(included_steps)
         print('Only one facet has been given the mechanisms will have varying '
               'colors and the potentials will change linestyle')

    _xticks=[]
    all_ens=[]
    lowest_xval,longest_mech=100,0
    for ifac,facet in enumerate(facets):
     #print(len(facets),ifac,colors)
     enstring='%s_vs_pot_%s'%(energy_type,facet)
     barenstring='%s_ddag_vs_pot_%s'%(energy_type,facet)
     for imech,single_mech_steps_in in enumerate(included_steps):
        if len(single_mech_steps_in) > longest_mech: longest_mech=len(single_mech_steps_in)
        #If facets are compared color by facet otherwise by mechanism
        if len(facets) > 1:
            color = colors[ifac]
    #        print(color)
        elif len(included_steps) > 1:
            color = colors[imech]
        else:
            color=colors[0]

        single_mech_steps=single_mech_steps_in
        added_ads=[None]*len(single_mech_steps_in)
        #If C2  is plotted in the same plot as C1 intermediates the name is given
        # with a plus and will be split here
        if np.any([['+' in i for i in single_mech_steps_in]]):
            single_mech_steps = [i.split('+')[0] for i in single_mech_steps_in]
            for istep,step in enumerate(single_mech_steps_in):
                if step.split('+')[-1] != single_mech_steps_in[istep]:
                    added_ads[istep] = step.split('+')[-1]
                else:
                    added_ads.append(None)

        #Recognize if adsorbate is given with a 2 infront e.g. 2CO
        if np.any([[i[0] == '2' for i in single_mech_steps_in if len(i)]]):
            single_mech_steps = [i.lstrip('2') for i in single_mech_steps]
            for istep,step in enumerate(single_mech_steps_in):
                if not len(step): continue
                if step[0] == '2':#single_mech_steps[istep]:
                    added_ads[istep] = step.lstrip('2')#split('+')
                else:
                    added_ads.append(None)

        if enstring not in alldata[single_mech_steps[0]].keys():
            add_CHE_and_energy_vs_RHE(alldata,facet,pH=pH)
            #add_vibrational_contribution(alldata,facet,barrier=[])


        for iph,ph in enumerate(pH):
         if isinstance(proton_donor,dict):
            pdonor=proton_donor[ph]
         else:
             pdonor=proton_donor


         for ipot,potential in enumerate(potentials):
            rhe_pot=potential-(V0_SHE-0.059*ph)
            IS_normalization=0
            for istep,current_intermediate in enumerate(single_mech_steps):
#                current_intermediate=single_mech_steps[istep]
                if istep < len(single_mech_steps)-1:
                    next_intermediate,inext_step=single_mech_steps[istep+1],istep+1
                    if next_intermediate == '':
                        next_intermediate=single_mech_steps[istep+2]
                        inext_step=istep+2
                if current_intermediate not in alldata.keys():
                    print('Could not find the intermediate ',current_intermediate)
                    continue
                if enstring not in alldata[current_intermediate]:
                    print('Intermdiate %s does not seem to have an energy with name %s'%(current_intermediate,enstring))
                    continue

#                print(alldata[current_intermediate].keys(),current_intermediate)
                _xticks.append(int(alldata[current_intermediate]['nHe']))

                En_IS=np.poly1d(alldata[current_intermediate][enstring])(potential)
                #Add CHE to endstate
                En_IS+=alldata[single_mech_steps[istep]]['nHe']*rhe_pot

                if added_ads[istep] is not None:
                    En_IS+=np.poly1d(alldata[added_ads[istep]][enstring])(potential)+\
                            alldata[added_ads[istep]]['nHe']*rhe_pot

                if normalize_to_IS and istep==0:
                    IS_normalization=-En_IS

                if current_intermediate == 'CO' and 'OCCO' in single_mech_steps:
                    if barenstring not in alldata['CO']:
                        alldata['CO'][barenstring] = {}

                    if barenstring in alldata['CO']:
                        if 'OCCO' not in alldata['CO'][barenstring].keys():
                            alldata['CO'][barenstring]['OCCO']= alldata['COCO'][barenstring]['OCCO']



                if annotate_intermediates:
                    stepout=current_intermediate
                    if added_ads[istep] is not None:
                        stepout+='+'+added_ads[istep]
                        if current_intermediate == added_ads[istep]:
                            stepout = '2'+current_intermediate

                    if labelsize is None:
                        plt.annotate(stepout,xy=(inext_step,ylim[1]-0.1*(imech+1)),color=color).draggable()
                    else:
                        plt.annotate(stepout,xy=(inext_step,ylim[1]-0.1*(imech+1)),color=color,fontsize=labelsize).draggable()

                #if istep == len(single_mech_steps)-1:
                #    plt.plot([istep+0.75,istep+1.25],[En_IS+IS_normalization,En_IS+IS_normalization],'-'+colors[imech])
                #    continue
                #print(istep,current_intermediate)
                if istep < len(single_mech_steps)-1:
                    if next_intermediate == '':
                        print('Skipping the step %i'%inext_step)
                    elif next_intermediate not in alldata.keys():
                        print('Could not find the intermediate for FS:',next_intermediate)
                        continue
                    else:
                        if enstring not in alldata[next_intermediate]:
                            print('Intermdiate %s_%s does not seem to have an energy'%(next_intermediate,facet))
                            continue

                    print(current_intermediate,next_intermediate)
                    En_FS=np.poly1d(alldata[next_intermediate][enstring])(potential)+\
                            alldata[next_intermediate]['nHe']*rhe_pot

                    if added_ads[inext_step] is not None:
                        En_FS+=np.poly1d(alldata[added_ads[inext_step]][enstring])(potential)+\
                                        alldata[added_ads[inext_step]]['nHe']*rhe_pot


                    #If no  barrier at all has  been calculated from the current intermediate
                    #Draw a straight line to the next intermediate
                    Eddag=None
                    if barenstring not in alldata[current_intermediate].keys():
                        #XXX: IS the arrow below ever needed?
#                        plt.arrow(istep+1.25,En_IS,0.5,En_FS-En_IS,head_width=0.0,length_includes_head=True, linewidth=0.1,color='r')
                        pass
                    #For multiple H adsorptions
                    elif next_intermediate in ['2H'] and current_intermediate in ['H']:
                        Eddag=np.poly1d(alldata['clean'][barenstring][next_intermediate]['base'])(potential)+alldata['H'][enstring]

                    #For a chemical step
                    elif next_intermediate in ['H','OCCO','CO2']:
                        #Eddag=np.poly1d(alldata['clean'][barenstring]['H']['base'])(potential)
                        Eddag=np.poly1d(alldata[current_intermediate][barenstring][next_intermediate]['base'])(potential)
                    #If the barrier between IS and FS has been calculated
                    elif next_intermediate in alldata[current_intermediate][barenstring]:
                     toads=next_intermediate
                     if pdonor in alldata[current_intermediate][barenstring][toads]:
                        Eddag=np.poly1d(alldata[current_intermediate][barenstring][toads][pdonor])(potential)


                    #print(facet,Eddag,single_mech_steps[istep])

                     #Add CHE to barriers (alkaline has CHE like IS, acid has CHE like FS)
                    if Eddag is not None:
                        Eddag+=alldata[current_intermediate]['nHe']*rhe_pot
                        if added_ads[istep] is not None and added_ads[inext_step] is not None:
                                Eddag+=np.poly1d(alldata[added_ads[istep]][enstring])(potential)


                        if pdonor == 'acid' and next_intermediate not in ['OCCO','CO2']:
                                Eddag+=rhe_pot

                        #Connect the states by lines (TODO: Maybe polynoms?)

                #plot everything
                #Plot IS
                if len(facets) > 1:
                    linestyle=linestyles[imech]
                elif len(potentials) > 1:
                    linestyle=linestyles[ipot]
                elif len(pH):
                    linestyle=linestyles[iph]

                current_xval,next_xval=istep,inext_step
                if add_potential_response_in_x:
                    current_xval=alldata[current_intermediate][enstring][0]+alldata[current_intermediate]['nHe']
                    next_xval=alldata[next_intermediate][enstring][0]+alldata[next_intermediate]['nHe']

                if current_xval < lowest_xval:
                    lowest_xval=current_xval

                #print(current_intermediate,current_xval)
#                plt.plot([istep+0.75,istep+1.25],
                plt.plot([current_xval-0.25,current_xval+0.25],
                        [En_IS+IS_normalization,En_IS+IS_normalization],linestyle=linestyle,
                        color=color,linewidth=4)
                all_ens.append(En_IS+IS_normalization)
                if istep == len(single_mech_steps)-1:
                    continue

                all_ens.append(En_FS+IS_normalization)

                #If the barrer to the FS has not been calculated
                #Draw a straight line to the next intermediate
                arrow_xlen=0.5+(next_xval-current_xval-1)
                straight_linestyle=linestyle
#                print(current_intermediate,barenstring,alldata[current_intermediate].keys())
#                print(barenstring,alldata[current_intermediate][barenstring])
                Barrier_found=0
                if linestyle=='--': straight_linestyle=':' #Hack because -- looks like a solid line
                if barenstring in alldata[current_intermediate]:
                 if next_intermediate in alldata[current_intermediate][barenstring]:# or next_intermediate == 'H':
                  if Eddag is not None:
                    Barrier_found=1
                    x_IS=current_xval+0.25
                    x_FS=next_xval-0.25
                    barx=np.sqrt((En_IS-Eddag)*(En_FS-Eddag)*(x_IS-x_FS)**2)+En_IS*x_FS-En_FS*x_IS+Eddag*(x_IS-x_FS)
                    barx/=En_IS-En_FS
                    if np.isnan(barx):
                        txtout=f'WARNING! Automatic detection of barrier x '
                        txtout+=f'for {current_intermediate} to '
                        txtout+=f'{next_intermediate} at {potential}V '
                        txtout+=f'failed, likely you are '
                        txtout+=f'barrierless. Check with check_barriers_with_line keyword.'
                        print(txtout)
                        barx=current_xval+(next_xval-current_xval)/2+0.08*(En_FS-En_IS)
                    #if add_potential_response_in_x:
                    # x_IS=current_xval+0.25
                    # x_FS=next_xval-0.25
                    # barx=np.sqrt((En_IS-Eddag)*(En_FS-Eddag)*(x_IS-x_FS)**2)+En_IS*x_FS-En_FS*x_IS+Eddag*(x_IS-x_FS)
                    # barx/=En_IS-En_FS
                    parpts=np.array([[current_xval+0.25,En_IS],
                     #   #[current_xval+(next_xval-current_xval)/2+0.04*(En_FS-En_IS),Eddag], # x-value should b/(2a) of the parabola function
                        [barx,Eddag], # x-value should b/(2a) of the parabola function
                        [next_xval-0.25,En_FS]])
                    #else:
                    # barx=np.sqrt((En_IS-Eddag)*(En_FS-Eddag)*(x_IS-x_FS)**2)+En_IS*x_FS-En_FS*x_IS+Eddag*(x_IS-x_FS)
                    # barx/=En_IS-En_FS
                    # if np.isnan(barx):
                     #    barx=current_xval+(next_xval-current_xval)/2+0.08*(En_FS-En_IS)
                     #print(barx)
                     #parpts=np.array([[current_xval+0.25,En_IS],
                     #   [barx,Eddag],
#                   #     [current_xval+(next_xval-current_xval)/2+0.08*(En_FS-En_IS),Eddag], # x-value should b/(2a) of the parabola function
                     #   [next_xval-0.25,En_FS]])

                    from general_tools import quad_fun
                    from scipy.optimize import curve_fit
                    coeff,dummy=curve_fit(quad_fun,parpts[:,0],parpts[:,1])
                    fitpts=np.linspace(current_xval+0.25,next_xval-0.25,30)#,include_endpoints=True)
                    fit=np.poly1d(coeff)(fitpts)+IS_normalization

                    # Emphasize chosen barriers based on emphasize_barrier input
                    emphasize=False
                    linewidth=1
                    if len(emphasize_barriers):
                        if isinstance(emphasize_barriers,dict):
                            for prop in emphasize_barriers:
                                if prop[:2] == 'pH':
                                    if abs(ph - float(prop[2:])) < 0.01:
                                        pairs=emphasize_barriers[prop]
                                else:
                                    print('pH given for emphasizing barriers not understood')
                        elif isinstance(emphasize_barriers,list):
                            pairs=emphasize_barriers

                        for pair in pairs:
                                if current_intermediate in pair and next_intermediate in pair:
                                    linewidth=6

                    # Plot the parabolic barriers
                    # XXX color of acidic needs to be updated
                    if pdonor == 'acid':
                            plt.plot(fitpts,fit,'--b',linewidth=linewidth)
#                            plt.plot([istep+1.35,istep+1.65],[Eddag+IS_normalization,Eddag+IS_normalization],'--b')
                    else:
                            #plt.plot(fitpts,fit,linestyles[imech],color=color)
                            plt.plot(fitpts,fit,linestyle,color=color,linewidth=linewidth)
                            if check_barriers_with_line:
                               plt.plot([next_xval-0.65,next_xval+0.35],[Eddag+IS_normalization,Eddag+IS_normalization],'--k') #+colors[imech])

                    all_ens.append(Eddag+IS_normalization)

                if not Barrier_found:
                    print(f'Barrier for step {current_intermediate} to {next_intermediate} in {pdonor} not found')
                    # Awful hack regarding the color for SelectCO2 deliverable
                    linecolor=colors[imech]
                    if pdonor == 'acid':
                            linecolor='b'
                    plt.arrow(current_xval+0.25,En_IS+IS_normalization,arrow_xlen,
                                En_FS-En_IS,head_width=0.0,length_includes_head=True,
                                linewidth=1.5,color=linecolor,linestyle=straight_linestyle)



    plt.ylim(ylim)
    if not add_potential_response_in_x:
        plt.xlim([-0.25,longest_mech-lowest_xval-0.75])
    else:
        #print(longest_mech-lowest_xval)
        plt.xlim([lowest_xval-0.5,longest_mech-lowest_xval+2.5])
        plt.xlabel(xlabel)


    if xticks is None:
        plt.xticks([])
    else:
        plt.xticks(np.unique(_xticks))
    if proton_donor=='base':
        donor='H$_2$O'
    elif proton_donor == 'acid':
        donor='H$_3$O$^+$'
    elif proton_donor == 'chemical':
        donor='chemical'

    elif isinstance(proton_donor,dict):
        donor=list(proton_donor.items())
#    print(','.join(['-'.join(i) for i in included_steps]))
    if title is None:
        plt.title(#','.join(['-'.join(i) for i in included_steps])+
            ', facet: '+facet+
            ', WF='+'-'.join([str(np.around(i,2)) for i in potentials])+
            ', pH='+'-'.join([str(np.around(i,1)) for i in pH])+
            ', proton donor: %s '%donor,fontsize=15)
        #'-'.join(included_steps))+
    else:
        plt.title(title)

    plt.ylabel(f'$\Delta${energy_type}$^\phi$ / eV')

    if annotate_reaction_conditions:
        sheout=','.join([str(np.around(i-V0_SHE,2)) for i in potentials])+'V$_{\mathrm{SHE}}$'
        phout = ','.join([str(np.around(i,1)) for i in pH])
        plt.annotate(f"Cu({facet})\n{sheout}\npH={phout}",(0.7,ylim[0]),ha='left',va='bottom',fontsize=23).draggable()

    plt.tight_layout()
    if view:
        if return_plt:
            return plt
        else:
            plt.show()
            return
    if potentials:
        #print(included_steps)
        if isinstance(included_steps[0],list):
            stepsout=['-'.join(steps) for steps in included_steps]
        else:
            stepsout=included_steps
        #print(stepsout)
        plt.savefig(outdir+'/FED_w_barrier_'+
                '_'.join(stepsout)+
                '_pot_'+'-'.join([str(np.around(i,2)) for i in potentials])+
                '_pH_'+'-'.join([str(np.around(i,2)) for i in pH])+
                '.'+outformat,transparent=True)
    else:
        plt.savefig(outdir+'/FED_w_barrier_'+'-'.join(included_steps)+'.'+outformat,transparent=True)
    plt.close()

def read_calculated_data(inputfile,facets=None,start_from_pkl=False,
                         indict=None,substrates=['Cu'],pklfile='results/parsed_data.pckl',
                         add_field=False,field_data=None,V0SHE=4.4,PZC=None,Capacitance=None):

    if not indict:
            sys.path.append('/Users/geokast/SelectCO2/endstates')
            from intermediates_dict import ads_and_electron

    else:
        ads_and_electron=indict
        #ads_and_electron={}

    if start_from_pkl:
        print('Reading data from %s'%pklfile)
        alldata=pckl.load(open(pklfile,'rb'))
        #if not indict:
        for ads in alldata.keys():
                ads_short=ads.lstrip('md-').lstrip('bdo-')
                if ads_short in ads_and_electron.keys():
                 ads_and_electron[ads_short].update(alldata[ads])
        return

    print('Reading data from %s'%inputfile)

    inlines = open(inputfile,'r').readlines()[1:]
    if isinstance(facets,str):
        facets=[facets]

    #Set up dictionary
    for facet in facets:
      for line in inlines:
         ads = line.split()[2]
         if line.split()[0] in substrates+['None']:
            if line.split()[1] == 'gas':
                ads+='_g'

            if ads in ['CO2_g','CH3COOH_g','CH2CO_g','OH_g']: continue
            if not indict and ads not in ads_and_electron.keys():
                    ads_and_electron[ads] = {}

      if add_field:
        collect_field_data(ads_and_electron,facet,field_data)

    for facet in facets:
      for line in inlines:
         ads = line.split()[2]
         if line.split()[0] not in substrates+['None']: continue
         if line.split()[1] not in [facet,'gas']: continue

         if line.split()[1] == 'gas':
                ads+='_g'

                if ads in ads_and_electron.keys():
                #if ads in ['CO2_g','CH3COOH_g','CH2CO_g','OH_g']: continue
                    ads_and_electron[ads]['E'] = float(line.split()[3])

         #elif line.split()[1] in [facet,'gas']:
         elif ads in ads_and_electron.keys():
             ads_and_electron[ads]['E_%s'%facet] = float(line.split()[3])
             #TODO: Here add field potential dependence!
             if add_field:
                 if not PZC or not Capacitance:
                     print('PZC and Capacitance have to be given to the read function!')
                     ads_and_electron[ads]['E_vs_pot_%s'%facet] = np.array([0,ads_and_electron[ads]['E_%s'%facet]])
                 else:
                     dedphi = ads_and_electron[ads]['dedq'][facet]*Capacitance
                     if PZC < 3: PZC=PZC+V0SHE
                     offset = ads_and_electron[ads]['E_%s'%facet]-PZC*dedphi
                     #print(ads,dedphi)
                     ads_and_electron[ads]['E_vs_pot_%s'%facet] = np.array([dedphi,offset])

             else:
                 ads_and_electron[ads]['E_vs_pot_%s'%facet] = read_beta_from_catmap_input(line)

         if line.split()[1] in [facet,'gas']:
             if ads in ads_and_electron.keys():
                 freq_inline=[None,None]
                 for isplit,splitline in  enumerate(line.split()):
                     if splitline[0] == '[':
                         freq_inline[0]=isplit
                     elif splitline[-1] == ']':
                         freq_inline[1]=isplit+1
                         break

                 if None not in freq_inline:
                     frequencies = [float(vib.replace(',','').replace('[','').replace(']',''))
                             for vib in line.split()[freq_inline[0]:freq_inline[1]]]

                 else:
                     print('No frequencies given for '+ads)
                     frequencies=[]

                 ads_and_electron[ads]['vibs_%s'%facet] = frequencies

             elif '-' in ads:
                 barads=ads.split('-')[0]
                 if barads not in ads_and_electron.keys() and barads != 'OC': continue
                 if ads.split('-')[1] == 'ele':
                     if barads == 'COCO':
                         bartoads='OCCO'

                     elif barads.replace('O','',1).replace('H','',1) in ads_and_electron.keys():
                         bartoads=barads.replace('O','',1).replace('H','',1)
                     elif barads[::-1].replace('O','',1).replace('H','',1) in ads_and_electron.keys():
                         bartoads=barads[::-1].replace('O','',1).replace('H','',1)
                         bartoads=bartoads[::-1]
                     else:
                         print('Could not determine the product of %s'%ads)
                         bartoads=ads

                 elif ads.split('-')[1] == 'H2O':
                     if 'H'+barads in ads_and_electron.keys():
                         bartoads='H'+barads
                     elif barads == 'OC':
                         barads,bartoads='CO','CHO'
                     elif barads == 'CO':
                         barads,bartoads='CO','COH'
                     elif barads+'H' in ads_and_electron.keys():
                         bartoads=barads+'H'
                     else:
                         print('Could not determine the product of %s'%ads,barads)
                         bartoads=ads

                 else:
                     print(ads,' special case check')
                     bartoads=barads[::-1]+barads


                 if 'E_ddag_%s'%facet not in ads_and_electron[barads].keys():
                     ads_and_electron[barads]['E_ddag_%s'%facet]={}
                     ads_and_electron[barads]['E_ddag_vs_pot_%s'%facet]={}
                     ads_and_electron[barads]['vibs_ddag_%s'%facet]={}
                 ads_and_electron[barads]['E_ddag_%s'%facet][bartoads] = float(line.split()[3])
                 ads_and_electron[barads]['E_ddag_vs_pot_%s'%facet][bartoads] =                                read_beta_from_catmap_input(line)

                 freq_inline=[None,None]
                 #da
                 for isplit,splitline in  enumerate(line.split()):
                     if splitline[0] == '[':
                         freq_inline[0]=isplit
                     elif splitline[-1] == ']':
                         freq_inline[1]=isplit+1
                         break

                 if None not in freq_inline:
                     frequencies = [float(vib.replace(',','').replace('[','').replace(']',''))
                             for vib in line.split()[freq_inline[0]:freq_inline[1]]]
                 #ads_and_electron[barads]['vibs_ddag_%s'%facet][bartoads] = read_vibrational_frequencies(ads,   line,vibfile,facet)

                 ads_and_electron[barads]['vibs_ddag_%s'%facet][bartoads] = frequencies
             else:
                 print(ads+' should be added to the dict')
    return ads_and_electron


#main()
