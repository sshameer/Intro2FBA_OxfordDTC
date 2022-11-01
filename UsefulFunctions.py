import matplotlib.pyplot as plt
import numpy as np
def plotBoundary(Model, solutions,thresh=1e-5,displayModelID=0,solutionNames=[]):
    bdryRxns = sorted([x for x in Model.boundary if any([abs(sol[x.id])>thresh for sol in solutions])],key=lambda x:abs(solutions[0][x.id]),reverse=True)
    width = 0.9/len(solutions)
    plt.ylabel('flux')
    for it,sol in enumerate(solutions):
        plt.bar(np.arange(len(bdryRxns))+(it)*width,[sol[x.id] for x in bdryRxns],width=width)
    if displayModelID:
        plt.xticks(range(len(bdryRxns)),[x.id for x in bdryRxns],rotation=90)
    else:
        plt.xticks(range(len(bdryRxns)),[x.name for x in bdryRxns],rotation=90)
    if solutionNames and len(solutionNames)==len(solutions):
        plt.legend(solutionNames)
    else:
        if solutionNames and len(solutionNames)!=len(solutions):
            print('Length of solutionNames does not match length of solutions.')
        plt.legend(['Solution '+str(x) for x in range(len(solutions))])
def plotATPBudget(Model,sol,thresh=1e-5,displayModelID=0):
    ATP_rxns=sorted([x for x in Model.metabolites.atp_c.reactions if abs(sol[x.id])>thresh],key=lambda x: abs(sol[x.id]*x.get_coefficient('atp_c')),reverse=True)
    cm = plt.get_cmap('tab20')
    # fig, [ax1, ax2] = plt.subplots(nrows=1, ncols=2)
    posSum=0
    negSum=0
    width = 0.1

    thisplot = plt.figure(figsize=(3,4))
    plts=[]
    for it,rxn in enumerate(ATP_rxns):
        if sol[rxn.id]*rxn.get_coefficient('atp_c')<0:
            plts+=[plt.bar(1, sol[rxn.id]*rxn.get_coefficient('atp_c'), bottom=negSum,color=cm(it/len(ATP_rxns)),width=width)]
            negSum+=sol[rxn.id]*rxn.get_coefficient('atp_c')
        else:
            plts+=[plt.bar(1, sol[rxn.id]*rxn.get_coefficient('atp_c'), bottom=posSum,color=cm(it/len(ATP_rxns)),width=width)]
            posSum+=sol[rxn.id]*rxn.get_coefficient('atp_c')

    plt.rcParams.update({'font.size': 10}) #sets a global fontsize
    plt.rcParams['xtick.major.size'] = 5 # adjusts tick line length and width
    plt.rcParams['xtick.major.width'] = 1
    plt.rcParams['ytick.major.size'] = 5
    plt.rcParams['ytick.major.width'] = 1
    plt.rcParams['axes.linewidth']=2 # makes axes line thicker
    plt.xlim(0.8,1.2)
    plt.ylabel("ATP produced/consumed (µmol/s)")
    handles, labels = plt.gca().get_legend_handles_labels()
    plt.axhline(0,linestyle="--",color="black")
    plt.tight_layout
    if displayModelID:
        plt.legend(plts,[x.id for x in ATP_rxns],bbox_to_anchor=(1,1))
    else:
        plt.legend(plts,[x.name for x in ATP_rxns],bbox_to_anchor=(1,1))
    plt.show()    
def plotNADBudget(Model,solution,thresh=1e-5,displayModelID=0):
    NADH_rxns = list(Model.metabolites.nadh_c.reactions)
    NADPH_rxns = list(Model.metabolites.nadph_c.reactions)
    nadh = Model.metabolites.nadh_c
    nadph = Model.metabolites.nadph_c
    NAD_rxns=sorted([x for x in NADH_rxns+NADPH_rxns if abs(sol[x.id])>thresh],key=lambda x: [abs(sol[x.id]*x.get_coefficient(nadh)) if nadh in x.metabolites else abs(sol[x.id]*x.get_coefficient(nadph))],reverse=True)
    cm = plt.get_cmap('tab20')
    # fig, [ax1, ax2] = plt.subplots(nrows=1, ncols=2)
    posSum=0
    negSum=0
    width = 0.1

    thisplot = plt.figure(figsize=(3,4))
    plts=[]
    for it,rxn in enumerate(NAD_rxns):
        nadhSol=0
        nadphSol=0
        if nadh in rxn.metabolites:
            nadhSol=sol[rxn.id]*rxn.get_coefficient('nadh_c')
        if nadph in rxn.metabolites:
            nadphSol=sol[rxn.id]*rxn.get_coefficient('nadph_c')
        thisSol = nadhSol+nadphSol 
        if thisSol<0:
            plts+=[plt.bar(1, thisSol, bottom=negSum,color=cm(it/len(NAD_rxns)),width=width)]
            negSum+=thisSol
        else:
            plts+=[plt.bar(1, thisSol, bottom=posSum,color=cm(it/len(NAD_rxns)),width=width)]
            posSum+=thisSol

    plt.rcParams.update({'font.size': 10}) #sets a global fontsize
    plt.rcParams['xtick.major.size'] = 5 # adjusts tick line length and width
    plt.rcParams['xtick.major.width'] = 1
    plt.rcParams['ytick.major.size'] = 5
    plt.rcParams['ytick.major.width'] = 1
    plt.rcParams['axes.linewidth']=2 # makes axes line thicker
    plt.xlim(0.8,1.2)
    plt.ylabel("NAD(P)H produced/consumed (µmol/s)")
    handles, labels = plt.gca().get_legend_handles_labels()
    plt.axhline(0,linestyle="--",color="black")
    plt.tight_layout
    if displayModelID:
        plt.legend(plts,[x.id for x in NAD_rxns],bbox_to_anchor=(1,1))
    else:
        plt.legend(plts,[x.name for x in NAD_rxns],bbox_to_anchor=(1,1))
    plt.show()

    
