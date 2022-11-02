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

def plotATPBudget(Model,solutions,thresh=1e-5,displayModelID=0):
    ATP_rxns=sorted([x for x in Model.metabolites.atp_c.reactions if any([abs(sol[x.id])>thresh for sol in solutions])],key=lambda x: abs(solutions[0][x.id]*x.get_coefficient('atp_c')),reverse=True)
    cm = plt.get_cmap('tab20')
    posSum=np.zeros(len(solutions))
    negSum=np.zeros(len(solutions))
    width = 0.7

    thisplot = plt.figure(figsize=(3,4))
    plts=[]
    for it,rxn in enumerate(ATP_rxns):
        rxnFluxes=[sol[rxn.id]*rxn.get_coefficient('atp_c') for sol in solutions]
        plts+=[plt.bar(np.arange(len(solutions)), rxnFluxes, bottom=[negSum[x] if rxnFluxes[x]<0 else posSum[x] for x in range(len(rxnFluxes))],color=cm(it/len(ATP_rxns)),width=width)]
        negSum=[negSum[x]+rxnFluxes[x] if rxnFluxes[x]<0 else negSum[x] for x in range(len(negSum))]
        posSum=[posSum[x]+rxnFluxes[x] if rxnFluxes[x]>0 else posSum[x] for x in range(len(posSum))]

    plt.rcParams.update({'font.size': 10}) #sets a global fontsize
    plt.rcParams['xtick.major.size'] = 5 # adjusts tick line length and width
    plt.rcParams['xtick.major.width'] = 1
    plt.rcParams['ytick.major.size'] = 5
    plt.rcParams['ytick.major.width'] = 1
    plt.rcParams['axes.linewidth']=2 # makes axes line thicker
    plt.xlim(-0.6,len(solutions)-0.4)
    plt.ylabel("ATP produced/consumed (µmol/s)")
    handles, labels = plt.gca().get_legend_handles_labels()
    plt.axhline(0,linestyle="--",color="black")
    plt.tight_layout
    if displayModelID:
        plt.legend(plts,[x.id for x in ATP_rxns],bbox_to_anchor=(1,1))
    else:
        plt.legend(plts,[x.name for x in ATP_rxns],bbox_to_anchor=(1,1))
    plt.show()    

def plotNADBudget(Model,solutions,thresh=1e-5,displayModelID=0):
    NADH_rxns = list(Model.metabolites.nadh_c.reactions)
    NADPH_rxns = list(Model.metabolites.nadph_c.reactions)
    nadh = Model.metabolites.nadh_c
    nadph = Model.metabolites.nadph_c
    NAD_rxns=sorted([x for x in NADH_rxns+NADPH_rxns if any([abs(sol[x.id])>thresh for sol in solutions])],key=lambda x: [abs(solutions[0][x.id]*x.get_coefficient(nadh)) if nadh in x.metabolites else abs(solutions[0][x.id]*x.get_coefficient(nadph))],reverse=True)
    cm = plt.get_cmap('tab20')
    
    posSum=np.zeros(len(solutions))
    negSum=np.zeros(len(solutions))
    width = 0.7

    thisplot = plt.figure(figsize=(3,4))
    plts=[]
    for it,rxn in enumerate(NAD_rxns):
        nadhSol=np.zeros(len(solutions))
        nadphSol=np.zeros(len(solutions))
        if nadh in rxn.metabolites:
            nadhSol=[sol[rxn.id]*rxn.get_coefficient('nadh_c') for sol in solutions]
        if nadph in rxn.metabolites:
            nadphSol=[sol[rxn.id]*rxn.get_coefficient('nadph_c') for sol in solutions]
        rxnFluxes=[nadhSol[x]+nadphSol[x] for x in range(len(solutions))]
        plts+=[plt.bar(np.arange(len(solutions)), rxnFluxes, bottom=[negSum[x] if rxnFluxes[x]<0 else posSum[x] for x in range(len(rxnFluxes))],color=cm(it/len(NAD_rxns)),width=width)]
        negSum=[negSum[x]+rxnFluxes[x] if rxnFluxes[x]<0 else negSum[x] for x in range(len(negSum))]
        posSum=[posSum[x]+rxnFluxes[x] if rxnFluxes[x]>0 else posSum[x] for x in range(len(posSum))]

    plt.rcParams.update({'font.size': 10}) #sets a global fontsize
    plt.rcParams['xtick.major.size'] = 5 # adjusts tick line length and width
    plt.rcParams['xtick.major.width'] = 1
    plt.rcParams['ytick.major.size'] = 5
    plt.rcParams['ytick.major.width'] = 1
    plt.rcParams['axes.linewidth']=2 # makes axes line thicker
    plt.xlim(-0.6,len(solutions)-0.4)
    plt.ylabel("NAD(P)H produced/consumed (µmol/s)")
    handles, labels = plt.gca().get_legend_handles_labels()
    plt.axhline(0,linestyle="--",color="black")
    plt.tight_layout
    if displayModelID:
        plt.legend(plts,[x.id for x in NAD_rxns],bbox_to_anchor=(1,1))
    else:
        plt.legend(plts,[x.name for x in NAD_rxns],bbox_to_anchor=(1,1))
    plt.show()

    
