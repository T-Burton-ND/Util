import cantera as ct
import argparse
import pandas as pd
import numpy as np

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-f',dest='yaml')
    #parser.add_argument('-T',dest='Temp')
    args = parser.parse_args()
    
    t_end = 100
    time_step = 1
    rule = 'css'

    # reactor setup
    gas = ct.Solution(args.yaml)
    r = ct.Reactor(contents=gas, energy='off',name='isothermal_reactor')
    #r = ct.IdealGasConstPressureReactor(contents=gas, energy='off', name='isothermal_reactor')
    sim = ct.ReactorNet([r])
    states = ct.SolutionArray(gas, extra=['t'])
    species = states.species_names
    rxns = gas.reactions()
    
    
    # use this block if you want to introduce random noise to the Ea values- could be useful for ensemble simulations
    # if you do want to use an ensemble, you would need to loop over the entire code block for each simulation
    """
    noise = np.random.normal(0,scale, len(rxns))
    temp = {}
    sim_bar = []
    # update barriers
    for cou, val in enumerate(gas.reactions()):
        temp = rxns[cou].input_data
        temp['rate-constant']['Ea'] = (rxns[cou].rate.activation_energy + noise[cou])
        gas.modify_reaction(cou, ct.Reaction.from_dict(temp, kinetics = gas))
    # if you run this block, you need to redefine the reactor (due to the modifications made to the gas and rxns pointers/objects) - it is easier to do it this way than modify the .yaml file Ea values
    r = ct.IdealGasConstPressureReactor(contents=gas, energy='off', name='isothermal_reactor')
    sim = ct.ReactorNet([r])
    states = ct.SolutionArray(gas, extra=['t'])
    """

    # cantera simulation
    states.append(r.thermo.state,t=sim.time)
    while sim.time <= t_end:
        sim.advance(sim.time + time_step)
        states.append(r.thermo.state, t=sim.time)

    # analysis of results                                                                                                                                                                                                  
    if rule == 'cnf':
        # parse depth from yaml name. If initial depth...                                                                                                                                                                                    
        if yaml_name.split("_")[-1].split('.')[0] == 0:
            # use every flux datapoint                                                                                                                                                                                                       
            spe_net = np.trapz(states.net_production_rates[:], dx=time_step, axis=0)
        # ignore first timestep                                                                                                                                                                                                              
        else: spe_net = np.trapz(states.net_production_rates[1:], dx=time_step, axis=0)
        net_states = zip(spe_net,species,range(len(species)))
        net_states = sorted(net_states,reverse=True)[:]

    elif rule == 'css':
        net_states = zip(states.X[-1,:],species,range(len(species)))
        net_states = sorted(net_states,reverse=True)[:]

    elif rule == 'final_mass':
        net_states = zip(states.Y[-1,:],species,range(len(species)))
        net_states = sorted(net_states,reverse=True)[:]

    elif rule == 'cum_conc':
        conc_net = np.trapz(states.x[:], dx=time_step, axis=0)
        net_states = zip(conc_net,species,range(len(species)))
        net_states = sorted(net_states,reverse=True)[:]

    # if you're interested in dumping to excel, this should do it (or close- I pulled this from a group member's code)
    with pd.ExcelWriter("test.xlsx") as writer:
        if rule == 'cnf': pd.DataFrame([spe_net], columns=species).to_excel(writer,sheet_name='cnf')
        elif rule == 'cum_conc': pd.DataFrame([conc_net], columns=species).to_excel(writer,sheet_name='cum_conc')
        pd.DataFrame(states.X, columns=species).to_excel(writer,sheet_name='concentration')
        pd.DataFrame(states.Y, columns=species).to_excel(writer,sheet_name='mass_percent')
        pd.DataFrame(states.net_rates_of_progress, columns=rxns).to_excel(writer,sheet_name='rxn_fluxes')
        pd.DataFrame(states.net_production_rates, columns=species).to_excel(writer,sheet_name='species_fluxes')



