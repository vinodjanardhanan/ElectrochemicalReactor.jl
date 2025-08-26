
function interface_kinetics!(du , u, p)
    
    T = p[:T]
    f = IdealGas.F/IdealGas.R/T
    
    cell = p[:cell]
    eco = cell.eChem    
    eParams = eco.ecp

    Ω = p[:res]    

    cd = u[1]
    ηa = u[2]
    ηc = u[3]
    

    
    αoa = eParams.α_ox_anode
    αra = eParams.α_rd_anode
    αoc = eParams.α_ox_cathode
    αrc = eParams.α_rd_cathode


    # du[1] = eParams.Ecell - (eco.Erev - ηa - abs(ηc) - ηΩ)
    du[1] = potential_balance(cell, u, Ω, eco)
    du[2]  = cd - eParams.i0a * (exp(αoa*f*ηa) - exp(-αra*f*ηa))    
    du[3]  = cd + eParams.i0c * (exp(αoc*f*ηc) - exp(-αrc*f*ηc))    

end

function potential_balance(p::FuelCell, u, Ω, eco::ElectrochemObject)
    cd = u[1]
    ηa = u[2]
    ηc = u[3]
    ηΩ = Ω * cd    
    return eco.ecp.Ecell - (eco.Erev - ηa - abs(ηc) - ηΩ)
end

function potential_balance(p::ElectrolysisCell, u,Ω, eco::ElectrochemObject)
    cd = u[1]
    ηa = u[2]
    ηc = u[3]
    ηΩ = Ω * cd
    return eco.ecp.Ecell - (eco.Erev + ηa + abs(ηc) + ηΩ)
end


function electrode_reaction(p::SOFC_H2) 
    T = p.anode.ch.gp.T
    pH2 = p.anode.conc[end][p.iH2]*IdealGas.R*T
    pH2O = p.anode.conc[end][p.iH2O]*IdealGas.R*T
    pO2 = p.cathode.conc[end][p.iO2]*IdealGas.R*T


    p.eChem.Erev = nernst_potential(H2Oxidation(), p.eChem.Estd, T, aH2=pH2/IdealGas.p_std, aH2O = pH2O/IdealGas.p_std, aO2=pO2/IdealGas.p_std)
    # p.eChem.Erev = nernst_h2(p.eChem.Estd, T, pH2=pH2, pH2O = pH2O, pO2=pO2)
    # The following is a call to a user defined function to calculate the exchange current density
    p.eChem.ecp = p.eChem.udf((pH2=pH2, pH2O = pH2O, pO2=pO2, T = T))        
    solve_electrochemistry(p)    
end

function electrode_reaction(p::SOEC_H2) 
    T = p.anode.ch.gp.T
    pH2 = p.cathode.conc[end][p.iH2]*IdealGas.R*T
    pH2O = p.anode.conc[end][p.iH2O]*IdealGas.R*T
    pO2 = p.anode.conc[end][p.iO2]*IdealGas.R*T

    p.eChem.Erev = nernst_potential(H2Oxidation(), p.eChem.Estd, T, aH2=pH2/IdealGas.p_std, aH2O = pH2O/IdealGas.p_std, aO2=pO2/IdealGas.p_std)
    # p.eChem.Erev = nernst_h2(p.eChem.Estd, T, pH2=pH2, pH2O = pH2O, pO2=pO2)
    p.eChem.ecp = p.eChem.udf((pH2=pH2, pH2O = pH2O, pO2=pO2, T = T))        
    solve_electrochemistry(p)    
end


function solve_electrochemistry(p::Cell)
    T = p.anode.ch.gp.T
    σ = (p.electrolyte.ic.sigma/T) * exp(-p.electrolyte.ic.E/IdealGas.R/T)
    Ω = p.electrolyte.t*100/σ # convert to ohm cm2    
    eobject = (cell=p, res=Ω, T = T)

    eprob = NonlinearProblem(interface_kinetics!, p.eChem.vsoln, eobject)
    p.eChem.vsoln = solve(eprob, NewtonRaphson())
    
end

#=
 Function to calculate the electrochemical fluxes at the interface between
 the electrode and the electrolyte 
- eco : ElectrochemObject
- anode : Electrode object
- cathode : Electrode object
=#
function calc_electrochemical_fluxes!(p::SOFC_H2)
    set_interface_flux!(p.anode, p.eChem.vsoln[1], p.iH2, 2)
    set_interface_flux!(p.anode, -p.eChem.vsoln[1], p.iH2O, 2)
    set_interface_flux!(p.cathode, p.eChem.vsoln[1], p.iO2, 4)
end

function calc_electrochemical_fluxes!(p::SOEC_H2)
    set_interface_flux!(p.cathode, -p.eChem.vsoln[1], p.iH2, 2)
    set_interface_flux!(p.cathode, p.eChem.vsoln[1], p.iH2O, 2)
    set_interface_flux!(p.anode, -p.eChem.vsoln[1], p.iO2, 4)
    
end

function electrochemical_species_index(sofc::SOFC_H2)
    sofc.iH2 = get_index("H2", sofc.anode.ch.gp.species)
    sofc.iH2O = get_index("H2O", sofc.anode.ch.gp.species)
    sofc.iO2 = get_index("O2", sofc.cathode.ch.gp.species)      
end

function electrochemical_species_index(soec::SOEC_H2)
    soec.iH2 = get_index("H2", soec.cathode.ch.gp.species)
    soec.iH2O = get_index("H2O", soec.cathode.ch.gp.species)
    soec.iO2 = get_index("O2", soec.anode.ch.gp.species)      
    
end

function electrochemical_species_index(sofc::SOFC_CO)
    sofc.iCO = get_index("CO", sofc.anode.ch.gp.species)
    sofc.iCO2 = get_index("CO2", sofc.anode.ch.gp.species)
    sofc.iO2 = get_index("O2", sofc.cathode.ch.gp.species)                    
end


function init_electrochemitry(sofc::SOFC_H2)    
    T = sofc.anode.ch.gp.T
    # aH2 = sofc.anode.conc[end][sofc.iH2]*IdealGas.R*T/IdealGas.p_std
    # aH2O = sofc.anode.conc[end][sofc.iH2O]*IdealGas.R*T/IdealGas.p_std
    # aO2 = sofc.cathode.conc[end][sofc.iO2]*IdealGas.R*T/IdealGas.p_std
    sofc.eChem.Estd = E0_H2(sofc.ch_anode.gp.thermo_obj, sofc.ch_cathode.gp.thermo_obj, T)    
    sofc.eChem.Erev = 0.0#nernst_potential(H2Oxidation(), sofc.eChem.Estd, T, aH2=aH2, aH2O = aH2O, aO2=aO2)        
    sofc.eChem.vsoln = [0.1, 0.1, -0.1] # initial guess for the solution
end

function init_electrochemitry(soec::SOEC_H2)    
    T = soec.anode.ch.gp.T
    # aH2 = soec.cathode.conc[end][soec.iH2]*IdealGas.R*T/IdealGas.p_std
    # aH2O = soec.cathode.conc[end][soec.iH2O]*IdealGas.R*T/IdealGas.p_std
    # aO2 = soec.anode.conc[end][soec.iO2]*IdealGas.R*T/IdealGas.p_std
    soec.eChem.Estd = E0_H2(soec.ch_anode.gp.thermo_obj, soec.ch_cathode.gp.thermo_obj, T)        
    soec.eChem.Erev = 0.0#nernst_potential(H2Oxidation(), soec.eChem.Estd, T, aH2=aH2, aH2O = aH2O, aO2=aO2)        
    soec.eChem.vsoln = [0.1, 0.1, -0.1] # initial guess for the solution
end


function init_electrochemitry(sofc::SOFC_CO)
    T = sofc.anode.ch.gp.T
    # aCO = sofc.anode.conc[end][sofc.iCO]*IdealGas.R*T/IdealGas.p_std
    # aCO2 = sofc.anode.conc[end][sofc.iCO2]*IdealGas.R*T/IdealGas.p_std
    # aO2 = sofc.cathode.conc[end][sofc.iO2]*IdealGas.R*T/IdealGas.p_std
    sofc.eChem.Estd = E0_CO(sofc.anode.ch.gp.thermo_obj, sofc.anode.ch.gp.T)    
    sofc.eChem.Erev = 0.0#nernst_potential(COxidation(), sofc.eChem.Estd, T, aCO=aCO, aO2=aO2, aCO2=aCO2)   
    sofc.eChem.vsoln = [0.1, 0.1, -0.1] # initial guess for the solution           
end

