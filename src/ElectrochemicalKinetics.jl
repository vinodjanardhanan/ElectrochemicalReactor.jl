
function interface_kinetics!(du , u, p)
    
    T = p[:T]
    f = IdealGas.F/IdealGas.R/T
    
    cell = p[:cell]
    eco = cell.core.eChem    
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
    T = p.core.anode.ch.gp.T
    pH2 = p.core.anode.conc[end][p.iH2]*IdealGas.R*T
    pH2O = p.core.anode.conc[end][p.iH2O]*IdealGas.R*T
    pO2 = p.core.cathode.conc[end][p.iO2]*IdealGas.R*T


    p.core.eChem.Erev = nernst_potential(H2Oxidation(), p.core.eChem.Estd, T, aH2=pH2/IdealGas.p_std, aH2O = pH2O/IdealGas.p_std, aO2=pO2/IdealGas.p_std)
    
    # p.core.eChem.Erev = nernst_h2(p.core.eChem.Estd, T, pH2=pH2, pH2O = pH2O, pO2=pO2)
    # The following is a call to a user defined function to calculate the exchange current density
    p.core.eChem.ecp = p.core.eChem.udf((pH2=pH2, pH2O = pH2O, pO2=pO2, T = T))        
    solve_electrochemistry(p)    
end

function electrode_reaction(p::SOEC_H2) 
    T = p.core.anode.ch.gp.T
    pH2 = p.core.cathode.conc[end][p.iH2]*IdealGas.R*T
    pH2O = p.core.cathode.conc[end][p.iH2O]*IdealGas.R*T
    pO2 = p.core.cathode.conc[end][p.iO2]*IdealGas.R*T

    p.core.eChem.Erev = nernst_potential(H2Oxidation(), p.core.eChem.Estd, T, aH2=pH2/IdealGas.p_std, aH2O = pH2O/IdealGas.p_std, aO2=pO2/IdealGas.p_std)
    # p.core.eChem.Erev = nernst_h2(p.core.eChem.Estd, T, pH2=pH2, pH2O = pH2O, pO2=pO2)
    p.core.eChem.ecp = p.core.eChem.udf((pH2=pH2, pH2O = pH2O, pO2=pO2, T = T))        
    solve_electrochemistry(p)    
end

function electrode_reaction(p::HTPEM) 
    T = p.core.anode.ch.gp.T
    pH2 = p.core.anode.conc[end][p.iH2]*IdealGas.R*T
    pH2O = p.core.cathode.conc[end][p.iH2O]*IdealGas.R*T
    pO2 = p.core.cathode.conc[end][p.iO2]*IdealGas.R*T    
    
    p.core.eChem.Erev = nernst_potential(H2Oxidation(), p.core.eChem.Estd, T, aH2=pH2/IdealGas.p_std, aH2O = pH2O/IdealGas.p_std, aO2=pO2/IdealGas.p_std)
    
    p.core.eChem.ecp = p.core.eChem.udf((pH2=pH2, pH2O = pH2O, pO2=pO2, T = T))        
    solve_electrochemistry(p)    
    
end


function solve_electrochemistry(p::Cell)
    T = p.core.anode.ch.gp.T
    σ = (p.core.electrolyte.ic.sigma/T) * exp(-p.core.electrolyte.ic.E/IdealGas.R/T)
    Ω = p.core.electrolyte.t*100/σ # convert to ohm cm2
    eobject = (cell=p, res=Ω, T = T)

    eprob = NonlinearProblem(interface_kinetics!, p.core.eChem.vsoln, eobject)
    p.core.eChem.vsoln = solve(eprob, NewtonRaphson())    

end

#=
 Function to calculate the electrochemical fluxes at the interface between
 the electrode and the electrolyte 
- eco : ElectrochemObject
- anode : Electrode object
- cathode : Electrode object
=#
function calc_electrochemical_fluxes!(p::SOFC_H2)
    set_interface_flux!(p.core.anode, p.core.eChem.vsoln[1], p.iH2, 2)
    set_interface_flux!(p.core.anode, -p.core.eChem.vsoln[1], p.iH2O, 2)
    set_interface_flux!(p.core.cathode, p.core.eChem.vsoln[1], p.iO2, 4)
end

function calc_electrochemical_fluxes!(p::SOEC_H2)
    set_interface_flux!(p.core.cathode, -p.core.eChem.vsoln[1], p.iH2, 2)
    set_interface_flux!(p.core.cathode, p.core.eChem.vsoln[1], p.iH2O, 2)
    set_interface_flux!(p.core.anode, -p.core.eChem.vsoln[1], p.iO2, 4)

end


function calc_electrochemical_fluxes!(p::HTPEM)
    set_interface_flux!(p.core.anode, p.core.eChem.vsoln[1], p.iH2, 2)
    set_interface_flux!(p.core.cathode, -p.core.eChem.vsoln[1], p.iH2O, 2)
    set_interface_flux!(p.core.cathode, p.core.eChem.vsoln[1], p.iO2, 4)    
end

function electrochemical_species_index(sofc::SOFC_H2)
    sofc.iH2 = get_index("H2", sofc.core.anode.ch.gp.species)
    sofc.iH2O = get_index("H2O", sofc.core.anode.ch.gp.species)
    sofc.iO2 = get_index("O2", sofc.core.cathode.ch.gp.species)      
end

function electrochemical_species_index(soec::SOEC_H2)
    soec.iH2 = get_index("H2", soec.core.cathode.ch.gp.species)
    soec.iH2O = get_index("H2O", soec.core.cathode.ch.gp.species)
    soec.iO2 = get_index("O2", soec.core.anode.ch.gp.species)      

end

function electrochemical_species_index(sofc::SOFC_CO)
    sofc.iCO = get_index("CO", sofc.core.anode.ch.gp.species)
    sofc.iCO2 = get_index("CO2", sofc.core.anode.ch.gp.species)
    sofc.iO2 = get_index("O2", sofc.core.cathode.ch.gp.species)                    
end


function electrochemical_species_index(htpem::HTPEM)
    htpem.iH2 = get_index("H2", htpem.core.anode.ch.gp.species)
    htpem.iH2O = get_index("H2O", htpem.core.cathode.ch.gp.species)
    htpem.iO2 = get_index("O2", htpem.core.cathode.ch.gp.species)          
end


function init_electrochemitry(sofc::SOFC_H2)    
    T = sofc.core.anode.ch.gp.T
    sofc.core.eChem.Estd = E0_H2(sofc.core.ch_anode.gp.thermo_obj, sofc.core.ch_cathode.gp.thermo_obj, T)    
    sofc.core.eChem.Erev = 0.0#nernst_potential(H2Oxidation(), sofc.core.eChem.Estd, T, aH2=aH2, aH2O = aH2O, aO2=aO2)        
    sofc.core.eChem.vsoln = [0.1, 0.1, -0.1] # initial guess for the solution
end

function init_electrochemitry(soec::SOEC_H2)    
    T = soec.core.anode.ch.gp.T    
    soec.core.eChem.Estd = E0_H2(soec.core.ch_anode.gp.thermo_obj, soec.core.ch_cathode.gp.thermo_obj, T)        
    soec.core.eChem.Erev = 0.0#nernst_potential(H2Oxidation(), soec.core.eChem.Estd, T, aH2=aH2, aH2O = aH2O, aO2=aO2)        
    soec.core.eChem.vsoln = [0.1, 0.1, -0.1] # initial guess for the solution
end


function init_electrochemitry(sofc::SOFC_CO)
    T = sofc.core.anode.ch.gp.T
    sofc.core.eChem.Estd = E0_CO(sofc.core.anode.ch.gp.thermo_obj, sofc.core.anode.ch.gp.T)    
    sofc.core.eChem.Erev = 0.0#nernst_potential(COxidation(), sofc.core.eChem.Estd, T, aCO=aCO, aO2=aO2, aCO2=aCO2)   
    sofc.core.eChem.vsoln = [0.1, 0.1, -0.1] # initial guess for the solution           
end


function init_electrochemitry(htpem::HTPEM)    
    T = htpem.core.anode.ch.gp.T
    htpem.core.eChem.Estd = E0_H2(htpem.core.ch_anode.gp.thermo_obj, htpem.core.ch_cathode.gp.thermo_obj, T)    
    htpem.core.eChem.Erev = 0.0#nernst_potential(H2Oxidation(), htpem.core.eChem.Estd, T, aH2=aH2, aH2O = aH2O, aO2=aO2)        
    htpem.core.eChem.vsoln = [0.1, 0.1, -0.1] # initial guess for the solution    
end
