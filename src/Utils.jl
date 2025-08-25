
"""
Function to create solution vector for channel object. 
"""
function create_soln_vector!(soln::Array{Float64}, ch::Channel)
    # mass fractions 
    
    mass_fracs = similar(ch.gp.mole_fractions)
    molefrac_to_massfrac!(mass_fracs, ch.gp.mole_fractions, ch.gp.thermo_obj.molwt)
    sps = length(soln) + 1
    append!(soln, mass_fracs)        
    spe = length(soln)
    ch.g_ptr = sps:spe
    # mass flux (ρu)
    ρ = density(ch.gp.mole_fractions, ch.gp.thermo_obj.molwt, ch.gp.T, ch.gp.p)    
    ch.m_ptr = length(soln)+1    
    push!(soln, ρ*ch.u)
    
end



"""
Function to create solution vector for electrode object. 
# Usage 
create_soln_vector!(soln::Array{Float64}, elc::Electrode; water=false)
-   soln : Array{Float64} solution vector 
-   ele : object of the type electrode  
-   water : A boolean field to specify whether water is formed in this electrode 
"""
function create_soln_vector!(soln::Array{Float64}, eltr::Electrode; water=false)
    Cb = eltr.ch.gp.p/IdealGas.R/eltr.ch.gp.T

    if eltr.mech != nothing
        covg = eltr.mech.sm.si.ini_covg        
    end

    mass_fracs = similar(eltr.mole_fracs)
    ρ = Cb * average_molwt(eltr.mole_fracs, eltr.ch.gp.thermo_obj.molwt)
    molefrac_to_massfrac!(mass_fracs, eltr.mole_fracs, eltr.ch.gp.thermo_obj.molwt)
    for i in 1:eltr.ncells
        # mass density of all species
        sps  =  length(soln)+1
        append!(soln,mass_fracs .* ρ)
        spe = length(soln)        
        push!(eltr.g_ptr, sps:spe)
        # density 
        push!(soln, ρ)
        push!(eltr.m_ptr, length(soln))
        #surface coverage
        if eltr.mech != nothing                    
            cvs = length(soln)+1
            # coverage            
            append!(soln,covg)
            cve = length(soln)                        
            push!(eltr.s_ptr, cvs:cve)            
        end
        
    end
end


#=
Get the mole fractions and calculate the concentrations in the channel
- ch : Channel object
- u : solution vector
=#
function update_properties!(ch::Channel, u)
    # n = length(ch.gp.species)
    # mass fractions 
    ch.gp.mole_fractions =  u[ch.g_ptr]    

    # mole fractions
    massfrac_to_molefrac!(ch.gp.mole_fractions, ch.gp.mole_fractions, ch.gp.thermo_obj.molwt)    

    # Total concentations 
    Cb = ch.gp.p/IdealGas.R/ch.gp.T
    # concentrations 
    ch.gp.conc = ch.gp.mole_fractions .* Cb

    #mass flux 
    ρu = u[ch.m_ptr]
    # println("From update ", u[ch.m_ptr])

    ρ = get_density(ch.gp.conc, ch.gp.thermo_obj.molwt)
    ch.u = ρu/ρ
    # println("Massflux update ", ρu, "\t", ρ, "\t", ch.u)

end


#=
Get the mole fractions in the finite volume cells of the electrode. 
=#
function update_properties!(el::Electrode, u)        

    for j in 1:el.ncells
        # mass density 
        el.mass_density[j] = u[el.g_ptr[j]]   
        ρ = sum(el.mass_density[j]) # total mass density
        # mass fractions
        el.mass_fracs = el.mass_density[j] /ρ
        # mole fractions
        massfrac_to_molefrac!(el.mole_fracs, el.mass_fracs, el.ch.gp.thermo_obj.molwt)
        
        Cb = ρ/average_molwt(el.mole_fracs, el.ch.gp.thermo_obj.molwt)
        p = Cb * IdealGas.R * el.ch.gp.T
        el.conc[j] = el.mole_fracs .* Cb
        if el.mech != nothing
            #update the surface reaction state 
            el.srs.p = p
            el.srs.T = el.ch.gp.T
            el.srs.mole_frac = el.mole_fracs
            el.srs.covg = u[el.s_ptr[j]]            
            SurfaceReactions.calculate_molar_production_rates!(el.srs, el.ch.gp.thermo_obj, el.mech)
            el.source[j] = el.srs.source # these are in mol/m2.s
        else
            el.source[j] = zeros(length(el.ch.gp.species))
        end
    end        
end


#=
Function to reset the upstream properties of the channel
- channel : Channel object
=#
function reset_upstream(channel::Channel, soln)
    channel.upstream.ρu = soln[channel.m_ptr]
    channel.upstream.u = channel.u
    # update the mass fractions
    molefrac_to_massfrac!(channel.upstream.mass_fracs, channel.gp.mole_fractions, channel.gp.thermo_obj.molwt)
end

#=
calculate the density based on concentrations
=#
function get_density(conc::Array{T,1}, molwt::Array{T,1}) where T    
    avg_mol_wt = sum((conc / sum(conc)) .* molwt)
    return sum(conc)*avg_mol_wt
end

#= the last cell face n'th cell flux is set to zero =#
set_flux_zeros!(el::Electrode) = el.jks[el.ncells] = zeros(length(el.ch.gp.species))

#= set the interfce flux based on the current density =#
set_interface_flux!(el::Electrode, cd::Float64, id::Int64, n::Int64) = el.jks[el.ncells][id] = cd*1e4*el.ch.gp.thermo_obj.molwt[id]/(n*IdealGas.F)

function E0_H2(anodeThermo, cathodeThermo, T::Float64)
    molwts = deepcopy(anodeThermo.molwt)
    append!(molwts, cathodeThermo.molwt)
    thermo_all = deepcopy(anodeThermo.thermo_all)
    append!(thermo_all, cathodeThermo.thermo_all)
    echemthermo = SpeciesThermoObj(molwts, thermo_all)
    return IdealGas.E0_H2(echemthermo, T)
end

#= function to create output streams =#
function create_output_streams(object::Cell, output_file_folder::String)
    anode_channel_stream = open(joinpath(output_file_folder, "anode_channel.csv"), "w")
    write_csv(anode_channel_stream, ["z", "T", "u"], object.ch_anode.gp.species)

    cathode_channel_stream = open(joinpath(output_file_folder, "cathode_channel.csv"), "w")
    write_csv(cathode_channel_stream, ["z", "T", "u"], object.ch_cathode.gp.species)

    anode_stream = open(joinpath(output_file_folder, "anode.csv"), "w")
    write_csv(anode_stream, ["z", "y"], object.anode.ch.gp.species)

    cathode_stream = open(joinpath(output_file_folder, "cathode.csv"), "w")
    write_csv(cathode_stream, ["z", "y"], object.cathode.ch.gp.species)

    echem_stream = open(joinpath(output_file_folder, "echem.csv"), "w")
    write_csv(echem_stream, ["z", "Estd", "Erev", "CD", "eta_a", "eta_c"])

    return CellStreams(anode_channel_stream, anode_stream, cathode_channel_stream, cathode_stream, echem_stream)
end

#=Function to close all output stream=#
function close_output_streams(os::CellStreams)
    close(os.anode_channel)
    close(os.cathode_channel)
    close(os.anode)
    close(os.cathode)
    close(os.echem)
end

#= main function call for writing output =#
function write_output(p::Cell, os::CellStreams, z::Float64)
    write_channel(os.anode_channel, z, p.ch_anode)
    write_channel(os.cathode_channel, z, p.ch_cathode)

    write_electrode(os.anode, z, p.anode)
    write_electrode(os.cathode, z, p.cathode)

    write_electrochemistry(os.echem, z, p.eChem)    
    
end

#= write the channel data =#
function write_channel(stream,z,channel)
    write_csv(stream, [z, channel.gp.T, channel.upstream.u], channel.gp.mole_fractions)    
end
#= write the electrode data =#
function write_electrode(stream, z, electrode::Electrode)
    δy = electrode.t/electrode.ncells
    y0 = 0.5*δy
    for i in 1:electrode.ncells
        conc = electrode.conc[i]
        electrode.mole_fracs = conc / sum(conc)
        write_csv(stream, [z, y0 + (i-1)*δy], electrode.mole_fracs)
        
    end    
end

#= write the electrochemistry data =#
function write_electrochemistry(stream, z, eco::ElectrochemObject )
    write_csv(stream, [z, eco.Estd, eco.Erev], eco.vsoln)
end