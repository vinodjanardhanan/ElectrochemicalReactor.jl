module ElectrochemicalReactor


using LightXML, RxnHelperUtils, ReactionCommons, Printf
using IdealGas, SurfaceReactions, DiffusionFlux, TransportProperties
using Sundials, DifferentialEquations, SteadyStateDiffEq
using NonlinearSolve, NLsolve

export SOFC_H2, SOEC_H2, SOFC_CO, EChemParams


@enum ElectrochemicalModel ElectrodeReaction=1 PorousElectrode=2


global axpos::Int64




include("Objects.jl")
include("Read.jl")
include("Utils.jl")
include("ElectrochemicalKinetics.jl")

# export EChemParams, get_cell_components
# export electrochemical_species_index, run_cell

export get_cell_components, run_cell


#=
Main function for the calculation of residuals
=#
function residual!(du, u, p, t)
    
    δy_anode = p.anode.t/p.anode.ncells
    δy_cathode = p.cathode.t/p.cathode.ncells
    
    anode_molwts = p.anode.ch.gp.thermo_obj.molwt    
    cathode_molwts = p.cathode.ch.gp.thermo_obj.molwt
    # eco = p[:eco]
    
    
    # Anode update 
    update_half_cell!(p.ch_anode, p.anode, u)
    calc_flux_electrode(p.anode, anode_molwts)

    # Cathode update 
    update_half_cell!(p.ch_cathode, p.cathode, u)
    calc_flux_electrode(p.cathode, cathode_molwts)    

    if electrode_kinetics == ElectrodeReaction        
        electrode_reaction(p)
        calc_electrochemical_fluxes!(p)
    end
    

    resid_channel!(du, u, p.anode.ch, p.electrolyte, p.δx)
    resid_electrode!(du, p.anode, δy_anode)

    resid_channel!(du, u, p.cathode.ch, p.electrolyte, p.δx)
    resid_electrode!(du, p.cathode, δy_cathode)

end


#=
Function to update the properties of the half cell
- ch : Channel object
- eltr : Electrode object
- soln : solution vector
=#
function update_half_cell!(ch::Channel, eltr::Electrode, soln)
    # update the channel properties
    update_properties!(ch, soln)
    # update the electrode properties
    update_properties!(eltr, soln)
end



#=
 Function to calculate the fluxes in the electrode
- electrode : Electrode object
- molwts : molecular weights of the species in the electrode
 =#
function calc_flux_electrode(electrode::Electrode, molwts)
    δy = electrode.t/electrode.ncells
    flux_interface!(electrode.ch.jks, electrode.ch.gp.conc, electrode.conc[1], electrode.ws, molwts, δy)          
    # calculate the mass fluxes within the electrode
    flux_dgm!(electrode.jks, electrode.conc, electrode.pm, electrode.ws, electrode.sp_trd, molwts, electrode.ch.gp.T, δy)        
    set_flux_zeros!(electrode)    
end


#=
function evaluation for species transduction in the channel
- du : residual vector
- u : solution vector 
- ch : Channel object
- ele : Electrolyte object
- δx : width of the finite volume cell
=#
function resid_channel!(du, u, ch::Channel, ele::Electrolyte, δx)    
    
    if axpos == 1
        for l in ch.g_ptr
            i = l - ch.g_ptr.start + 1
            du[l] = ch.upstream.mass_fracs[i] - u[l]
        end
        # println("convection ", ch.upstream.ρu, "\t", u[ch.m_ptr], "\t", δx)
        du[ch.m_ptr] = ch.upstream.ρu - u[ch.m_ptr]        
    else
        ρ = get_density(ch.gp.conc, ch.gp.thermo_obj.molwt)    
        ch.u = u[ch.m_ptr]/ρ
        
        diff_tot = 0.0
        for l in ch.g_ptr
            i = l - ch.g_ptr.start + 1
            convection = (ch.upstream.ρu * ch.upstream.mass_fracs[i] - u[ch.m_ptr]*u[l])/δx
            diff = -ch.jks[i]*ele.w/(ch.h*ch.w) 
            du[l] = (convection+diff)/ρ

            diff_tot += diff
            
        end
        

        # println("convection ", ch.upstream.ρu, "\t", u[ch.m_ptr], "\t", δx)
        convection = (ch.upstream.ρu - u[ch.m_ptr])/δx
        # momentum
        du[ch.m_ptr] =  convection + diff_tot
        
    end
    
end

#=
Function evaluation of species transport in the electrode
# Usage
    resid_electrode!(du, elc, δy)
- du : residual vector
- elc : Electrode object
- δy : width of the finite volume cell
=#
function resid_electrode!(du, elc, δy)
    ng = length(elc.ch.gp.species)    
    for j in 1:elc.ncells
        # sum_flux = 0.0
        if j == 1 # Interface with the channel
            diffusion = -(elc.ch.jks - elc.jks[j])/δy
            source = elc.source[j][1:ng]
            du[elc.g_ptr[j]] = (diffusion + source )/elc.pm.ϵ
            sum_flux = -sum(diffusion + source)            
        else
            diffusion = (elc.jks[j-1] - elc.jks[j])/δy
            source = elc.source[j][1:ng]
            du[elc.g_ptr[j]] = (diffusion + source )/elc.pm.ϵ
            sum_flux = sum(diffusion + source)            
        end
        # Surface coverage integration
        if elc.mech != nothing
            du[elc.s_ptr[j]] = (elc.source[j][ng+1:end] .* elc.mech.sm.si.site_coordination) / (elc.mech.sm.si.density*1e4)
        end
        du[elc.m_ptr[j]] = sum_flux/elc.pm.ϵ
    end
end



#=
Function to read the component properties from the input file
- input_file : path to the input file
- lib_dir : path to the library directory
=#
function get_cell_components(input_file::String, lib_dir::String)
    objects, solver_ctrl = read_input(input_file, lib_dir)    
    return objects, solver_ctrl
end


function run_cell(em::F, cell::Cell, output_file_folder=".") where {F <: Function}
    cell.eChem.udf = em
    soln = prepare_solution(cell.ch_anode, cell.anode, cell.cathode, cell.ch_cathode)
    electrochemical_species_index(cell)
    update_properties!(cell.anode, soln)
    update_properties!(cell.cathode, soln)
    init_electrochemitry(cell)
    cell_streams = create_output_streams(cell, output_file_folder)
    run_cell(cell, soln, cell_streams)
end

function run_cell(object::Cell, soln::Array{Float64,1}, cell_streams::CellStreams) 
    
    ncells = object.ncells
    t_span = (0, 1e5)    
    svrcntrl = object.slvr_cntrl

    println("Axial position in the channel")
    for i in 1:ncells+1                
        global axpos = i
        prob = ODEProblem(residual!, soln, t_span, object)            
        sol = solve(prob, CVODE_BDF(), reltol=svrcntrl.reltol, abstol=svrcntrl.abstol, save_everystep = false)
        soln = sol.u[end]

        reset_upstream(object.ch_anode, soln)
        reset_upstream(object.ch_cathode, soln)
        #axial position
        z = (axpos - 1)* object.δx
        write_output(object, cell_streams, z)

        @printf("%4e\n",z)
    end
    
end


end

