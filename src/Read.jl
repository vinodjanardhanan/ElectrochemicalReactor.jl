
#=
This function reads the input files and create the data structures required for the electrochemical cell simulation.    
# Usage 
    read_input(input_file::String, lib_dir::String) 
-   input_file : path to the input file
-   lib_dir : path to the library directory
=#

function read_input(input_file::String, lib_dir::String) where F

    thermo_file = get_path(lib_dir, "therm.dat")
        
    xmldoc = parse_file(input_file)
    xmlroot = root(xmldoc)

    

    
    
    global electrode_kinetics = get_electrode_kinetics(xmlroot)
    
    
    cell_length = get_value_from_xml(xmlroot, "length")
    cellT = get_value_from_xml(xmlroot, "cellT")
    ncells = Integer(get_value_from_xml(xmlroot, "ncells"))



    # read the anode channel
    ac = xmlroot["anode_channel"]
    ch_anode = get_flow_channel(ac[1], cellT, thermo_file)
    # read the cathode channel
    cc = xmlroot["cathode_channel"]
    ch_cathode = get_flow_channel(cc[1], cellT, thermo_file)

    # read the anode electrode 
    electrode = xmlroot["anode"]
    anode = get_electrode(electrode[1], ch_anode, lib_dir)

    # read the cathode electrode 
    electrode = xmlroot["cathode"]
    cathode = get_electrode(electrode[1], ch_cathode, lib_dir)

    # read the electrolyte properties 
    elec = xmlroot["electrolyte"]
    electrolyte = get_electrolyte(elec[1])

    solver_ctrl = get_solver_controls(xmlroot)
    δx = cell_length/ncells
    
    
        
    objects = (ch_anode = ch_anode, # anode channel
        anode=anode, # anode 
        cathode=cathode, # cathode 
        ch_cathode=ch_cathode, # cathode channel
        electrolyte=electrolyte, # electrolyte 
        δx=δx, # axial discretization
        ncells=ncells # number of axial cells         
    )
    return objects, solver_ctrl

end

#=
Function to prepare the solution vector for the electrochemical cell simulation.
    The solution vector is arranged in the fashion, channel followed by the elctrode,
    which forms a half cell. In this arrangement the first cell in the electrode
    is the interface cell with the channel. 

# Usage 
    prepare_solution(ch_anode, anode, cathode, ch_cathode)
-   ch_anode : Channel object for the anode
-   anode : Electrode object for the anode
-   cathode : Electrode object for the cathode
-   ch_cathode : Channel object for the cathode
=#
function prepare_solution(ch_anode, anode, cathode, ch_cathode)
    # create the solution vector
    soln = Array{Float64,1}()   
    # anode half cell
    prepare_soln_half_cell!(soln, ch_anode, anode)
    # cathode half cell
    prepare_soln_half_cell!(soln, ch_cathode, cathode)
    
    return soln
    
end

#=
Function to prepare the solution vector for the half cell.
- soln : solution vector
- ch : Channel object
- eltr : Electrode object
=#
function prepare_soln_half_cell!(soln, ch, eltr)
    create_soln_vector!(soln, ch)
    create_soln_vector!(soln, eltr)
end


    


function get_electrode_kinetics(xmlroot::XMLElement)
    rxn_model = attribute(xmlroot, "model")
    if rxn_model == "ElectrodeReaction"
        return ElectrodeReaction
    elseif rxn_model == "PorousElectrode"
        return PorousElectrode
    else
        error("kinetics model not specified in the input file")
    end
end

"""
Function to read the channel data
# Usage 
    get_flow_channel(channel::XMLElement, T, l, thermo_file)
-   channel : XMLElement
-   T : Temperature (K)
-   thermo_file : path to therm.dat file 
"""
function get_flow_channel(channel::XMLElement, T, thermo_file)
    channel_species = get_collection_from_xml(channel, "gasphase")
    thermo_obj = create_thermo(channel_species, thermo_file)
    mole_fracs = get_molefraction_from_xml(channel, thermo_obj.molwt, channel_species)
    u = get_value_from_xml(channel, "u")
    p = get_value_from_xml(channel, "p")
    h = get_value_from_xml(channel, "height")
    w = get_value_from_xml(channel, "width")
    conc = mole_fracs*(p/IdealGas.R/T)
    gp = Gasphase(channel_species, mole_fracs, conc, thermo_obj, T, p)        
    jks = Array{Float64,1}(undef, length(channel_species))
    mass_fracs = similar(mole_fracs)
    molefrac_to_massfrac!(mass_fracs, mole_fracs, thermo_obj.molwt)
    ρ = density(mole_fracs, thermo_obj.molwt, T, p)
    upstream = Upstream(u, ρ*u, mass_fracs)
    return Channel(gp, u, h, w, jks, upstream)
end


"""
Function to read the electrode data 
"""
function get_electrode(electrode::XMLElement, channel::Channel, lib_dir)
    t = get_value_from_xml(electrode, "thickness")    
    n_cells = Integer(get_value_from_xml(electrode, "ncells"))
    ϵ = get_value_from_xml(electrode, "porosity")
    τ = get_value_from_xml(electrode, "tortuosity")        
    pore_dia = get_value_from_xml(electrode, "pore_dia")
    part_dia = get_value_from_xml(electrode, "particle_dia")
    AbyV = get_value_from_xml(electrode, "AsV") # default value is 1e3
    props = Properties(ϵ, τ, pore_dia, part_dia)
    mole_fracs = get_molefraction_from_xml(electrode, channel.gp.thermo_obj.molwt, channel.gp.species)    
    mass_fracs = similar(mole_fracs)
    mech = get_text_from_xml(electrode, "surf_mech")
    if mech != nothing
        file_path = get_path(lib_dir, mech)                
        surf_mech = SurfaceReactions.compile_mech(file_path, channel.gp.thermo_obj, channel.gp.species)
        covg = surf_mech.sm.si.ini_covg        
        surf_conc = similar(covg)
        surf_rxn_rate = zeros(length(surf_mech.sm.reactions))
        nspecies = length(channel.gp.species) + length(covg)
        source = zeros(nspecies)
        conc_all = zeros(nspecies)
        surf_state = SurfaceRxnState(channel.gp.T, channel.gp.p, mole_fracs, covg, surf_conc, surf_rxn_rate, source, conc_all)
    else
        surf_mech = nothing
        surf_state = nothing
    end

    # create work space 
    n_sp = length(channel.gp.species)
    Dkn_e = Array{Float64, 1}(undef, n_sp)
    D_ij_e = Matrix(undef, n_sp, n_sp)
    Dkl_DGM = Matrix(undef, n_sp, n_sp)
    Dkm = Array{Float64,1}()
    ws = WorkSpace(Dkn_e, D_ij_e, Dkl_DGM, Dkm)
    tr_file = get_path(lib_dir, "transport.dat")
    sp_trd = create_transport_data(channel.gp.species, tr_file)
    jks = Array{Array{Float64,1},1}(undef, n_cells)    
    ρY = Array{Array{Float64,1},1}(undef, n_cells)
    conc = Array{Array{Float64,1},1}(undef, n_cells)
    source = Array{Array{Float64,1},1}(undef, n_cells)
    
    
    effective_coefficients!(ws, props, sp_trd, channel.gp.p, channel.gp.T, channel.gp.thermo_obj.molwt)
    return Electrode(channel, t, n_cells, props, surf_mech, surf_state, ws, sp_trd, AbyV, mole_fracs, mass_fracs, ρY, conc, source, jks)
end


"""
Function to read the electrolyte data
"""
function get_electrolyte(elec::XMLElement)
    t = get_value_from_xml(elec, "thickness")
    w = get_value_from_xml(elec, "width")
    sigma = get_value_from_xml(elec, "sigma")
    E = get_value_from_xml(elec, "E")
    ic = IonicConductivity(sigma, E)
    return Electrolyte(t, w, ic)
end


"""
Function to get the solver controls
"""
function get_solver_controls(xmlroot::XMLElement)
    solver = xmlroot["solver"]
    abstol = get_value_from_xml(solver[1], "abstol")
    reltol = get_value_from_xml(solver[1], "reltol")
    return Solver(abstol, reltol)
end

