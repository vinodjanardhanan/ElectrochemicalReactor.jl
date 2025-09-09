"""
Structure definition of Gasphase 
"""
mutable struct Gasphase 
    species::Array{String,1}
    mole_fractions::Array{Float64,1}
    conc::Array{Float64,1}
    thermo_obj::SpeciesThermoObj
    T::Float64
    p::Float64
end

mutable struct Upstream     
    u::Float64
    ρu::Float64
    mass_fracs::Array{Float64,1}
end
"""
Definition of flow channel
-   gp : Gasphase(speices, mole_fractions, thermo_obj, T, p)
-   u : velocity 
-   h : height of the flow channel 
-   jks : inteface fluxes 
-   m_ptr : pointer to the position of ρu in the soln vector 
-   g_ptr : pointer to the starting point of gasphase concentration in the solution vector 
"""
mutable struct Channel
    gp::Gasphase
    u::Float64    
    h::Float64
    w::Float64    
    jks::Array{Float64}
    upstream::Upstream
    m_ptr::Int64
    g_ptr::UnitRange{Int64}    
    Channel(gp,u,h,w,jk, upstream) = new(gp,u,h,w,jk, upstream, 0, 0:0)
end


"""
Definition of electrode 
-   ch : Channel object
-   t : thickness of the electrode
-   ncells : number of cells along the thickness of the electrode
-   pm : Porous media properties defined in DiffusionFlux.jl
-   mech : mechanism object of present 
-   ws : Work space for flux calculation defined in DiffusionFlux.jl
-   sp_trd : Species transport data 
-   conc : storage space for concentrations 
-   mole_fracs : storage space for mole fractions
-   jks : Fluxes 
-   m_ptr : pointer to the position of ρu  in the solution vector
-   g_ptr : pointer to the position of gasphase concentration in the solution vector 
-   s_ptr : pointer to the position of surface coverage in the solition vector 
"""
mutable struct Electrode
    ch::Channel
    t::Float64
    ncells::Int64    
    pm::Properties
    mech::Union{Nothing,SurfaceMechDefinition}
    srs::Union{Nothing, ReactionState}
    ws::WorkSpace
    sp_trd
    AbyV::Float64
    mole_fracs::Array{Float64,1 }    
    mass_fracs::Array{Float64,1 }
    mass_density::Array{Array{Float64,1},1}
    conc::Array{Array{Float64,1},1}
    p::Array{Float64,1}
    source::Array{Array{Float64,1},1}
    jks::Array{Array{Float64,1},1}    
    g_ptr::Array{UnitRange{Int64}}
    m_ptr::Array{Int64}
    s_ptr::Array{UnitRange{Int64}} 
    Electrode(ch,t,ncells,pm,mech, srs, 
    ws,sp_trd, AbyV ,mole_fracs, mass_fracs, 
    mass_density, conc,p,source,jks) = new(ch,t,ncells,pm,mech,srs, ws,sp_trd, AbyV,mole_fracs,mass_fracs,mass_density,conc,p,source,jks,[],[],[])
end

"""
Arrhenius parameters to determine the conductivity of the ionic conducting phase 
-   sigma : Prefac
-   E : activation energy
"""
struct IonicConductivity 
    sigma::Float64
    E::Float64
end

mutable struct Solver 
    abstol::Float64
    reltol::Float64    
end
Solver() = Solver(1e-6, 1e-6)


"""
Definition of electrolyte 
-   t : electrolyte thickness
-   ic : ionic conductivity
"""
struct Electrolyte  
    t::Float64
    w::Float64    
    ic::IonicConductivity
end


mutable struct EChemParams
    Ecell::Union{Float64,StepRangeLen}
    i0a::Float64
    i0c::Float64
    α_ox_anode::Float64
    α_rd_anode::Float64
    α_ox_cathode::Float64
    α_rd_cathode::Float64
    EChemParams() = new(0.0, 0.0, 0.0, 0.5, 0.5, 0.5, 0.5)    
    EChemParams(Ecell, i0a, i0c, α_ox_anode, α_rd_anode, α_ox_cathode, α_rd_cathode) = 
    new(Ecell, i0a, i0c, α_ox_anode, α_rd_anode, α_ox_cathode, α_rd_cathode)    
end


mutable struct ElectrochemObject
    ecp::EChemParams    
    vsoln::Array{Float64,1}
    Erev::Float64
    Estd::Float64       
    udf::Function    
    ElectrochemObject() = new(EChemParams(), zeros(Float64, 3), 0.0, 0.0, x-> x)
end




abstract type Cell end
abstract type FuelCell <: Cell end
abstract type ElectrolysisCell <: Cell end

mutable struct CellCore 
    ch_anode::Channel
    anode::Electrode
    ch_cathode::Channel
    cathode::Electrode    
    electrolyte::Electrolyte
    δx::Float64
    ncells::Int64
    eChem::ElectrochemObject
    slvr_cntrl::Solver
    CellCore(ch_anode, anode, ch_cathode, cathode, electrolyte, δx, ncells, slvr_cntrl) = 
    new(ch_anode, anode, ch_cathode, cathode, electrolyte, δx, ncells, ElectrochemObject(), slvr_cntrl)               
end

mutable struct SOFC_H2 <: FuelCell
    core::CellCore
    iH2::Int64
    iH2O::Int64
    iO2::Int64
    SOFC_H2(core) = new(core, 0, 0, 0)
end

mutable struct SOEC_H2 <: ElectrolysisCell
    core::CellCore
    iH2::Int64
    iH2O::Int64
    iO2::Int64
    SOEC_H2(core) = new(core, 0, 0, 0)
end

mutable struct HTPEM <: FuelCell
    core::CellCore
    iH2::Int64
    iH2O::Int64
    iO2::Int64    
    HTPEM(core) = new(core, 0, 0, 0)
end

mutable struct SOFC_CO <: FuelCell
    core::CellCore
    iCO::Int64
    iCO2::Int64
    iO2::Int64
    SOFC_CO(core) = new(core, 0, 0, 0)
end
   

struct CellStreams
    anode_channel::IO
    anode::IO
    cathode_channel::IO
    cathode::IO
    echem::IO
end