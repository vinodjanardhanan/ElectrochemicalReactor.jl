using ElectrochemicalReactor, IdealGas 
using Test

@testset "ElectrochemicalReactor.jl" begin
    if Sys.isapple() || Sys.islinux()
        lib_dir ="lib/"
    elseif Sys.iswindows()
        lib_dir = "lib\\"
    end


    @testset "Testing SOFC" begin
        input_file = joinpath("sofc", "sofc.xml")
        function i0a(p)
            Ades = 5.59e19 # s cm2/mol 
            Edes = 88.12e3 # J/mol
            Γ = 2.6e-9 #mol/cm2
            γ0 = 0.01
            i0 = 8.5 # A/cm2
            pre_fac =  (Ades/10000)*(Γ*10000)^2
            pstarH2 = pre_fac *  sqrt(2*π*IdealGas.R*p.T*0.002) * exp(-Edes/(IdealGas.R*p.T))/γ0
            Nr = (p.pH2/pstarH2)^0.25 * (p.pH2O/IdealGas.p_std)^0.75
            Dr = 1 + (p.pH2/pstarH2)^0.5
            return i0 * Nr/Dr            
        end

        function i0c(p)
            Ao2 = 4.9e8 # atm 
            Eo2 = 200e3     # J/mol
            pstar_o2 = IdealGas.p_std * Ao2 * exp(-Eo2/(IdealGas.R*p.T))
            Nr = (p.pO2/pstar_o2)^0.25 
            Dr = 1 + (p.pO2/pstar_o2)^0.5
            return 2.8 * Nr/Dr                 
        end

        function echemModel(p)            
            Ecell = 0.5            
            α_ox_anode = 1+0.5
            α_rd_anode = 0.5
            i0anode = i0a(p)
            i0cathode = i0c(p)            
            α_ox_cathode = 0.5
            α_rd_cathode = 0.5 
            return EChemParams(Ecell, i0anode, i0cathode, α_ox_anode, α_rd_anode, α_ox_cathode, α_rd_cathode)
        end
        objects, solver_ctrl = get_cell_components(input_file, lib_dir)
        sofc_core = CellCore(objects.ch_anode, objects.anode, objects.ch_cathode, objects.cathode, objects.electrolyte, objects.δx, objects.ncells, solver_ctrl)
        sofc = SOFC_H2(sofc_core)        
        run_cell(echemModel, sofc, "sofc")        
    end


    @testset "Testing SOEC" begin
        input_file = joinpath("soec", "soec.xml")
        function echemModel(arg)
            # println(arg)
            Ecell = 1.5
            i0a = 10000e-5
            i0c = 10000e-3
            α_ox_anode = 0.5
            α_rd_anode = 0.5
            α_ox_cathode = 0.5
            α_rd_cathode = 0.5 
            EChemParams(Ecell, i0a, i0c, α_ox_anode, α_rd_anode, α_ox_cathode, α_rd_cathode)
        end
        objects, solver_ctrl = get_cell_components(input_file, lib_dir)
        soec_core  = CellCore(objects.ch_anode, objects.anode, objects.ch_cathode, objects.cathode, objects.electrolyte, objects.δx, objects.ncells, solver_ctrl)
        soec = SOEC_H2(soec_core)
        run_cell(echemModel, soec, "soec")        
    end



    @testset "Testing HTPEM" begin
        input_file = joinpath("htpem", "htpem.xml")
        function echemModel(arg)
            # println(arg)
            Ecell = 0.6
            i0a = 10000e-5
            i0c = 10000e-8
            α_ox_anode = 0.5
            α_rd_anode = 0.5
            α_ox_cathode = 0.5
            α_rd_cathode = 0.5 
            EChemParams(Ecell, i0a, i0c, α_ox_anode, α_rd_anode, α_ox_cathode, α_rd_cathode)
        end
        objects, solver_ctrl = get_cell_components(input_file, lib_dir)
        htpem_core  = CellCore(objects.ch_anode, objects.anode, objects.ch_cathode, objects.cathode, objects.electrolyte, objects.δx, objects.ncells, solver_ctrl)
        htpem = HTPEM(htpem_core)
        run_cell(echemModel, htpem, "htpem")        
    end

end
