using ElectrochemicalCell
using Test

@testset "ElectrochemicalCell.jl" begin
    if Sys.isapple() || Sys.islinux()
        lib_dir ="lib/"
    elseif Sys.iswindows()
        lib_dir = "lib\\"
    end


    @testset "Testing SOFC" begin
        input_file = joinpath("sofc", "sofc.xml")
        function echemModel(arg)
            # println(arg)
            Ecell = 0.6
            i0a = 10000e-3
            i0c = 10000e-5
            α_ox_anode = 0.5
            α_rd_anode = 0.5
            α_ox_cathode = 0.5
            α_rd_cathode = 0.5 
            EChemParams(Ecell, i0a, i0c, α_ox_anode, α_rd_anode, α_ox_cathode, α_rd_cathode)
        end
        objects, solver_ctrl = get_cell_components(input_file, lib_dir)
        sofc = SOFC_H2(objects.ch_anode, objects.anode, objects.ch_cathode, objects.cathode, objects.electrolyte, objects.δx, objects.ncells)
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
        soec = SOEC_H2(objects.ch_anode, objects.anode, objects.ch_cathode, objects.cathode, objects.electrolyte, objects.δx, objects.ncells)
        run_cell(echemModel, soec, "soec")        
    end

end
