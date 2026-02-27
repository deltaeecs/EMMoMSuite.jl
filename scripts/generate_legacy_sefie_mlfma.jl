using MoM_AllinOne
using DataFrames, CSV
using MoM_Visualizing

# Set up paths
const MOM_ALLINONE_DIR = "f:/OneDrive/MoM/MoM_AllinOne"
const OUTPUT_DIR = joinpath(@__DIR__, "../test_results/legacy_baseline")

if !isdir(OUTPUT_DIR)
    mkpath(OUTPUT_DIR)
end

function run_legacy_sefie_mlfma()
    println("Running Legacy Case: SEFIE_MLFMA_Jet")
    
    # Reset parameters
    setPrecision!(Float32)
    SimulationParams.SHOWIMAGE = false 

    # Parameters
    filename = joinpath(MOM_ALLINONE_DIR, "meshfiles", "jet_100MHz.nas")
    meshUnit = :m
    frequency = 1e8
    ieT = :EFIE
    sbfT = :RWG
    vbfT = :nothing
    solverT = :gmres
    rtol = 1e-3
    restart = 50
    
    # Source: Legacy default
    # PlaneWave(π/2, 0, 0f0, 1f0) -> Theta=90, Phi=0, Alpha=0, V=1.
    source = PlaneWave(π/2, 0, 0f0, 1f0)

    # Observation
    θs_obs = LinRange{Precision.FT}(-π, π, 721)
    ϕs_obs = LinRange{Precision.FT}(0, π/2, 2)

    # Run Solver
    # We use @eval Main to run the script in global scope as it expects
    result = @eval Main begin
        filename = $filename
        meshUnit = $(QuoteNode(meshUnit))
        frequency = $frequency
        ieT = $(QuoteNode(ieT))
        sbfT = $(QuoteNode(sbfT))
        vbfT = $(QuoteNode(vbfT))
        solverT = $(QuoteNode(solverT))
        rtol = $rtol
        restart = $restart
        source = $source
        θs_obs = $θs_obs
        ϕs_obs = $ϕs_obs
        
        # Include fast solver logic manually to avoid SAI preconditioner bug
        # include(joinpath($MOM_ALLINONE_DIR, "src", "fast_solver.jl"))
        
        using LinearAlgebra

        # Intel MKL 可以带来更好的性能
        # using MKL, MKLSparse

        # 更新参数
        inputParameters(;frequency = frequency, ieT = ieT)
        updateVSBFTParams!(;sbfT = sbfT, vbfT = vbfT)

        # 网格文件读取
        meshData, εᵣs   =  getMeshData(filename; meshUnit=meshUnit);

        # 基函数生成
        ngeo, nbf, geosInfo, bfsInfo =  getBFsFromMeshData(meshData; sbfT = sbfT, vbfT = vbfT)

        # 设置介电参数
        setGeosPermittivity!(geosInfo, 2(1-0.0002im))

        # 八叉树构建
        # octree  =   Octree(geosInfo; nmin = 50)
        nLevels, octree = getOctreeAndReOrderBFs!(geosInfo, bfsInfo; leafCubeEdgel = Precision.FT(0.23Params.λ_0), nInterp = 4)

        # 预计算
        # preCalMLFMAInfo!(octree, geosInfo, bfsInfo)

        # 层数
        # nLevels     =   length(octree.levels)

        # 叶层
        leafLevel   =   octree.levels[nLevels];
        # 计算近场矩阵CSC
        ZnearCSC     =   calZnearCSC(leafLevel, geosInfo, bfsInfo);

        # 构建矩阵向量乘积算子
        Zopt  =   MLFMAIterator(ZnearCSC, octree, geosInfo, bfsInfo);

        ## 根据近场矩阵和八叉树计算 SAI 左预条件
        # Zprel    =   sparseApproximateInversePl(ZnearCSC, leafLevel)
        # Skip SAI due to threading bug
        Zprel = LinearAlgebra.I

        # 激励向量
        V    =   getExcitationVector(geosInfo, nbf, source);

        # 求解
        ICoeff, ch   =   solve(Zopt, V; solverT = solverT, Pl = Zprel, rtol = rtol, restart = restart);

        using LinearAlgebra
        @info norm(ZnearCSC) norm(V)
        # RCS
        radarCrossSection(θs_obs, ϕs_obs, ICoeff, geosInfo) 
    end
    
    # Extract results
    RCS = result[3]
    
    # Save results
    output_file = joinpath(OUTPUT_DIR, "SEFIE_MLFMA_Jet.csv")
    
    df = DataFrame(
        Theta_Rad = θs_obs,
        Theta_Deg = rad2deg.(θs_obs),
        RCS_Phi0_dBsm = 10log10.(RCS[:, 1]),
        RCS_Phi90_dBsm = 10log10.(RCS[:, 2])
    )
    
    CSV.write(output_file, df)
    println("Saved results to $output_file")
end

run_legacy_sefie_mlfma()
