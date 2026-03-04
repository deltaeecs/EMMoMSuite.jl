"""
Phase 20 端口体系测试

TDD RED 阶段：测试先于实现，初次运行时预期全部失败。

涵盖：
  20.2 LumpedPort
  20.3 DiscretePort / CoaxPort
  20.4 DifferentialPairPort
  20.5 SParameterData + Touchstone I/O
  20.1 WavePort（矩形波导 TE10 解析验证）
"""

using Test
using LinearAlgebra
using EMSuite   # 通过主模块导入 Ports 相关类型

# ─────────────────────────────────────────────────────────────────────────────
# 辅助函数
# ─────────────────────────────────────────────────────────────────────────────

"""生成一个简单的 3 节点线段网格（用于 LumpedPort 单元测试）"""
function _simple_wire_mesh()
    # 3 节点，2 边的线段（近似矩形导线）
    # 这里借用 TriangleMesh 的两个退化三角形
    # 实际测试只需检查矩阵维度与对角添加正确性
    return nothing  # placeholder
end

# ─────────────────────────────────────────────────────────────────────────────
# 20.2 LumpedPort
# ─────────────────────────────────────────────────────────────────────────────

@testset "LumpedPort" begin

    @testset "构造函数" begin
        # 基本构造：电压源，50 Ω 端口阻抗，位于 edge_idx=5
        port = LumpedPort(1, 5, 50.0 + 0.0im, :voltage_source, 1.0 + 0.0im, 0.0 + 0.0im)
        @test port.id == 1
        @test port.edge_idx == 5
        @test port.impedance ≈ 50.0 + 0.0im
        @test port.port_type == :voltage_source
    end

    @testset "工厂函数" begin
        # lumped_port 简便构造
        port_v = lumped_port(1, 5; impedance=50.0, type=:voltage_source, voltage=1.0)
        @test port_v.id == 1
        @test port_v.edge_idx == 5
        @test real(port_v.impedance) ≈ 50.0
        @test port_v.port_type == :voltage_source

        port_l = lumped_port(2, 7; impedance=50.0, type=:load)
        @test port_l.port_type == :load
    end

    @testset "阻抗矩阵贡献" begin
        # 4×4 单位矩阵，在 edge_idx=2 添加 100+50j Ω 负载
        Z = Matrix{ComplexF64}(I, 4, 4)
        port = lumped_port(1, 2; impedance=100.0 + 50.0im, type=:load)
        assemble_lumped_port_impedance!(Z, port)
        # 对角项应增加 port.impedance
        @test Z[2, 2] ≈ 1.0 + (100.0 + 50.0im)
        # 其他对角项不变
        @test Z[1, 1] ≈ 1.0 + 0.0im
        @test Z[3, 3] ≈ 1.0 + 0.0im
        # 非对角项不变
        @test Z[1, 2] ≈ 0.0 + 0.0im
    end

    @testset "激励向量贡献" begin
        # 电压源端口在 edge_idx=3，V_src=2.0V，于 4 元素向量中
        V = zeros(ComplexF64, 4)
        port = lumped_port(1, 3; impedance=50.0, type=:voltage_source, voltage=2.0 + 0.0im)
        add_port_excitation!(V, port)
        @test V[3] ≈ 2.0 + 0.0im   # 仅该边被激励
        @test V[1] ≈ 0.0
        @test V[2] ≈ 0.0
    end

    @testset "端口电压/电流/功率提取" begin
        # I = [0, 0.01+0j, 0, 0]，port at edge_idx=2，Z=50Ω
        I_coeff = ComplexF64[0.0, 0.01, 0.0, 0.0]
        port    = lumped_port(1, 2; impedance=50.0, type=:load)
        V_p = port_voltage(port, I_coeff)
        I_p = port_current(port, I_coeff)
        P_p = port_power(port, I_coeff)

        @test I_p ≈ 0.01 + 0.0im
        @test V_p ≈ 0.01 * 50.0 + 0.0im        # V = Z * I
        @test P_p ≈ 0.5 * real(V_p * conj(I_p)) # P = 0.5 Re(V I*)
        @test P_p ≈ 0.5 * (0.01)^2 * 50.0
    end

    @testset "匹配负载 S11 = 0" begin
        # 单端口匹配条件：Z_efie = Z_ref = 50Ω → (Z-Z0)/(Z+Z0) = 0
        Z_efie = ComplexF64[50.0][:, :]   # 1×1，端口阻抗 = Z_ref
        port   = lumped_port(1, 1; impedance=50.0, type=:load)
        s_mat  = extract_s_matrix([port], Z_efie; Z_ref=50.0)
        @test abs(s_mat[1, 1]) < 1e-10
    end

    @testset "单端口 S11 = (Z-Z0)/(Z+Z0)" begin
        # Z_ant = 100Ω, Z_ref = 50Ω → S11 = 50/150 = 1/3
        Z_efie = ComplexF64[100.0][:, :]
        port   = lumped_port(1, 1; impedance=50.0, type=:load)
        s_mat  = extract_s_matrix([port], Z_efie; Z_ref=50.0)
        @test s_mat[1, 1] ≈ (100.0 - 50.0) / (100.0 + 50.0)
    end

    @testset "2端口 S 矩阵对称性（互易网络）" begin
        # 对称 2×2 Z 矩阵 → 对称 S 矩阵
        Z_efie = ComplexF64[100.0 10.0; 10.0 80.0]
        p1 = lumped_port(1, 1; impedance=50.0, type=:load)
        p2 = lumped_port(2, 2; impedance=50.0, type=:load)
        s_mat = extract_s_matrix([p1, p2], Z_efie; Z_ref=50.0)
        @test s_mat[1, 2] ≈ s_mat[2, 1]   atol=1e-12  # 互易性
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# 20.3 DiscretePort + CoaxPort
# ─────────────────────────────────────────────────────────────────────────────

@testset "DiscretePort" begin
    @testset "构造函数" begin
        port = DiscretePort(1, 3, 5, :R, 50.0 + 0.0im)
        @test port.id == 1
        @test port.node_positive == 3
        @test port.node_negative == 5
        @test port.element_type == :R
        @test real(port.value) ≈ 50.0
    end
end

@testset "CoaxPort" begin
    @testset "Z_c 解析公式" begin
        # 标准 50Ω 同轴：b/a = e^(2π*50/η₀) ≈ e^(2π*50/376.73) ≈ e^0.8327 ≈ 2.3
        # 反过来：给定 a=1mm, b=2.3mm, ε_r=1
        a = 1e-3; b = exp(2π * 50.0 / 376.730313) * a   # ≈ 2.3 mm
        port = CoaxPort(1, a, b, 1.0)
        Zc   = coax_impedance(port)
        @test abs(Zc - 50.0) < 0.01   # < 0.01Ω 误差（< 0.02%）
    end

    @testset "有填充介质" begin
        # ε_r=4 时 Z_c 减半
        a = 1e-3; b = exp(2π * 50.0 / 376.730313) * a
        port2 = CoaxPort(1, a, b, 4.0)
        Zc2   = coax_impedance(port2)
        @test abs(Zc2 - 25.0) < 0.01
    end

    @testset "短路时 Z_c → 0（b→a 极限）" begin
        port3 = CoaxPort(1, 1.0e-3, 1.0e-3 * (1 + 1e-9), 1.0)
        Zc3   = coax_impedance(port3)
        @test abs(Zc3) < 0.001  # ln(b/a) ≈ 1e-9 → Z → 0
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# 20.4 DifferentialPairPort
# ─────────────────────────────────────────────────────────────────────────────

@testset "DifferentialPairPort" begin
    @testset "构造函数" begin
        pair = DifferentialPairPort(1, 1, 2)
        @test pair.id == 1
        @test pair.port_positive == 1
        @test pair.port_negative == 2
    end

    @testset "变换矩阵酉性 M·M† = I" begin
        # 2 端口差分对（1 对），M 是 2×2 矩阵
        M2 = mixed_mode_transform_matrix([DifferentialPairPort(1, 1, 2)])
        @test M2 * M2' ≈ I   atol=1e-12
    end

    @testset "4 端口变换矩阵酉性" begin
        pairs = [DifferentialPairPort(1, 1, 2), DifferentialPairPort(2, 3, 4)]
        M4 = mixed_mode_transform_matrix(pairs)
        @test size(M4) == (4, 4)
        @test M4 * M4' ≈ I   atol=1e-12
    end

    @testset "对称传输线 SCD = 0" begin
        # 完全对称 2 端口 → S_cc = S_dd 对角；S_cd = 0
        # 构造对称 single-ended S 矩阵
        S_se = ComplexF64[0.0 1.0; 1.0 0.0]   # 理想传输（S21=S12=1, S11=S22=0）
        pairs = [DifferentialPairPort(1, 1, 2)]
        S_mm = convert_to_mixed_mode(S_se, pairs)
        # 对称情况下 S_cd = S_mm[1,2] 和 S_mm[2,1] 应全为 0
        @test abs(S_mm[1, 2]) < 1e-12
        @test abs(S_mm[2, 1]) < 1e-12
    end

    @testset "sdd/scc/scd 提取一致性" begin
        pairs = [DifferentialPairPort(1, 1, 2)]
        S_se  = ComplexF64[0.1 0.9; 0.9 0.1]
        S_mm  = convert_to_mixed_mode(S_se, pairs)
        @test sdd_matrix(S_mm, pairs) isa Matrix
        @test scc_matrix(S_mm, pairs) isa Matrix
        @test scd_matrix(S_mm, pairs) isa Matrix
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# 20.5 SParameterData + Touchstone I/O
# ─────────────────────────────────────────────────────────────────────────────

@testset "SParameterData" begin
    @testset "构造函数" begin
        freqs = [1e9, 2e9, 3e9]
        S3 = [ComplexF64[0.1 0.9; 0.9 0.1] for _ in freqs]
        data = SParameterData(freqs, S3, [50.0, 50.0], ["test comment"])
        @test length(data.frequencies) == 3
        @test size(data.s_matrices[1]) == (2, 2)
        @test data.port_impedances == [50.0, 50.0]
    end

    @testset "Touchstone s1p round-trip" begin
        mktempdir() do dir
            path = joinpath(dir, "test.s1p")
            freqs = [1e9, 2e9]
            S = [ComplexF64[-0.1 + 0.05im;;] for _ in freqs]
            data_out = SParameterData(freqs, S, [50.0], String[])
            write_touchstone(data_out, path)
            @test isfile(path)
            data_in = read_touchstone(path)
            @test length(data_in.frequencies) == 2
            @test data_in.frequencies[1] ≈ 1e9
            @test data_in.s_matrices[1][1, 1] ≈ -0.1 + 0.05im   atol=1e-6
        end
    end

    @testset "Touchstone s2p round-trip" begin
        mktempdir() do dir
            path = joinpath(dir, "test.s2p")
            freqs = [1e9, 10e9]
            S = [ComplexF64[0.0 1.0; 1.0 0.0] for _ in freqs]
            data_out = SParameterData(freqs, S, [50.0, 50.0], ["# 2-port ideal transmission"])
            write_touchstone(data_out, path)
            data_in = read_touchstone(path)
            @test data_in.s_matrices[1] ≈ S[1]   atol=1e-10
            @test data_in.s_matrices[2] ≈ S[2]   atol=1e-10
        end
    end

    @testset "多频点 round-trip 精度" begin
        mktempdir() do dir
            path = joinpath(dir, "test_multi.s2p")
            freqs = collect(range(1e9, 10e9, 5))
            # 随机复数 S 矩阵（归一化）
            rng_S = [ComplexF64[0.1+0.2im 0.8-0.1im;
                                0.8-0.1im 0.1+0.2im] for _ in freqs]
            data_out = SParameterData(freqs, rng_S, [50.0, 50.0], String[])
            write_touchstone(data_out, path)
            data_in = read_touchstone(path)
            for i in 1:5
                @test data_in.s_matrices[i] ≈ rng_S[i]   atol=1e-10
            end
        end
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# 20.1 WavePort — 矩形波导 TE10 模式解析验证
# ─────────────────────────────────────────────────────────────────────────────

@testset "WavePort" begin

    @testset "构造函数" begin
        port = WavePort(1; mode=:TE10, a=0.02286, b=0.01016,
                        de_embed_length=0.0, reference_impedance=50.0)
        @test port.id == 1
        @test port.mode == :TE10
        @test port.a ≈ 0.02286
    end

    @testset "矩形波导 TE10 β 误差 < 0.1%" begin
        # X 波段波导 WR-90：a=22.86mm, b=10.16mm
        a = 22.86e-3; b = 10.16e-3
        freq = 10e9
        c0   = 299792458.0
        k    = 2π * freq / c0
        # 解析 β_TE10 = √(k² − (π/a)²)
        beta_analytic = sqrt(k^2 - (π/a)^2)

        port  = WavePort(1; mode=:TE10, a=a, b=b, de_embed_length=0.0)
        modes = compute_port_modes(port, freq)
        beta_computed = modes.beta[1]

        rel_err = abs(beta_computed - beta_analytic) / abs(beta_analytic)
        @test rel_err < 1e-3   # < 0.1%
    end

    @testset "截止频率以下 β 为纯虚数" begin
        a = 22.86e-3; b = 10.16e-3
        freq_cut = 299792458.0 / (2 * a)   # f_c = c0 / (2a) ≈ 6.556 GHz
        freq     = 0.9 * freq_cut           # 低于截止
        port     = WavePort(1; mode=:TE10, a=a, b=b, de_embed_length=0.0)
        modes    = compute_port_modes(port, freq)
        # β 应为纯虚
        @test abs(real(modes.beta[1])) < 1e-6
        @test abs(imag(modes.beta[1])) > 1.0
    end

    @testset "去嵌入相位差 = 2βL" begin
        a = 22.86e-3; b = 10.16e-3; freq = 10e9
        c0 = 299792458.0; k = 2π*freq/c0
        beta = sqrt(complex(k^2 - (π/a)^2))
        L    = 0.05   # 5 cm
        # 构造 S11=0, S21=exp(-j*beta*L*2) 的单端口（transmission line）
        S_raw = ComplexF64[0.0 exp(-im*beta*L); exp(-im*beta*L) 0.0]
        port1 = WavePort(1; mode=:TE10, a=a, b=b, de_embed_length=L)
        port2 = WavePort(2; mode=:TE10, a=a, b=b, de_embed_length=L)
        S_deembed = de_embed_s_matrix(S_raw, [port1, port2], freq)
        # 去嵌入后 S21 ≈ 1 (相位完全补偿)
        @test abs(S_deembed[1, 2]) ≈ 1.0   atol=1e-10
        @test abs(S_deembed[2, 1]) ≈ 1.0   atol=1e-10
        @test abs(S_deembed[1, 1]) < 1e-10
    end
end
