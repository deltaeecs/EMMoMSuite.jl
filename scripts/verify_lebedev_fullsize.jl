# 用原版管线（generate_dataset_on_pkpt 的完整 50ρ x 500pos 数据集）验证 pk=13 -> pt=27
# 训练式插值权重在其自身留出集上的表现，并计时。不写 h5 文件。
#
# 用法: julia --project=. scripts/verify_lebedev_fullsize.jl

using EMMoMSuite
using EMMoMSuite.FastAlgorithms.MLFMA.Interpolation: truncation_kernel
using LinearAlgebra, SparseArrays, Statistics, Printf

import EMMoMSuite.FastAlgorithms.Lebedev.LebedevSortedPoints as LSP
import EMMoMSuite.FastAlgorithms.Lebedev.pinv2interpW as PW
import EMMoMSuite.FastAlgorithms.Lebedev.dataset_generator as DG

pk, pt = 13, 27
nInterp = 9
rel_l = EMMoMSuite.Utilities.find_zero_bisection(x -> truncation_kernel(x) - (pk + 1) / 2, 0)

t0 = time()
T, P = DG.generate_dataset_on_pkpt(pk, pt, rel_l)
t1 = time()
@printf("数据集生成(50ρ x 500pos): %.1f s, tArray %s, pArray %s\n",
    t1 - t0, size(T), size(P))

flag = trunc(Int, 0.8 * size(T, 2))
xx2D = vcat(real(T[:, 1:flag]), imag(T[:, 1:flag]))
yy2D = vcat(real(P[:, 1:flag]), imag(P[:, 1:flag]))
Tte, Pte = T[:, (flag+1):end], P[:, (flag+1):end]

τt, τp = (pk - 1) ÷ 2, (pt - 1) ÷ 2
tnodes, pnodes = LSP.get_t_nodes(τt), LSP.get_t_nodes(τp)
w = PW.interpWeightsInitial(tnodes, pnodes; nInterp = nInterp)
t2 = time()
PW.pinv2W!(w, nInterp, xx2D, yy2D)
t3 = time()
@printf("逐行 pinv 拟合: %.1f s, W 非零元=%d (每行 %.1f)\n",
    t3 - t2, nnz(w), nnz(w) / size(w, 1))

ε = mean(PW.acc(w, Tte, Pte))
@printf("完整数据集下 现行管线(实部约束) 自身留出集平均相对误差: %.3e\n", ε)
@printf("（saveInterpW2file 初始 fileε=1.0，无绝对精度阈值：ε<1.0 的烂矩阵第一次也会被保存）\n")
