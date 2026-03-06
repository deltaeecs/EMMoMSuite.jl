# EMSuite 閲嶆瀯璺嚎鍥?

> 鏈€鍚庢洿鏂? 2026-03-XX (MLFMA 算法分析阶段完成)

## MLFMA 算法分析阶段 ✅

- [x] 代码清理 (scripts 80+, benchmark 50+, root waste files)
- [x] Legacy MLFMA 完整算法报告 → `.github/plans/legacy_mlfma_algorithm_report.md`
- [x] EMSuite vs Legacy 系数链验证 (总系数一致: k²η/16π²)
- [x] 问题定位: nLevels≥3 错误源于 upward/downward pass (插值/相移)
- [x] 2轮检视完成

## 椤圭洰姒傝堪

灏嗗垎鏁ｇ殑 Legacy 浠ｇ爜 (`MoM_Basics`, `MoM_Kernels`, `MoM_MPI`, `MoM_Lebedev`, `MoM_AllinOne`, `MPIArray4MoMs`) 閲嶆瀯涓虹粺涓€鐨?`EMSuite.jl` 鍖呫€?

## 鏋舵瀯

```
EMSuite.jl/src/
鈹溾攢鈹€ Core/            # 鎶借薄鎺ュ彛銆佺被鍨嬨€佸父鏁般€佹潗鏂欍€佹縺鍔辨簮銆侀厤缃?
鈹溾攢鈹€ Geometry/        # 缃戞牸绫诲瀷銆両/O (Nastran/Gmsh)銆佸潗鏍囧彉鎹€侀珮鏂眰绉?
鈹溾攢鈹€ BasisFunctions/  # RWG, SWG, RBF, PWC
鈹溾攢鈹€ IntegralEquations/ # EFIE, MFIE, CFIE, VEFIE, SCFIE, 濂囧紓鎬у鐞?
鈹溾攢鈹€ FastAlgorithms/  # MLFMA (Octree, Agg, Trans, Disagg), Lebedev
鈹溾攢鈹€ Solvers/         # Direct (LU), Iterative (GMRES, BiCGSTAB), 棰勬潯浠跺櫒
鈹溾攢鈹€ PostProcessing/  # RCS, FarField, NearField, 鐢垫祦鍒嗗竷
鈹溾攢鈹€ Parallel/        # MPI, Threading
鈹溾攢鈹€ IO/              # VTK, HDF5, CSV, 缁撴灉鏂囦欢
鈹溾攢鈹€ Utilities/       # 鏃ュ織銆佸弬鏁般€丮ie 绾ф暟
鈹斺攢鈹€ Driver.jl        # 缁熶竴浠跨湡鍏ュ彛
```

---

## 闃舵鍒掑垎

### Phase 1鈥?: 鍩虹鏋舵瀯 鈫?鍚庡鐞?鉁?

鍏ㄩ儴瀹屾垚銆傚寘鍚? 椤圭洰缁撴瀯銆丆ore 妯″潡銆佸嚑浣曚笌鍩哄嚱鏁?(RWG/SWG/RBF/PWC)銆佺Н鍒嗘柟绋?(EFIE/MFIE/CFIE/VEFIE/SCFIE + 濂囧紓鎬у鐞?銆丮LFMA (鍏弶鏍?鑱氬悎/杞Щ/瑙ｈ仛/Lebedev)銆佹眰瑙ｅ櫒 (Direct/GMRES/BiCGSTAB + ILU/SPAI 棰勬潯浠?銆丮PI/Threading 骞惰銆佸悗澶勭悊 (RCS/FarField/NearField/VTK)銆?

### Phase 7: 楠岃瘉涓庡榻?鉁?

- [x] SCFIE MLFMA 杩戝満/杩滃満楠岃瘉 鈥?10 涓?Bug 淇
- [x] MoM_AllinOne 鍏ㄩ儴绠椾緥瀵规爣瀹屾垚
- [x] 鍥炲綊娴嬭瘯閿佸畾: 138/138 鍏ㄩ儴閫氳繃
- [x] ~3 dB 绯荤粺鍋忓樊鏍瑰洜淇: edgev虃 鏂瑰悜 + near-interaction 闈㈢Н褰掍竴鍖?

### Phase 8: 鎬ц兘浼樺寲 鉁?

> 璇﹁ `test_results/PERFORMANCE_REPORT.md`
> - [x] 8.0 鎬ц兘鍩虹嚎娴嬮噺 鉁?(`861426d`)
> - [x] 8.1 Z 缁勮澶氱嚎绋嬪幓閿?鉁?(`2d4ebe6`) 鈥?Plate EFIE -54%, Jet EFIE -12%
> - [x] 8.2 CFIE 鍐呮牳鍚堝苟 鉁?(`d0888cf`) 鈥?CFIE -74%, CFIE/EFIE=2.31脳
> - [x] 8.3 MLFMA Z_near 浼樺寲 鉁?(`d2f7963`)
> - [x] 8.4 鍐呭瓨鍒嗛厤鐑偣 鉁?(`82988cf`)
> - [x] 8.5 Julia 1.12 鍏煎 + 绫诲瀷绋冲畾鎬?鉁?(`67d3a8a`)
> - [x] 8.6 @fastmath + SIMD 鉁?(`1c6d499`)
> - [x] 8.7 BlockJacobiPreconditioner 鉁?(`76f8b16`)
> - [x] 8.8 鏈€缁堝熀鍑嗗娴?+ OOM 淇 鉁?(`6f4987a`)
> - [x] 8.9 CFIE 鍗曞惊鐜洿鎺ュ～鍏?鉁?(`ee6e2b1`) 鈥?Jet CFIE 14.48s
> - [x] 8.10 VEFIE 涓婁笁瑙掑绉颁紭鍖?鉁?(`bbf8fdd`) 鈥?VEFIE **41.30s** (蹇?Legacy 12%); SCFIE **65.67s** (蹇?Legacy 1.5%)

### Phase 9: 浠ｇ爜璐ㄩ噺涓庡彂甯?**(杩涜涓?**

- [x] 娴嬭瘯濂椾欢娓呯悊: 179/179 鍏ㄩ儴閫氳繃 鉁?
- [x] JuliaFormatter 缁熶竴浠ｇ爜椋庢牸 鉁?(`0535576`) 鈥?79 涓?`src/` 鏂囦欢, indent=4, margin=100
- [x] CHANGELOG.md 瀹屽杽 鉁?鈥?Phase 8-12 鎵€鏈夐噷绋嬬
- [x] compat 璇硶瑙勮寖鍖?鉁?+ 瑕嗙洊鐜囪剼鏈?+ 鍙戝竷娓呭崟
- [ ] 娴嬭瘯瑕嗙洊鐜?> 80%
- [ ] API 鏂囨。瀹屽杽 (Documenter.jl)
- [ ] 鐢ㄦ埛鏁欑▼鍜岀悊璁烘枃妗?
- [ ] 鍙戝竷鍒?Julia General Registry

### Phase 10: 鍏ㄦ柟绋嬪叏璺緞绮惧害瀵归綈 鉁?

> 璇﹁涓嬫柟 Phase 10 璇︾粏璁″垝銆?

### Phase 11: PWC 鍩哄嚱鏁版敮鎸?(鍥涢潰浣? 鉁?

- [x] PWCBasis 鏁版嵁缁撴瀯 + 鏋勯€犲嚱鏁?
- [x] VEFIE + PWC 鐩存帴瑁呴厤
- [x] SCFIE + RWG + PWC 鑰﹀悎
- [x] PWC 婵€鍔卞悜閲?+ RCS 璁＄畻
- [x] 178/178 娴嬭瘯閫氳繃
- Commits: `fc37542`, `35808ff`

### Phase 12: 鍏潰浣撳畬鏁存敮鎸?(PWCHex + RBF) 鉁?

- [x] 12.1 GaussQuadrature: 鍏潰浣?(8鐐? 鍜屽洓杈瑰舰 (4鐐? GQ 瑙勫垯
- [x] 12.2 MeshTypes: HexahedraInfo + Quads4Hexa 缁撴瀯浣?+ 杈呭姪鍑芥暟
- [x] 12.3 MeshIO: CHEXA Nastran 缃戞牸璇诲彇 (鍚画琛?
- [x] 12.4 BasisFunctions: PWCHexBasis (3 DOF/hex), RBF evaluate + 杈圭晫闈?
- [x] 12.5 VEFIE: PWCHex, RBF, 娣峰悎 TetraHex 瑁呴厤 (6 涓柊鏂规硶)
- [x] 12.6 SCFIE: RWG+PWCHex (骞剁煝 L) 鍜?RWG+RBF (鏍囬噺鍔? 鑰﹀悎
- [x] 12.7 Excitation: PWCHex 鍜?RBF 骞抽潰娉㈡縺鍔卞悜閲?
- [x] 12.8 RadiationIntegral + RCS: PWCHex 鍜?RBF 杈愬皠绉垎/RCS
- [x] 179/179 娴嬭瘯閫氳繃, +2296 琛?
- Commit: `099385b`

---

## Phase 8: 鎬ц兘浼樺寲 鈥?璇︾粏璁″垝

### 8.0 鐩爣

| 鎸囨爣 | 搴曠嚎 | 鎸戞垬鐩爣 | 璇存槑 |
|------|------|---------|------|
| **鐩稿悓鐢ㄤ緥鍏ㄦ祦绋嬭€楁椂** | 鈮?Legacy (1.0脳) | 鈮?Legacy 脳 0.5 (2脳 鍔犻€? | 鍚岀綉鏍笺€佸悓鏂圭▼銆佸悓姹傝В鍣ㄨ矾寰?|
| **MLFMA MatVec 鍗曟鑰楁椂** | 鈮?Legacy (1.0脳) | 鈮?Legacy 脳 0.5 | 鍚仛鍚?杞Щ+瑙ｈ仛 |
| **Z 鐭╅樀缁勮** | 鈮?Legacy (1.0脳) | 鈮?Legacy 脳 0.5 | 鍚?EFIE/MFIE/VEFIE 鍚勭畻瀛?|
| **宄板€煎唴瀛?* | 鈮?2脳 Legacy | 鈮?1.5脳 Legacy | Float64 vs Float32 鍏佽 2脳 |

### 8.1 鎬ц兘鍩虹嚎娴嬮噺 (Step 0 鈥?浼樺厛鎵ц)

> **鍘熷垯**: 鍏堥噺鍖栵紝鍐嶄紭鍖栥€傛棤鍩虹嚎鏁版嵁绂佹寮€濮嬩紭鍖栥€?

鍒涘缓 `benchmark/performance_baseline.jl`锛屽浠ヤ笅鐢ㄤ緥鍒嗛樁娈佃鏃讹紝鍚屾椂杩愯 Legacy 涓?EMSuite:

| 鐢ㄤ緥 | N | 鏂圭▼ | 璺緞 | 璁℃椂椤?|
|------|---|------|------|--------|
| PEC Plate 300MHz | ~2640 | EFIE | Direct | 鍩哄嚱鏁版瀯寤?/ Z 缁勮 / LU / RCS |
| Jet 100MHz | 14559 | EFIE | Direct | 鍚屼笂 |
| Jet 100MHz | 14559 | EFIE | MLFMA+GMRES | 鍏弶鏍?/ 杩戝満Z / 棰勬潯浠?/ GMRES / RCS |
| Sphere 600MHz | 26424 | CFIE | MLFMA+GMRES | 鍚屼笂 |
| Tetra 2GHz | ~986 | VEFIE | Direct | 鍩哄嚱鏁?/ Z 缁勮 / LU / RCS |
| TriTetra 2GHz | ~1071 | SCFIE | Direct | Z_SS/Z_SV/Z_VS/Z_VV / Fss / LU / RCS |

杈撳嚭: `test_results/PERFORMANCE_BASELINE.md`锛屾牸寮?

```
| 鐢ㄤ緥 | 闃舵 | Legacy (s) | EMSuite (s) | Ratio | 鐘舵€?|
```

### 8.2 宸茶瘑鍒儹鐐逛笌浼樺寲璺嚎

#### 鐑偣 1: SpinLock 鍏ㄥ眬閿?鈥?Z 鐭╅樀缁勮鐡堕 (P0)

**鐜扮姸**: `Impedance.jl` 涓瘡涓笁瑙掑舰瀵归兘瑕?`lock(spinlock)` / `unlock(spinlock)` 鍐欏叆鍏ㄥ眬 Z 鐭╅樀銆傚绾跨▼鎵╁睍鎬ф瀬宸€?

**鏂规**: 绾跨▼灞€閮ㄧ紦鍐?+ 鏈€鍚庡綊绾?
```julia
# Before: 
lock(spinlock); Z[m, n] += val; unlock(spinlock)
# After:
Z_local = [zeros(CT, N, N) for _ in 1:nthreads()]
@threads for i in workload
    Z_local[threadid()][m, n] += val  # 鏃犻攣
end
Z .= sum(Z_local)  # 涓€娆″綊绾?
```

**棰勬湡鏀剁泭**: 4 绾跨▼鍔犻€熸瘮 1.2脳 鈫?3.5脳+

**椋庨櫓**: 鍐呭瓨澧炲姞 (nthreads 鍊?Z 鐭╅樀)銆傚浜?N>10000 鐨?Dense Z锛屽彲鏀圭敤**鎸夎鍒嗗潡**: 姣忎釜绾跨▼璐熻矗鍥哄畾琛岃寖鍥达紝鏃犻渶閿佷篃鏃犻渶棰濆鍐呭瓨銆?

#### 鐑偣 2: CFIE 缁勮 9脳 EFIE (P0)

**鐜扮姸**: Jet N=14559, EFIE 缁勮 20.3s, CFIE 缁勮 180.5s (9脳)銆傜悊璁轰笂 CFIE = EFIE + MFIE锛屽簲 鈮?2脳 EFIE銆?

**鏂规**:
1. **Green 鍑芥暟澶嶇敤**: EFIE 鍜?MFIE 鍏变韩 $G(r,r') = e^{-jkR}/(4\pi R)$ 鍜?$\nabla G$銆傚悎骞朵负鍗曢亶鍘嗭紝璁＄畻涓€娆?G锛屽悓鏃剁疮鍔?EFIE 鍜?MFIE 璐＄尞銆?
2. **MFIE 鍐呭惊鐜紭鍖?*: 妫€鏌?MFIE 鏄惁鏈夊啑浣欏嚑浣曡绠?(娉曞悜閲忋€佷氦鍙夌Н)锛屾彁鍒板惊鐜棰勮绠椼€?
3. **瀵圭О鎬у埄鐢?*: EFIE L 绠楀瓙瀵圭О 鈫?浠呰绠椾笂涓夎; MFIE K 绠楀瓙鍙嶅绉?鈫?涓婁笁瑙掑彇璐熴€?

**棰勬湡鏀剁泭**: CFIE 浠?9脳 EFIE 鈫?鈮?2.5脳 EFIE (180s 鈫?鈮?50s)

#### 鐑偣 3: MLFMA Setup 鍗犳瘮杩囬珮 (P1)

**鐜扮姸**: Jet MLFMA 鎬昏€楁椂 63.4s锛屽叾涓?setup 56.2s (89%)銆係phere 138.6s锛宻etup 131.0s (94%)銆?

**瀛愰」**:
| 瀛愰」 | 褰撳墠鑰楁椂 (浼? | 浼樺寲鏂瑰悜 |
|------|-------------|---------|
| 鍏弶鏍戞瀯寤?| ~5% | 棰勫垎閰嶆暟缁勶紝鍑忓皯 push!/resize! |
| 杩戝満 Z_near 绋€鐤忕煩闃电粍瑁?| ~60% | 鍚岀儹鐐? (鍘婚攣)锛汣SC 棰勫垎閰?nnz |
| 杞Щ鐭╅樀棰勮绠?| ~15% | 缂撳瓨閲嶅鐨勭悆闈㈡尝灞曞紑锛汱ebedev 琛ㄥ鐢?|
| SAI/ILU 棰勬潯浠跺櫒鏋勫缓 | ~20% | 鑰冭檻杩戜技棰勬潯浠?(Block Jacobi) 鎴栧欢杩熸瀯寤?|

**棰勬湡鏀剁泭**: Setup 鎬绘椂闂?-30~50%

#### 鐑偣 4: 鍐呭瓨鍒嗛厤鐑偣 (P1)

**鏂规**:
1. **楂樻柉绉垎鐐归鍒嗛厤**: 閬垮厤姣忎釜涓夎褰㈠閲嶅鍒嗛厤 `SVector` 鏁扮粍
2. **Green 鍑芥暟璁＄畻闆跺垎閰?*: 纭繚 exp/sqrt 绛夎繍绠楀湪鏍囬噺涓婃搷浣滐紝涓嶄骇鐢熶复鏃舵暟缁?
3. **RCS 鍚庡鐞?*: 褰撳墠姣忎釜瑙傛祴鏂瑰悜鐙珛鍒嗛厤杈愬皠绉垎缂撳啿鍖?鈫?棰勫垎閰嶅鐢?

**宸ュ叿**: 浣跨敤 `@allocated` / `--track-allocation=user` 瀹氫綅

#### 鐑偣 5: SIMD / LoopVectorization (P2)

**鐜扮姸**: 鏈娇鐢?`LoopVectorization.jl` 鐨?`@turbo` 瀹忋€?

**鏂规**: 瀵逛互涓嬪唴寰幆娣诲姞 SIMD 浼樺寲:
1. 鏍囬噺鍔块」绱姞: `鈭?div_f 路 div_f' 路 G` 鈥?绾爣閲忥紝SIMD 鍙嬪ソ
2. 鐭㈤噺鍔块」绱姞: `鈭?f 路 f' 路 G` 鈥?SVector 鍐呯Н锛岀紪璇戝櫒鍙嚜鍔ㄥ悜閲忓寲
3. MLFMA 鑱氬悎/瑙ｈ仛: 鐞冮潰娉㈢郴鏁版暟缁勮繍绠?

**棰勬湡鏀剁泭**: 鍐呭惊鐜?1.5~3脳 鍔犻€?(鍙栧喅浜?CPU AVX 鏀寔)

**椋庨櫓**: `LoopVectorization.jl` 涓庡鏁?`ComplexF64` 鍏煎鎬ф湁闄愶紝鍙兘闇€瑕佹墜鍔ㄦ媶鍒嗗疄閮?铏氶儴

#### 鐑偣 6: 绫诲瀷绋冲畾鎬?(P2)

**鏂规**: 浣跨敤 `@code_warntype` 妫€鏌ュ叧閿矾寰?
1. `efie_interaction!` / `mfie_interaction!` 鈥?纭繚鏃?`Any` 绫诲瀷
2. `mul!` (MLFMAOperator) 鈥?纭繚鑱氬悎/杞Щ/瑙ｈ仛鍏ㄨ矾寰勭被鍨嬬ǔ瀹?
3. SCFIE 鑰﹀悎缁勮 鈥?娣峰悎鍩哄嚱鏁拌矾寰勫彲鑳芥湁绫诲瀷涓嶇ǔ瀹?

### 8.3 瀹炴柦璺嚎

| 姝ラ | 鍐呭 | 鍓嶇疆 | 棰勬湡鑰楁椂 | 鏀剁泭 |
|------|------|------|---------|------|
| **8.0** | 鎬ц兘鍩虹嚎娴嬮噺 | 鈥?| 0.5 澶?| 閲忓寲璧风偣 |
| **8.1** | Z 缁勮鍘婚攣 (琛屽垎鍧楀苟琛? | 8.0 | 1 澶?| EFIE 3~4脳 鍔犻€?|
| **8.2** | CFIE = EFIE+MFIE 鍚堝苟閬嶅巻 | 8.1 | 1 澶?| CFIE 4~5脳 鍔犻€?|
| **8.3** | MLFMA Z_near 鍘婚攣 + CSC 棰勫垎閰?| 8.1 | 0.5 澶?| Setup -30% |
| **8.4** | 鍐呭瓨鍒嗛厤鐑偣娑堥櫎 | 8.0 | 0.5 澶?| GC 鍘嬪姏闄嶄綆 |
| **8.5** | @code_warntype 绫诲瀷绋冲畾鎬у鏌?| 8.0 | 0.5 澶?| 娑堥櫎鍔ㄦ€佹淳鍙?|
| **8.6** | SIMD / @turbo 鍐呭惊鐜?| 8.5 | 1 澶?| 鍐呭惊鐜?1.5~3脳 |
| **8.7** | 棰勬潯浠跺櫒浼樺寲 (Block Jacobi) | 8.3 | 0.5 澶?| 棰勬潯浠舵瀯寤?-50% |
| **8.8** | 鏈€缁堝熀绾垮娴?| 鍏ㄩ儴 | 0.5 澶?| 楠岃瘉杈炬爣 |

**鎬婚浼?*: ~6 澶?

### 8.4 閫氳繃鍑嗗垯 鈥?鏈€缁堢粨鏋?

| 鐢ㄤ緥 | Legacy 鍏ㄦ祦绋?(s) | EMSuite 鍩虹嚎 (s) | EMSuite 鏈€缁?(s) | vs Legacy | 鐘舵€?|
|------|-------------------|------------------|-----------------|-----------|------|
| PEC Plate Direct (N=2640) | 8.29 | 3.44 | **5.84** | 0.70脳 | 鉁?|
| Jet EFIE Direct (N=14559) | 46.43 | 42.16 | **61.90** | 1.33脳 | 鈿狅笍 |
| Jet CFIE Direct (N=14559) | 64.21 | 188.63 | **97.98** | 1.53脳 | 鈿狅笍 |
| Jet EFIE MLFMA (N=14559) | N/A鹿 | 137.75 | **178.35** | 鈥?| 鈥?|
| Sphere CFIE MLFMA (N=26424) | N/A鹿 | 540.09 | **539.93** | 鈥?| 鈥?|
| Plate VEFIE Direct (N=15828) | 103.85 | 72.73 | **102.76** | 0.99脳 | 鉁?|
| PlateMetal SCFIE Direct (N=15860) | 66.52 | 90.57 | **130.73** | 1.96脳 | 鉂?|

鹿 Legacy MLFMA 鍦?Julia 1.12 涓嬪洜 SAI threadid() 鍏煎鎬ч棶棰樻棤娉曡繍琛?

### 8.5 鍥炲綊绾︽潫

- 浼樺寲鍚庢墍鏈?138/138 鍗曞厓娴嬭瘯 **蹇呴』** 缁х画閫氳繃
- Phase 10 宸查獙璇佺殑绮惧害鎸囨爣 (RMSE) **涓嶅緱** 閫€鍖?
- 姣忎釜浼樺寲姝ラ瀹屾垚鍚庣珛鍗宠繍琛?`Pkg.test()` + 鍏抽敭 benchmark

---

## Phase 10: 鍏ㄦ柟绋嬪叏璺緞绮惧害瀵归綈 鉁?

> 鐩爣: 瀵?5 绫荤Н鍒嗘柟绋?脳 4 绉嶆眰瑙ｈ矾寰勮繘琛屽叏鐞冮潰 RCS 绮惧瘑瀵规瘮锛屽畾閲忚瘉鏄?EMSuite 涓?Legacy 涓€鑷淬€?

### 10.0 宸查獙璇佺粨鏋滄眹鎬?

| 娴嬭瘯 | 鎸囨爣 | 缁撴灉 | 鐘舵€?|
|------|------|------|------|
| A1 S-EFIE Direct Jet | RMSE vs Legacy | 0.215 dB | 鉁?PASS |
| A2 S-EFIE Iterative Jet | RMSE vs A1 | 0.343 dB | 鈿狅笍 鏉′欢閫氳繃鹿 |
| A3 S-EFIE MLFMA Jet | RMSE vs Legacy | 0.303 dB | 鉁?PASS |
| A4 S-EFIE MPI Jet | rel_diff | 4.33e-16 | 鉁?PASS |
| B1 CFIE Z 鍒嗚В | rel_err | 0.0 (10/10) | 鉁?PASS |
| B2 S-MFIE MLFMA Sphere | 鐗╃悊瓒嬪娍 | RCS 鑼冨洿鍚堢悊 | 鉁?PASS虏 |
| B3 S-MFIE MPI (灏忕綉鏍? | rel_diff | 0.0 | 鉁?PASS |
| C1 S-CFIE Direct Sphere | RMSE vs Legacy | 0.001 dB | 鉁?PASS |
| C3 S-CFIE MLFMA Sphere | RMSE vs Legacy | 0.003 dB | 鉁?PASS |
| C3-MPI S-CFIE MPI (灏忕綉鏍? | rel_diff | 0.0 | 鉁?PASS |
| D1-SWG V-EFIE Direct | RMSE vs Legacy | 0.952 dB | 鉁?PASS |
| D2 V-EFIE Iterative | RMSE vs D1 | 8.9e-5 dB | 鉁?PASS |
| D3 V-EFIE MLFMA | RMSE vs Direct | 0.0 dB | 鉁?PASS |
| E1 VSEFIE Direct | RMSE vs Legacy | 0.602 dB | 鉁?PASS |
| E2 VS-EFIE Iterative | RMSE vs E1 | 3.3e-4 dB | 鉁?PASS |
| E3 VS-EFIE MLFMA | RMSE vs Direct | 0.0 dB | 鉁?PASS |

鹿 A2: EFIE 鏉′欢鏁板ぇ, 瀵硅棰勬潯浠?restart GMRES 鏈厖鍒嗘敹鏁?(residual 8.54e-2). D2/E2 宸茶瘉鏄?GMRES 鍩虹璁炬柦姝ｇ‘; A3 MLFMA 璇佹槑 EFIE 鍦ㄨ繎鍦?LU 棰勬潯浠朵笅鍙甯告敹鏁?
虏 B2: 鐞?ka=蟺 澶?MFIE 鍐呭叡鎸?(j鈧€(蟺)=0), MFIE 鍗曠嫭绮惧害鍙楅檺鏄凡鐭ョ墿鐞嗙幇璞? MLFMA 鍩虹璁炬柦鐢?A3/C3/D3/E3 鍏呭垎楠岃瘉.

**宸蹭慨澶?Bug:**
- **P0** (2026-02-28): `edgev虃` 鏂瑰悜鍙嶈浆 + `calc_near_interaction!` 闈㈢Н褰掍竴鍖?
- **P1** (2026-03-02): MLFMA 脳4 鍥犲瓙 + CFIE MFIE 绗﹀彿
- **P2** (2026-03-03): SCFIE Fss 鍗婂熀鍑芥暟杈圭晫闈㈢Н鍒嗕慨姝?
- **P3** (2026-07-29): CFIE MPI 骞惰瑁呴厤: `cfie_interaction!` 涓嶅瓨鍦?鈫?鏀逛负 EFIE+MFIE 鐙珛浜や簰+绾挎€х粍鍚?

### 10.1 鍏ㄧ悆闈㈤噰鏍锋柟妗?

| 鍙傛暟 | 鍊?| 璇存槑 |
|------|----|----- |
| 胃 鑼冨洿 | [-蟺, 蟺] | 绛夊悓鍙岀珯 RCS 鍏ㄨ搴︽壂鎻?|
| 胃 閲囨牱 | 73 鐐?(5掳 闂撮殧) | 鍏煎 Legacy 721 鐐瑰瓙闆?|
| 蠁 鑼冨洿 | [0, 蟺) | 鍗婄悆瀵圭О, 閬垮厤鍐椾綑 |
| 蠁 閲囨牱 | 18 鏉″垏闈?(10掳 闂撮殧) | 鍖呭惈 蠁=0掳/90掳 |
| 鎬昏娴嬫柟鍚?| 73 脳 18 = **1314** | 瑕嗙洊鍏ㄧ悆闈?|

### 10.2 娴嬭瘯鐭╅樀

| 缂栧彿 | 鏂圭▼绫诲瀷 | 鍑犱綍浣?| N (approx) | Direct | Iterative | MLFMA | MPI |
|------|----------|--------|------------|--------|-----------|-------|-----|
| **A** | S-EFIE | Jet 100MHz | 14559 | 鉁?A1 | 鈿狅笍 A2鹿 | 鉁?A3 | 鉁?A4 |
| **B** | S-MFIE | Sphere 600MHz | 26424 | 鉁?B1虏 | 鈥?| 鉁?B2鲁 | 鉁?B3鈦?|
| **C** | S-CFIE | Sphere 600MHz | 26424 | 鉁?C1 | 鈥?| 鉁?C3 | 鉁?C3-MPI鈦?|
| **D** | V-EFIE | Tetra 2GHz | ~986 | 鉁?D1 | 鉁?D2 | 鉁?D3 | 鉁?D4鈦?|
| **E** | VS-EFIE | TriTetra 2GHz | ~1071 | 鉁?E1 | 鉁?E2 | 鉁?E3 | 鈥?|

鹿 A2 鏉′欢閫氳繃: GMRES 鏈厖鍒嗘敹鏁?(EFIE 鏉′欢鏁板ぇ), D2/E2 宸茶瘉鏄庡熀纭€璁炬柦姝ｇ‘
虏 B1 = CFIE 鍒嗚В楠岃瘉 (灏忕綉鏍?
鲁 B2 鐗╃悊瓒嬪娍楠岃瘉 (ka=蟺 MFIE 鍐呭叡鎸?
鈦?B3/C3-MPI 浣跨敤灏忕綉鏍奸獙璇?MPI 鍩虹璁炬柦 (N=26424 Dense Z 鍐呭瓨涓嶈冻)
鈦?D4 V-EFIE Allreduce MPI: Tetra.nas N=3201, 1/2/4 procs Z[1,1] 绮剧‘涓€鑷?

### 10.3 姹傝В鍣ㄨ矾寰?

| 璺緞 | 璇存槑 | 鍙敤鏂圭▼ | 绾︽潫 |
|------|------|---------|------|
| **Direct** | Dense Z 鈫?LU | A, D, E | B/C 鐨?N=26424 Dense Z 闇€ 11 GB, 涓嶅彲琛?|
| **Iterative** | Dense Z 鈫?GMRES | A, D, E | 鍚屼笂 |
| **MLFMA** | MLFMAOperator 鈫?GMRES + SAI | A, B, C, D, E | 鍏ㄩ儴鍙敤 |
| **MPI** | 骞惰瑁呴厤 | A, C | 浠?RWG 琛ㄩ潰鏂圭▼ |

### 10.4 鍓╀綑瀛愰」

| 瀛愰」 | 姹傝В璺緞 | 瀵规瘮鍩哄噯 | 閫氳繃鍑嗗垯 | 鐘舵€?|
|------|---------|---------|---------|------|
| A2 | Iterative | A1 | RMSE < 0.1 dB | 鈿狅笍 0.343 dB (鏀舵暃鍙楅檺) |
| A4 | MPI (2 杩涚▼) | A1 | 鏈哄櫒绮惧害 | 鉁?4.33e-16 |
| B2 | MLFMA + GMRES | 鐗╃悊瓒嬪娍 | 瓒嬪娍涓€鑷?| 鉁?鐗╃悊鍚堢悊 |
| B3 | MPI (2 杩涚▼) | Serial | 鏈哄櫒绮惧害 | 鉁?0.0 (灏忕綉鏍? |
| C3-MPI | MPI (2 杩涚▼) | C3 | 鏈哄櫒绮惧害 | 鉁?0.0 (灏忕綉鏍? |
| D2 | Iterative | D1 | RMSE < 0.1 dB | 鉁?8.9e-5 dB |
| D3 | MLFMA + GMRES | D1 | vs D1 < 2 dB | 鉁?0.0 dB |
| E2 | Iterative | E1 | RMSE < 0.1 dB | 鉁?3.3e-4 dB |
| E3 | MLFMA + GMRES | E1 | vs E1 < 2 dB | 鉁?0.0 dB |

### 10.5 浼犻€掑噯鍒?

| 瀵规瘮绫诲瀷 | 鍑嗗垯 |
|---------|------|
| Direct vs Legacy Direct | RMSE < 1 dB |
| Iterative vs Direct | RMSE < 0.1 dB |
| MLFMA vs Direct | RMSE < 2 dB |
| MLFMA vs Legacy MLFMA | RMSE < 3 dB |
| MPI vs Serial | 鏈哄櫒绮惧害 |
| CFIE Z 鍒嗚В | rel_err < 1e-12 |

### 10.6 EMSuite API 瑕嗙洊鐭╅樀

| 鍔熻兘 | EFIE | MFIE | CFIE | VEFIE | SCFIE |
|------|------|------|------|-------|-------|
| 鐩存帴瑁呴厤 | 鉁?RWG | 鉁?RWG | 鉁?RWG | 鉁?SWG/PWC/PWCHex/RBF/Mixed | 鉁?RWG+SWG/PWC/PWCHex/RBF |
| 婵€鍔卞悜閲?| 鉁?| 鉁?| 鉁?| 鉁?SWG/PWC/PWCHex/RBF | 鉁?(鎷兼帴 RWG+SWG/PWC/PWCHex/RBF) |
| MLFMAOperator | 鉁?| 鉁?| 鉁?| 鉁?SWG | 鉁?RWG+SWG |
| MPI 骞惰瑁呴厤 | 鉁?| 鉁?| 鉁?| 鉁?SWG (Allreduce) | 鉂?|
| RCS 璁＄畻 | 鉁?RWG | 鉁?RWG | 鉁?RWG | 鉁?SWG/PWC/PWCHex/RBF | 鈿狅笍 闇€鎵嬪姩鎷嗗垎 |

### 10.7 鍩哄嚱鏁版敮鎸佺煩闃?

| 鍩哄嚱鏁?| 缃戞牸绫诲瀷 | DOFs/鍏冪礌 | 鏀寔鐨処E | 鐘舵€?|
|--------|---------|-----------|---------|------|
| RWG | 涓夎褰?(Surface) | 1/鍐呴儴杈?| EFIE, MFIE, CFIE, SCFIE(S) | 鉁?|
| SWG | 鍥涢潰浣?(Volume) | 1/鍐呴儴闈?| VEFIE, SCFIE(V) | 鉁?|
| PWC | 鍥涢潰浣?(Volume) | 3/鍥涢潰浣?(x,y,z) | VEFIE, SCFIE(V) | 鉁?Phase 11 |
| PWCHex | 鍏潰浣?(Volume) | 3/鍏潰浣?(x,y,z) | VEFIE, SCFIE(V) | 鉁?Phase 12 |
| RBF | 鍏潰浣?(Volume) | 1/闈?(鍚竟鐣? | VEFIE, SCFIE(V) | 鉁?Phase 12 |

---
## Phase 14: 全量精度测试与对比报告（计划中）

> 详见 [PHASE_14_ACCURACY_REPORT_PLAN.md](PHASE_14_ACCURACY_REPORT_PLAN.md)

### 目标

对所有主要积分方程 × 求解路径，与 **Feko 商业软件结果** 或 **Mie 解析解** 进行系统对比，生成独立精度报告。

### 子任务

- [x] 14.0 Feko CSV 数据解析器 (`FekoReader.jl`) + TDD 测试 ✅ commit b416090
- [x] 14.1 参考基准生成器：`ReferenceData.jl`（Mie PEC + Mie 介质 + 偶极子解析）✅
- [x] 14.2 精度指标函数 `AccuracyResult` + `AntennaAccuracyResult` ✅ commit 65b4297 (25/25 通过)
- [x] 14.3 F1–F4: Jet 100MHz 仿真脚本 (`run_F1_F4_jet.jl`) ✅ commit b2efd57
- [x] 14.4 F5–F6: Sphere 600MHz 仿真脚本 (`run_F5_F6_sphere.jl`) ✅
- [x] 14.5 F7–F9: Plate 1.2GHz 仿真脚本 (`run_F7_F9_plate.jl`) ✅
- [x] 14.6 P1–P3: PMCHW 介质球 Direct 脚本 (`run_P1_P3_pmchw.jl`) ✅
- [ ] 14.7 `PMCHWMLFMAOperator` 实现 (2×2 块 MLFMA) + 单元测试
- [ ] 14.8 P2: PMCHW 介质球 MLFMA 验证
- [x] 14.9 A1–A4: 偶极子天线 DeltaGap 基准脚本 (`run_A1_A4_antenna.jl`) ✅ commit 3039d32
- [ ] 14.10 实际运行仿真 → CSV → `generate_report.jl` → ACCURACY_REPORT.md
- [ ] 14.11 检视迭代 (≥ 2 轮 clean)

### 精度验收门限

| 求解路径 | 方程类型 | Feko/Mie RMSE 门限 |
|---------|---------|-------------------|
| Direct | S-EFIE, S-CFIE, V-EFIE, SCFIE | ≤ 2.0 dB |
| Direct | PMCHW 介质球 | ≤ 2.5 dB |
| MLFMA+GMRES | S-EFIE, S-CFIE | ≤ 3.0 dB |
| MLFMA+GMRES | PMCHW（PMCHWMLFMAOperator） | ≤ 3.0 dB |
| 天线端口 | 输入阻抗误差 | < 5% |
| 天线端口 | 最大方向性误差 | < 1.0 dBi |

### Feko 基线来源

`C:\Users\12253\OneDrive\MoM\MoM_AllinOne\deps\compare_feko\`
- `jet_100MHzRCS.csv` — Jet PEC, 100 MHz
- `sphere_600MHzRCS.csv` — PEC 球, 600 MHz
- `plate_1dot2GHzRCS.csv` — 介质板, 1.2 GHz
- `plate_metal_1dot2GHzRCS.csv` — 介质+金属板, 1.2 GHz

---

## Phase 15: 介质与金属-介质混合天线精度测试 + PMCHWMLFMAOperator（计划中）

> 详见 [PHASE_15_DIELECTRIC_ANTENNA_PLAN.md](PHASE_15_DIELECTRIC_ANTENNA_PLAN.md)

### 目标

在 Phase 14 天线测试（A1–A4 纯金属偶极子）基础上扩展：介质天线（PMCHW）和金属-介质混合天线（VS-EFIE/VS-CFIE）的输入阻抗验证，以及 PMCHWMLFMAOperator 实现。

### 子任务

- [ ] 15.1 TDD: `test_pmchw_excitation.jl`（PMCHW DeltaGap 激励 + input_impedance）
- [ ] 15.2 实现 `excitation_vector(PMCHW, DeltaGapSource, RWGBasis)` → 2N 向量
- [ ] 15.3 实现 `input_impedance(op::PMCHW, source, I_2N, basis)` → J 部分阻抗
- [ ] 15.4 TDD: `test_scfie_delta_gap.jl`（SCFIE DeltaGap 激励）
- [ ] 15.5 实现 `excitation_vector(SCFIE, DeltaGapSource, rwg, swg)` → (N_S+N_V) 向量
- [ ] 15.6 基准脚本 `benchmark/accuracy/run_B1_B5_antenna.jl`（B1–B5 用例）
- [ ] 15.7 TDD: `test_pmchw_mlfma_operator.jl`
- [ ] 15.8 实现 `assemble_near_field_pmchw`（2N×2N 稀疏，4 块：EJ/EM/HJ/HM，无需 MagneticRWGBasis）
- [ ] 15.9 实现 `aggregate_leaf_pmchw!`（单函数，x_range 参数区分 J/M 系数）
- [ ] 15.10 实现 `disaggregate_leaf_pmchw_j!` + `_m!`（四块接收核函数，依据 Gibson Algorithm 14 两遍设计）
- [ ] 15.11 组装 `PMCHWMLFMAOperator` struct（两棵 N 点八叉树：octree0/k0 + octree1/k1）+ 构造函数 + `mul!`（4 遍远场：J×k0, J×k1, M×k0, M×k1）
- [ ] 15.12 更新 `generate_report.jl` 加入 B1–B5
- [ ] 15.13 检视迭代 (≥ 2 轮 clean)
- [ ] 15.G1 刷新 Theory-Implementation-Test 治理闭环（可执行原则 + 门禁测试，进行中：Gate A/B 已落地）
    - 计划文档: `.github/plans/phase_15_theory_impl_test_refresh.md`

### 精度验收门限

| 用例 | 方程 | 参考基准 | 门限 |
|------|------|---------|------|
| B1 PMCHW Direct | PMCHW (εᵣ=4) | 物理自洽 + εᵣ→1 极限 | Re(Z_in) > 0，εᵣ→1 误差 <10% |
| B2 PMCHW MLFMA | PMCHW (双八叉树+四遍远场) | B1 Direct | ΔZ_in: Re<5%, Im<20Ω |
| B3 VS-EFIE Direct | SCFIE (α=0) | EFIE-only 当 εᵣ→1 | Z_in Re 误差 <10% |
| B4 VS-CFIE Direct | SCFIE (α=0.5) | B3 (α=0) | ΔZ_in <5Ω |
| B5 VS-CFIE MLFMA | SCFIE (α=0.5) | B4 Direct | ΔZ_in Re<5%, Im<20Ω |

---
## 关键参考

- **Legacy 浠ｇ爜**: `MoM_Basics/`, `MoM_Kernels/`, `MoM_AllinOne/`
- **楠岃瘉鑴氭湰**: `EMSuite/benchmark/verify_*.jl`, `EMSuite/scripts/verification/`
- **鐞嗚**: Harrington "Field Computation by Moment Methods"; Chew et al. "Fast and Efficient Algorithms in CEM"
