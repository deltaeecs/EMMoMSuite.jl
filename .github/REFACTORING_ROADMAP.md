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
- [x] 娴嬭瘯瑕嗙洊鐜?> 80%
- [x] API 鏂囨。瀹屽杽 (Documenter.jl)
- [x] 鐢ㄦ埛鏁欑▼鍜岀悊璁烘枃妗?
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
- [x] 14.7 PMCHW block/operator shell + MLFMA backend 主线实现（Dense shell 先行，MLFMA 后接） + 单元测试
- [x] 14.8 P2: PMCHW 介质球 MLFMA 验证
- [x] 14.9 A1–A4: 偶极子天线 DeltaGap 基准脚本 (`run_A1_A4_antenna.jl`) ✅ commit 3039d32
- [x] 14.10 实际运行仿真 → CSV → `generate_report.jl` → ACCURACY_REPORT.md
- [x] 14.11 检视迭代 (≥ 2 轮 clean)

### 精度验收门限

| 求解路径 | 方程类型 | Feko/Mie RMSE 门限 |
|---------|---------|-------------------|
| Direct | S-EFIE, S-CFIE, V-EFIE, SCFIE | ≤ 2.0 dB |
| Direct | PMCHW 介质球 | ≤ 2.5 dB |
| MLFMA+GMRES | S-EFIE, S-CFIE | ≤ 3.0 dB |
| MLFMA+GMRES | PMCHW（block/operator shell + MLFMA backend） | ≤ 3.0 dB |
| 天线端口 | 输入阻抗误差 | < 5% |
| 天线端口 | 最大方向性误差 | < 1.0 dBi |

### Feko 基线来源

`C:\Users\12253\OneDrive\MoM\MoM_AllinOne\deps\compare_feko\`
- `jet_100MHzRCS.csv` — Jet PEC, 100 MHz
- `sphere_600MHzRCS.csv` — PEC 球, 600 MHz
- `plate_1dot2GHzRCS.csv` — 介质板, 1.2 GHz
- `plate_metal_1dot2GHzRCS.csv` — 介质+金属板, 1.2 GHz

---

## Phase 15: 介质与金属-介质混合天线精度测试 + PMCHW Transmission Block/Operator 架构（计划中）

> 详见 [PHASE_15_DIELECTRIC_ANTENNA_PLAN.md](PHASE_15_DIELECTRIC_ANTENNA_PLAN.md)

### 目标

在 Phase 14 天线测试（A1–A4 纯金属偶极子）基础上扩展：介质天线（PMCHW）和金属-介质混合天线（VS-EFIE/VS-CFIE）的输入阻抗验证，并把 PMCHW 顶层实现迁移到 Bempp 风格 block/operator shell；Dense backend 先行，MLFMA 作为 backend 收编，H-matrix 延后。

### 子任务

- [x] 15.1 TDD: `test_pmchw_excitation.jl`（PMCHW DeltaGap 激励 + input_impedance）
- [x] 15.2 实现 `excitation_vector(PMCHW, DeltaGapSource, RWGBasis)` → 2N 向量
- [x] 15.3 实现 `input_impedance(op::PMCHW, source, I_2N, basis)` → J 部分阻抗
- [x] 15.4 TDD: `test_scfie_delta_gap.jl`（SCFIE DeltaGap 激励）
- [x] 15.5 实现 `excitation_vector(SCFIE, DeltaGapSource, rwg, swg)` → (N_S+N_V) 向量
- [x] 15.6 基准脚本 `benchmark/accuracy/run_B1_B5_antenna.jl`（B1–B5 用例）
- [x] 15.7 TDD: `test_pmchw_mlfma_operator.jl`
- [x] 15.8 实现 `assemble_near_field_pmchw`（2N×2N 稀疏，4 块：EJ/EM/HJ/HM，无需 MagneticRWGBasis）
- [x] 15.9 实现 `aggregate_leaf_pmchw!`（单函数，x_range 参数区分 J/M 系数）
- [x] 15.10 实现 `disaggregate_leaf_pmchw_j!` + `_m!`（四块接收核函数，依据 Gibson Algorithm 14 两遍设计）
- [x] 15.11 完成 PMCHW block/operator shell 与 Dense backend 验证，再把现有 `PMCHWMLFMAOperator` 收编为 MLFMA backend（兼容 facade 可保留）
- [x] 15.12 更新 `generate_report.jl` 加入 B1–B5
- [x] 15.13 检视迭代（连续 3 轮无新发现后通过）
- [x] 15.14 Alternative formulation：N-Muller dense baseline + 正式测试
- [x] 15.15 同一 dielectric sphere 上比较 PMCHW / N-Muller conditioning、直接解与 GMRES 行为
- [x] 15.G1 刷新 Theory-Implementation-Test 治理闭环（函数级开发合同 + 连续 3 轮文档检视）
    - 当前执行文档: `.github/plans/phase_15_theory_impl_test_refresh.md`
    - 强制路线: Gate D0（Dense PMCHW -> Mie 基线复验）→ Gate D1（近场逐元素）→ Gate B（单 pass）→ Gate C（端到端）
    - 强制要求: 禁止新增临时诊断脚本，诊断结论必须转正式测试或直接删除
    - 最新定位: `k1` 共享 core 的红灯已压缩到 leaf far-path 近距阈值；`leaf_size=0.10` 下 `near_range=4` 红、`near_range=7` 绿，正式构造函数已据此上调 near-range 标尺
    - 当前状态: 正式 `GD2`、`Gate B`、`Gate C` 已随 near-range 修正转绿；共享 core 的旧红灯保留为诊断回归，用于证明阈值机制而非重新作为主线阻塞
    - B2 当前状态: 输入阻抗实部误差 `3.885%`，已满足 `<5%` 验收门限
    - 新增中尺度回归: `N=540` 时 `Gate C` 仍为 green，但 `GMRES(200)` 对 `LU` 的 B2 偏差放大；dense+GMRES 与 MLFMA+GMRES 一致，后续应把问题归入求解器/预条件阶段，而非 MLFMA matvec 正确性阶段
    - 新增接口修复: `BlockJacobiPreconditioner` 已支持非连续索引块，`MLFMAOperator` / `PMCHWMLFMAOperator` 可直接构造 `BlockJacobiPreconditioner(op)` 与兼容签名 `BlockJacobiPreconditioner(op, basis)`
    - 最新边界修正: 上述“归入求解器/预条件阶段”的判断只对 `GMRES(200)` 低精度阶段成立；当 `restart` 提升到 `300/1000` 时，dense GMRES 已接近 `LU`，而 PMCHW MLFMA 仍保留约 `5.8%` 阻抗偏差，因此下一轮重点必须转向**长 Krylov 下的 PMCHW MLFMA 算子保真度**
    - 额外约束: 现有 `Diagonal` / `ILU(Z_near)` / `BlockJacobi(op)` 在 `N=540` 夹具上都比无预条件更差；继续扩大 `near_range` 则会触发 OOM，因此这两条都不能作为主线修复方案
    - 外部调研新增结论: Bempp 的 transmission/FMM 主线是 block multitrace + `singular_part + FMM evaluator`，FreeFEM/BEMTool/Htool 主线是 cluster-tree H-matrix；二者都不同于本地 `PMCHWMLFMAOperator` 的“双八叉树 + 四遍远场”专用实现，因此后续应把问题进一步拆到 operator/block 级误差，而不是只在整体 2N `mul!` 上调参
    - 新增 Gate S 主线: 必须补做 `dense weak` / `dense strong` / `MLFMA weak` / `MLFMA strong` 四路对照；只有在 `dense strong` 已回收而 `MLFMA strong` 仍失配时，才允许把剩余长 Krylov 偏差正式归因到 MLFMA 压缩误差
- [x] 15.G2 冻结 Phase 15 架构主线为 Bempp 风格 block/operator shell
    - 顶层 transmission 语义固定为四块 `EJ/EM/HJ/HM` block/operator shell，而不是继续扩张单体 `PMCHWMLFMAOperator`
    - `weak_form` / `strong_form` / mass-matrix-aware solve 归入 shell；Dense / MLFMA / H-matrix 统一视为 backend
    - 执行顺序冻结为：Dense shell → strong-form 分界（Gate S）→ MLFMA backend 收编 → 后续才允许评估 H-matrix backend
    - FreeFEM/Htool 路线被正式延后到 shell 与 Gate S 稳定之后，禁止在当前阶段直接跳入 H-matrix 主线
- [x] 15.G3 落地 PMCHW Dense block/operator shell 最小实现
    - 新增 `PMCHWBlockOperator` 与 `DensePMCHWBackend`，显式承载四块 `EJ/EM/HJ/HM` 语义
    - `weak_form` / `strong_form` 已在 shell 层落位；默认 `strong_form` 现已接入真实 RWG surface Gram 的 2N block pairing
    - 新增 `test/test_pmchw_operator_shell.jl`，验证 block split、默认 strong-form pairing、matvec 与 GMRES 兼容性
    - 当前完成的是 Dense shell + Dense pairing 基线；`Gate S` 的正式中尺度四路对照与 MLFMA backend 收编仍未完成
- [x] 15.G4 落地 Gate S 的 dense 半边回归
    - 新增 `strong_form_rhs` 与 `recover_trial_coefficients`，固定 strong-form 的 shell 层 RHS/解向量变换语义
    - 新增 `test/test_pmchw_gate_s_dense.jl`，在 `N=540` 夹具上正式比较 `dense weak` / `dense strong` / `dense LU`
    - 当前证据已确认 `dense strong` 在同一 GMRES 设置下优于 `dense weak`，但完整 Gate S 仍缺 `MLFMA weak` / `MLFMA strong`
- [x] 15.G5 完成 PMCHW MLFMA backend 的 shell 收编小夹具验证
    - 新增 `MatrixFreePMCHWBackend`，允许把 `PMCHWMLFMAOperator` 作为 matrix-free backend 接入 `PMCHWBlockOperator`
    - 新增 strong-form 变换包装算子，用于非矩阵 backend 的 `strong_form(shell)`
    - 新增 `test/test_pmchw_operator_shell_mlfma.jl`，验证 MLFMA backend 在 shell 下的 weak/strong 求解链可执行，并与 Dense shell 在小夹具上保持近似一致
    - 当前尚未完成的是 **中尺度** `MLFMA weak / strong` 两路的稳定正式回归；现阶段只有 Dense 半边在 `N=540` 夹具上已被正式锁定
- [x] 15.G6 锁定中尺度四路 Gate S 专门回归
    - 新增 `test/test_pmchw_gate_s_mlfma_medium.jl`，在同一 `N=540` 夹具上统一比较 `dense weak` / `dense strong` / `MLFMA weak` / `MLFMA strong`
    - 当前结果：`9/9 pass`，运行时间约 `9m05s`
    - 关键指标：`relw = 3.45e-3`、`rels = 3.18e-3`、`zgw = 2.34e-5`、`zgs = 5.75e-5`
    - 当前边界：该测试已把 Gate S 的中尺度四路主对照正式化，但因运行成本较高，暂不并入默认 `test/runtests.jl`
- [x] 15.G7 接通 B2 天线基准脚本的 PMCHW MLFMA 实跑路径
    - `benchmark/accuracy/run_B1_B5_antenna.jl` 的 B2 已从占位状态改为 shell + MatrixFree backend + strong-form GMRES
    - 单项复验结果：`B2_PMCHW_MLFMA` 已通过，`Z_in = +0.000 + j(+0.186) Ω`，相对 Direct 参考的虚部误差约 `0.01 Ω`
    - 近零参考实部下新增 `re_ref_floor_ohm` 判定下限，避免球面 delta-gap 场景中因 `Re(Z_ref)≈0` 触发伪失败
    - 当前边界：`15.6` 仍未整体勾选，因为 B5 仍待实现；但 B2 已不再是脚本中的空白占位
- [x] 15.G8 验证天线基准脚本默认入口可执行
    - 已复跑 `benchmark/accuracy/run_B1_B5_antenna.jl` 默认入口，`B1 / B3 / B4` 全部通过
    - 当前脚本状态：默认入口绿色，`B2` 单项绿色，剩余缺口仅为 `B5`
    - 这一步锁定了 Phase 15 天线基准框架本身已可用，后续实现重点可收束到 `B5` 或 alternative formulation / backend 误差预算
- [x] 15.G9 完整 B1–B5 天线基准全绿
    - 已复跑 `benchmark/accuracy/run_B1_B5_antenna.jl B1 B2 B3 B4 B5`，当前结果 `PASS: 6 / 6`
    - `B5` 已接通到 `SCFIE + MLFMAOperator + GMRES` 路径，并与 `B4 Direct` 在当前 TriTetra 夹具上保持一致
    - `benchmark/accuracy/generate_report.jl` 已修正 B5 标签并可正确汇总 `B1`–`B5`
    - 当前边界：`15.13` 检视迭代仍未完成；Phase 15 的下一实现重点已从天线基准脚本接线转向 alternative formulation / backend 误差预算
- [x] 15.G10 Alternative formulation 子流启动：N-Muller dense baseline
    - 已新增 `src/IntegralEquations/NMuller.jl`，并由 `src/IntegralEquations/IntegralEquations.jl` 与 `src/EMSuite.jl` 正式导出
    - 已新增 `test/test_nmuller.jl`，当前结果 `15/15 pass`
    - 已新增 `test/test_nmuller_comparison.jl` 与 `benchmark/compare_pmchw_nmuller_sphere.jl`，在共享球夹具上记录 formulation 行为差异
    - 当前观测：`cond(NMuller)/cond(PMCHW) = 1.4837e-2`，同一 GMRES 预算下 `NMuller` 的 residual 与 relative solution error 均低于 `PMCHW`
    - 执行文档：`.github/plans/phase_15_nmuller_dense_baseline.md`
- [x] 15.G11 N-Muller DeltaGap / input_impedance + PMCHW MLFMA 显式 budget 接口
    - 已为 `NMuller` 新增 `DeltaGapSource` 激励与 `input_impedance` 路径，形成最小诊断性馈电链；其物理端口语义仍待与 formulation-specific 参考实现对齐
    - 已新增 `test/test_nmuller_excitation.jl`，当前结果 `10/10 pass`
    - 已为 `PMCHWMLFMAOperator` 新增 `PMCHWMLFMAErrorBudget`，显式暴露 `near_range`、`L_min` 与 `leaf_size_eff` 控制，同时保持旧构造入口兼容
    - `test/test_pmchw_mlfma_operator.jl` 已新增 budget 回归，当前结果 `PMCHWMLFMAOperator budget interface preserves defaults and exposes overrides | 8/8 pass`
    - 最新诊断：`benchmark/compare_pmchw_nmuller_impedance.jl` 显示当前 N-Muller DeltaGap / input_impedance 组合与 PMCHW `Z_in` 相差 7-8 个数量级，说明这条链路仍处于语义排查阶段，暂不能升格为正式阻抗验收门
- [x] 15.G12 PMCHW MLFMA budget sweep 基准入口
    - 已新增 `benchmark/compare_pmchw_mlfma_budget.jl`，固定 small / medium 夹具比较不同 `PMCHWMLFMAErrorBudget` 配置
    - 基准输出项覆盖 `near_range`、`leaf_size_eff`、`L_min`、`nnz_near / near_density`、matvec 相对误差、强形式 `Z_in` 误差与构造/求解耗时
    - 对应结果会落盘到 `test_results/accuracy/PMCHW_MLFMA_budget_sweep_<preset>.csv`，用于后续冻结正式预算门
    - 首轮 `medium` 观测：budget 能把 `near_density` 从 `1.0000` 压到 `0.2350`–`0.3007`，同时把 `rel_matvec` 控制在 `7.3e-4`–`7.7e-4`；但固定 `GMRES(100)` 下 `Z_in` 误差几乎不变，当前说明 solver 截止仍是主导因素
- [x] 15.G13 Phase 15 检视迭代闭环（3 轮 clean）
    - Round 1（架构/算法审查）: 核对 shell 职责边界、四块语义、strong-form 路径、backend 收编一致性，未发现新的 Phase 15 阻塞问题
    - Round 2（关键回归）: `test_pmchw_operator_shell.jl`、`test_pmchw_operator_shell_mlfma.jl`、`test_pmchw_gate_s_dense.jl`、`test_nmuller.jl`、`test_nmuller_comparison.jl` 全部通过
    - Round 3（中尺度门禁复验）: `test_pmchw_gate_s_mlfma_medium.jl`、`test_pmchw_mlfma_budget_medium.jl`、`test_nmuller_planewave_gmres_trajectory_medium.jl` 全部通过
    - 残余风险（非 Phase 15）: `docs/setup_docs.jl`、`test/test_legacy_parity.jl` 仍存在语法错误，已记录为跨阶段遗留，不阻塞 Phase 15 关闭
- [x] 15.G13 PMCHW MLFMA medium budget 专门门禁
        - 已新增 `test/test_pmchw_mlfma_budget_medium.jl`，固定 `N=540` 夹具比较 `default / loose_near / fixed_leaf_0p04_nr9` 三组预算
        - 当前门禁锁定三类事实：
            1. budget 会显著改变 `near_density`
            2. default budget 的 matvec fidelity 优于更稀疏预算
            3. 固定 `GMRES(100)` 下三组预算的 `Z_in` 与最终残差几乎不变，说明当前偏差主要不由 budget 决定
        - 该回归运行成本较高，暂不并入默认 `runtests.jl`
- [x] 15.G14 PMCHW MLFMA 代表性预算长 Krylov 对照
    - 已新增 `benchmark/compare_pmchw_mlfma_budget_krylov.jl`，固定 `N=540` medium 夹具、固定三组代表性预算 `default / loose_near / fixed_leaf_0p04_nr9`
    - `benchmark/compare_pmchw_mlfma_budget.jl` 现默认只跑代表性三组预算；如需恢复探索集，显式传入 `full`
    - 已固定两档 Krylov 对照：`short = restart=100, maxiter=100, reltol=1e-4`，`long = restart=300, maxiter=600, reltol=1e-6`
    - 当前 `medium` 结果已落盘到 `test_results/accuracy/PMCHW_MLFMA_budget_krylov_medium.csv`
    - 关键结论：
        1. `short` 档下三组预算相对 dense 的 `Z_in` gap 仍只有亚欧姆量级（实部 `0.03`–`0.27Ω`，虚部 `0.18`–`0.54Ω`），预算影响仍被 Krylov 截止误差淹没
        2. 刷新后的 `long` 档三组预算相对 dense 的 `Z_in` gap 也都只剩亚欧姆量级：`default = 0.03815Ω / 0.00236Ω`，`loose_near = 0.48107Ω / 0.02029Ω`，`fixed_leaf_0p04_nr9 = 0.06282Ω / 0.07629Ω`
        3. 因此 G14 现已与 BG2、Arnoldi、GMRES 轨迹诊断重新对齐：当前 long-Krylov 主线不再表现为 budget 主导的大分叉
- [x] 15.G15 PMCHW MLFMA medium long-Krylov budget 专门门禁
    - 已新增 `test/test_pmchw_mlfma_budget_krylov_medium.jl`，固定 `N=540` 夹具、固定 long Krylov 配置 `restart=300, maxiter=600, reltol=1e-6`
    - 当前门禁只保留最有信息量的两组预算：`default` 与 `loose_near`
    - K 接收链修复后，旧门限已被证明过时；最新回归结果为 `15.BG2 PMCHW MLFMA medium long-Krylov budget gate | 16/16 pass | 29m08.0s`
    - 当前门禁锁定四类事实：
        1. dense strong 在同一 long Krylov 配置下已明显逼近 `LU`（`Re` 误差 `6.78%`, `Im` 误差 `4.28Ω`）
        2. `default` budget 相对 dense strong 的 `Z_in` gap 已收缩到 `0.03815Ω / 0.00236Ω`
        3. `loose_near` budget 虽然仍更差，但 gap 也仅 `0.48107Ω / 0.02029Ω`
        4. 因此 receive 修复后，long-Krylov 剩余边界已不能再归因为“budget 主导的 backend fidelity 大分叉”
- [x] 15.G16 PMCHW medium block/operator fidelity 诊断
    - 已新增 `benchmark/compare_pmchw_block_fidelity_medium.jl`，固定 `N=540` 夹具、`default / loose_near` 两组预算、`J_only / M_only` 两类物理 probe
    - 诊断语义已冻结为：strong-form 对照必须比较同一物理 probe；输入先过 `trial_pairing`，输出再映回测试空间，禁止直接把 weak-space probe 送入 `strong_form`
    - 结果已落盘到 `test_results/accuracy/PMCHW_block_fidelity_medium.csv`
    - 当前结论：
        1. 映回物理空间后，`strong` 与 `weak` 的 block fidelity 指标一致，因此当前 medium 边界不是 strong-form 专属失真
        2. `J_only` probe 在 `default` budget 下总体误差仅 `3.26e-5`，继续维持高保真
        3. 修复 PMCHW K 接收链后，`M_only` probe 在 `default` budget 下已降到 `2.98e-4`，`loose_near` 下也仅 `2.03e-3`
        4. `M×k0` 与 `M×k1` 单 pass 的 `E-row` 误差也分别降到 `2.87e-4 / 9.89e-4`（`default`）与 `1.06e-3 / 3.00e-3`（`loose_near`），说明此前 medium 边界确由 K-type receive 公式错误驱动
    - 后续主线：基于已修复的 `M` 通道接收链，重新评估 medium long-Krylov budget gap 是否仍由 backend fidelity 主导，而不是继续停留在旧的 `M_only -> E-row` 故障结论
- [x] 15.G17 PMCHW medium block fidelity 专门门禁
    - 已新增 `test/test_pmchw_block_fidelity_medium.jl`，固定 `N=540` 夹具，把 medium block fidelity 诊断升级为正式专门回归
    - 当前结果：K 接收链修复后已切换为精度上界门禁；最新回归结果为 `15.BF1 PMCHW medium block fidelity gate | 26/26 pass | 4m42.9s`
    - 当前门禁锁定三类事实：
        1. 映回物理空间后，`strong` 与 `weak` 的 fidelity 指标一致，因此当前 medium block gap 不属于 strong-only 问题
        2. `J_only` 与 `M_only` 在 `default` budget 下都已进入 `1e-4 ~ 1e-5` 级别，`loose_near` 下 `M_only` 仍只到 `2e-3` 量级
        3. `M×k0` 与 `M×k1` 两条 pass 在 `E-row` 上均已恢复到 `1e-3` 级别以内/附近，且 `loose_near` 相比 `default` 仍会稳定恶化，保留了 budget 灵敏度
    - 后续主线：把 medium long-Krylov 的剩余边界重新与 BG1/BG2 对照，判断 receive 修复后是否还需要进一步收缩 near/far budget 或 GMRES/预条件链
- [x] 15.G18 Receive 修复后的 budget 主线复核
    - 已复核 `test_results/accuracy/PMCHW_MLFMA_budget_sweep_medium.csv` 与 `test_results/accuracy/PMCHW_MLFMA_budget_krylov_medium.csv`
    - 结论已更新为：短 Krylov 结论保留，但旧的 long-Krylov 大 gap 归因被推翻：
        1. 短 Krylov `GMRES(100)` 下，代表性 budget 的 `Z_in` spread 仍很小，主导项仍是 solver 截止
        2. medium budget sweep 仍保持 `default > loose > fixed` 的 `near_density` 梯度，且 `default` matvec fidelity 最好（`3.20e-4 < 7.30e-4 < 7.70e-4`）
        3. long Krylov 下现行 BG2 值已收缩到：`default = 0.03815Ω / 0.00236Ω`，`loose_near = 0.48107Ω / 0.02029Ω`
        4. 因此 Phase 15 剩余边界已从“budget 主导”转向“更深 Krylov 轨迹 / formulation / conditioning 放大”
- [x] 15.G19 PMCHW medium Arnoldi 子空间保真诊断
    - 已新增 `benchmark/compare_pmchw_krylov_subspace_medium.jl`，按 dense strong-form 的 Arnoldi 子空间方向比较 `dense strong` 与 `MLFMA strong` 的 matvec
    - 当前结果已落盘到 `test_results/accuracy/PMCHW_krylov_subspace_medium.csv`
    - 当前结论：
        1. 在前 12 个 Arnoldi 方向上，`default` budget 的最大 `rel_total` 仅 `8.12e-5`，`loose_near` 也仅 `4.26e-4`
        2. `E-row` 误差与总体误差同量级，`H-row` 略大但仍仅 `3.48e-4`（default）/ `2.85e-3`（loose）
        3. 因此 BG2 的 long-Krylov 大 gap 不能简单归因于“早期 Krylov 子空间方向上的 fast matvec 已明显失真”
    - 后续主线：继续向更深 Krylov 子空间、重启/正交化轨迹或 formulation 条件数放大效应收缩，而不是回退到随机 probe 或已修复的 M-pass receive 链
- [x] 15.G20 PMCHW medium GMRES 轨迹诊断
    - 已新增 `benchmark/compare_pmchw_gmres_trajectory_medium.jl`，固定 checkpoint `100 / 300 / 600`，直接比较 dense strong 与 MLFMA strong 的解轨迹
    - 当前结果已落盘到 `test_results/accuracy/PMCHW_gmres_trajectory_medium.csv`
    - 当前结论：
        1. `default` 路径在 `300` 步时与 dense strong 几乎重合（`gap = 0.00620Ω / 0.00566Ω`），到 `600` 步仍仅 `0.03815Ω / 0.00236Ω`
        2. `loose_near` 到 `600` 步的解向量偏差增大到 `1.52e-2`，但阻抗 gap 仍只到 `0.48107Ω / 0.02029Ω`
        3. 因此 receive 修复后，long-Krylov 主线不再表现为“大尺度阻抗分叉”，而更像是深层轨迹上的次级放大问题
- [x] 15.G21 PMCHW vs N-Muller medium dense GMRES 轨迹对照
    - 已新增专门测试 `test/test_nmuller_gmres_trajectory_medium.jl`，当前结果 `25/25 pass`
    - 已新增 benchmark `benchmark/compare_pmchw_nmuller_gmres_trajectory_medium.jl`
    - 当前结果已落盘到 `test_results/accuracy/PMCHW_NMuller_gmres_trajectory_medium.csv`
    - 当前结论：
        1. 在 `50 / 100 / 150 / 200 / 250` 全部 checkpoint 上，`NMuller` 的 `rel_vs_LU` 与 `resnorm` 都持续优于 `PMCHW`
        2. `PMCHW` 在 `250` 步时仍为 `rel=1.000066e+00`, `res=7.078038e-01`，而 `NMuller` 已收敛到 `rel=7.808047e-02`, `res=4.174158e-02`
        3. 因此当前剩余边界更像是 `PMCHW formulation / conditioning` 问题，而不是 N-Muller 支线本身或已修复的 PMCHW fast backend 主干
- [x] 15.G22 PMCHW vs N-Muller medium plane-wave dense GMRES 轨迹对照
    - 已新增专门测试 `test/test_nmuller_planewave_gmres_trajectory_medium.jl`，当前结果 `25/25 pass`
    - 已新增 benchmark `benchmark/compare_pmchw_nmuller_planewave_gmres_trajectory_medium.jl`
    - 当前结果已落盘到 `test_results/accuracy/PMCHW_NMuller_planewave_gmres_trajectory_medium.csv`
    - 当前结论：
        1. 把 RHS 从 random probe 切到 `PlaneWave` 物理激励后，`NMuller` 仍在全部 checkpoint 上同时优于 `PMCHW` 的 `rel_vs_LU` 与相对残差
        2. `PMCHW` 在 `250` 步时仍为 `rel=1.002186e+00`, `rres=1.779276e-03`，而 `NMuller` 已在 `67` 步内停机并稳定到 `rel=1.149702e-02`, `rres=9.583083e-06`
        3. 因而“剩余问题更像 PMCHW formulation / conditioning，而不是 repaired backend”这一判断已经从随机 RHS 扩展到物理激励 RHS
- [x] 15.G23 PMCHW medium plane-wave dense weak/strong 轨迹对照
    - 已新增专门测试 `test/test_pmchw_gate_s_planewave_trajectory_medium.jl`，当前结果 `31/31 pass`
    - 已新增 benchmark `benchmark/compare_pmchw_gate_s_planewave_trajectory_medium.jl`
    - 当前结果已落盘到 `test_results/accuracy/PMCHW_gate_s_planewave_trajectory_medium.csv`
    - 当前结论：
        1. 在 `PlaneWave` 物理激励下，`strong` 相比 `weak` 的 `rel_vs_LU` 只带来边际改善，且 5 个 checkpoint 上改善比例都不足 1%
        2. 到 `250` 步时，`weak = rel 1.002186e+00 / rres 1.779276e-03`，`strong = rel 1.001109e+00 / rres 1.868816e-03`，两者都没有收回到接近 LU 的区域
        3. 因而 PMCHW 当前剩余问题不能再简单归结为“只差 strong-form”；在 plane-wave 工况下，主边界仍落在更深的 formulation / conditioning 机制
- [x] 15.G24 PMCHW vs N-Muller medium plane-wave dense restart 扫描
    - 已新增专门测试 `test/test_pmchw_nmuller_planewave_restart_sweep_medium.jl`，当前结果 `22/22 pass`
    - 已新增 benchmark `benchmark/compare_pmchw_nmuller_planewave_restart_sweep_medium.jl`
    - 当前结果已落盘到 `test_results/accuracy/PMCHW_NMuller_planewave_restart_sweep_medium.csv`
    - 当前结论：
        1. `IterativeSolvers.gmres` 默认 `restart = min(20, size(A,2))`；此前 plane-wave dense 轨迹里 `PMCHW` 的坏终点确有明显 restart 主导成分
        2. 对同一 medium 夹具，`PMCHW` 从 `restart=20` 提升到 `restart=250` 后，`rel_vs_LU` 从 `1.002186e+00` 降到 `5.186365e-02`，相对残差从 `1.779276e-03` 降到 `8.501219e-06`
        3. 但在同一 `restart=250` 下，`NMuller` 仍保持 `rel=7.800091e-03`, `rres=8.610863e-06`, `iters=44`；因此新边界应修正为“PMCHW 先有强 restart 敏感性，剥离后仍保留明显 formulation gap”，而不是单纯把全部剩余问题都归为 backend 或 strong-form

### 当前主线说明

- **主交付 formulation 仍是 PMCHW**：Phase 15 当前所有 medium 级精度与收敛治理都以 PMCHW shell + MLFMA backend 为主线。
- **N-Muller 当前定位为 dense 对照基线**：只用于比较 conditioning / GMRES 行为，帮助区分“formulation 问题”与“backend 问题”；它不是当前阻塞主线。
- **当前剩余问题的最新归因**：已从早期的 `M` pass receive 错误、再到中期的“budget 主导 long-Krylov 大分叉”，收缩为 **更深 Krylov 轨迹 / formulation 条件数放大**。
- **最新证据边界**：上述归因不仅在 random-RHS dense 对照下成立，在 `PlaneWave` 物理激励 RHS 下也继续成立。
- **对 PMCHW 自身的新分界**：`PlaneWave` 工况下的 dense `weak/strong` 对照表明，strong-form 不是当前 medium 主问题的单独解；它只带来边际改善，无法把 PMCHW 拉回 LU 邻域。
- **对 restart 机制的新分界**：默认 `restart=20` 会显著夸大 PMCHW plane-wave dense 轨迹失真；但把 `restart` 提升到 `250` 后，PMCHW 虽已把相对残差收回到 `~1e-5`，相对 LU 误差仍显著高于 `NMuller`。
- **本段收尾结论**：Phase 15 当前这段 dense plane-wave 归因子流可先收口到 `G22/G23/G24`；后续若继续推进，应直接进入显式 full-restart / Arnoldi 级诊断，而不再重复 weak-vs-strong 或默认-restart 现象复测。
- **下一执行顺序**：
    1. 继续在 PMCHW 主线上做更深 Krylov / restart / Arnoldi 诊断，并把默认 restart 与 full-restart 行为明确拆开；
    2. 用 N-Muller 只做同夹具 dense 轨迹对照；
    3. 暂不扩展 N-Muller 天线语义，不让支线反客为主。

- **工作区清理状态（2026-03-08）**：已删除 `scripts/` 下 5 个一次性诊断脚本（4 个 `tmp_*` 探针 + 1 个 `diagnose_pmchw_farfield_blocks.jl`），正式测试、benchmark、CSV 结果与计划文档全部保留；下一阶段若需要继续追踪，应直接在现有专门回归入口上扩展，而不是重新堆积临时脚本。

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
