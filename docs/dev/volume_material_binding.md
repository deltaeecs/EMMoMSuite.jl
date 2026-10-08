# 体积分四面体材料机制设计（S3）

状态：设计稿（未实现）。前置：Phase 19.5 `BoundMesh`（src\Geometry\MeshMaterialBind.jl）、
Phase 19.1 gmsh 四面体网格（src\Geometry\TetMeshing.jl）、GmshIO .msh v4.1 读取、
benchmark\cases 的界面材料模型（`Interface` / `MaterialContext` / 工厂分派）。

## 1. 目标

把 benchmark\cases 已定稿的「界面两侧材料」模型（`Interface(surface, plus, minus,
closed, flip)`）推广到**体积分（VIE / 体等效原理）**路径：几何不再是单张闭曲面，
而是被 `Physical Volume` 划分成若干材料区域的四面体网格，每个四面体通过
`mesh.tags` 绑定一种材料。设计原则与 S1/S2 一致：

1. 材料描述与求解器分派解耦：唯一下游分支点是工厂函数；
2. 复用既有地基，不新造绑定机制 —— `BoundMesh{MeshT, ModelT}` 已经实现
   tag→材料绑定（`bind_materials` / `validate_bindings` / `element_material`），
   本设计只补「命名区域 → tag」与「多区域」两层；
3. 行为可验证：与解析参照（介质球 Mie 体散射）或 PMCHW 表面解做交叉对比。

## 2. 数据链路

```
.geo (Physical Volume per region)
  → gmsh 3D mesh (TetrahedraMesh, tet_tags = physical volume id, GmshIO.jl:181)
  → bind_materials(mesh, Dict(tag => material))          # BoundMesh 已有
  → VolumeMaterialContext(k0, eta0, regions, layout)     # 新：工厂产出
  → SWG basis + 体积分算子装配 → LU/GMRES
  → 体等效流 J_v = jω(ε-ε₀)E（后处理按 region 派发）
  → RCS / 远场（复用 radarCrossSection 既有路径）
```

### 2.1 几何与 tag 约定

- .geo 中每个材料区域显式声明 `Physical Volume("region_name") = {ids};`；
  无 Physical Volume 的四面体 tag 默认 0（`TetrahedraMesh` 构造器，
  src\Geometry\MeshTypes.jl:39-41），必须由 `validate_bindings` 拦截。
- 传导/介质体外表面（供表面-体积耦合或边界条件）用 MeshBoundary.jl
  `extract_boundary` 从四面体网格取闭曲面，继承相邻四面体的 region 归属，
  避免 .geo 里重复维护 Physical Surface。
- 每个四面体的材料唯一：`element_material(bm, idx)`（KeyError 即报错），
  不允许重叠区域 —— 由 gmsh 保证，读入后校验 tag 唯一性。

### 2.2 材料模型

复用 benchmark\cases\CaseRunner.jl 的 `Material` 抽象：
`PEC()`（体内 PEC 区域 → 退化为表面 EFIE，不进体积分）、
`Dielectric(eps_r, mu_r)`（→ 体积分区域）、`AIR` 常量（背景介质）。
cases.toml 侧沿用 S2 的 `[material.<name>]` 注册表，新增
`[[case.region]]` 块（见 §4）。映射到求解器材料模型用现有的
`Isotropic(eps_r*eps0*(1-0im), mu_r)`（Phase 19.5 示例用法）。

## 3. 求解器接口（新 VolumeMaterialContext + 工厂）

```julia
Base.@kwdef struct VolumeMaterialContext
    k0::Float64; eta0::Float64
    bm::BoundMesh                      # TetrahedraMesh + tag→Material
    layout::Symbol                     # :j_only | :jm | :volume
    interfaces::Vector{Interface}      # 可选：并存的表面界面
end

build_volume_operator(bm, freq; ie = "auto", alpha = 0.5) -> (op, ctx)
```

`ie = "auto"` 规则（纯材料驱动，与 S1 工厂同构）：
- 全部区域 PEC → 表面 EFIE（抽取外表面，走现有路径）；
- 恰好一个介质区域且用户显式 `ie="PMCHW"` → 现有 PMCHW 表面路径（回归保护）；
- 存在介质区域 → SWG 体积分算子（`:volume`）。

后处理派发沿 `postprocess_rcs` / `postprocess_farfield` 的模式扩展
`:volume` 分支：体等效电流密度逐四面体积分后再进远场/RCS 求积。

## 4. cases.toml 扩展（草案）

```toml
[material.fr4]
eps_r = 4.4
eps_r_im = -0.088        # tanδ = 0.02

[[case]]
name = "V1_dielectric_box_VIE"
geo  = "box_1m.geo"
dim  = 3                  # 体四面体
freq = 3.0e8
mie_radius = 0.5          # dielectric box 无解析解 → 用 Mie 球仅作趋势参照
# 或对 dielectric sphere 用 mie_dielectric_bistatic_rcs_dBsm 精确对比

[[case.region]]
surface = "body"          # Physical Volume 标签
material = "fr4"
```

解析规则：`region` 列表替换（而非叠加）`interface` 列表；两者互斥，
同时出现报错。

## 5. 需要实现/补齐的部件

| 部件 | 现状 | 工作量 |
|---|---|---|
| Physical Volume → tet tags | GmshIO.jl:181 已读 tet_tags；需确认 v4.1 entities 映射（ Physical Volume id 还是体实体 tag）并加单测 | 小 |
| `bind_materials` 命名区域入口 | 现仅 `Dict{Int}`；加 `bind_regions(mesh, Dict{String})`：读 .msh 时顺带返回 region 名称→tag 映射 | 小 |
| SWG 基函数 | MeshTypes.jl:223 有 SWG 符号存取痕迹，需确认完整 SWGBasis 是否存在；若无则实现（低阶：每四面体 1 个常量磁流偶极子基） | 中 |
| 体积分算子装配 | 无（主要缺口）。非均匀介质偶极/SWG 装配 + 奇异积分处理 | 大 |
| 体等效流后处理 | radarCrossSection/farField 接受体流分布的变体 | 中 |
| 交叉验证 | 介质球 VIE vs `mie_dielectric_bistatic_rcs_dBsm`（src\Accuracy\ReferenceData.jl:180）；同网格收敛阶报告 | 小 |

## 6. 验收标准

1. `V1_dielectric_sphere_VIE`（r=0.15 m，εᵣ=4，600 MHz，mesh_size=0.02 与
   G10 同密）：RCS RMSE ≤ 0.5 dB（对齐 G10 的 0.044 dB 量级留裕度）；
2. 全 PEC 区域用例经 auto 规则回落到现有 EFIE 路径，rcs.csv 与 G1 基线逐位一致；
3. `validate_bindings` 对缺 tag 网格报错而非静默放行；
4. cases.toml 旧格式（S2 interface 写法）行为不变。
