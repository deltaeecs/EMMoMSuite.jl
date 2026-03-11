"""
    Accuracy — Phase 14 精度对比工具模块

包含：
- `FekoReader`:      解析 MoM_AllinOne compare_feko/*.csv 文件
- `ReferenceData`:   Mie PEC/介质 + 偶极子解析参考
- `AccuracyMetrics`: AccuracyResult / AntennaAccuracyResult 指标计算
"""
module Accuracy

include("FekoReader.jl")
include("ReferenceData.jl")
include("AccuracyMetrics.jl")
include("BenchmarkReportData.jl")

using .FekoReader: read_feko_rcs, split_phi_cuts
using .ReferenceData: mie_pec_rcs_dBsm, mie_pec_bistatic_rcs_dBsm,
                      mie_dielectric_rcs_dBsm, mie_dielectric_bistatic_rcs_dBsm,
             dipole_halfwave_Zin_analytic, dipole_resonant_Zin_analytic,
             dipole_halfwave_farfield_analytic, extract_sphere_radius
using .AccuracyMetrics: AccuracyResult, AntennaAccuracyResult,
                        compute_rcs_accuracy, compute_antenna_accuracy,
                        print_accuracy_report

export read_feko_rcs, split_phi_cuts
export mie_pec_rcs_dBsm, mie_pec_bistatic_rcs_dBsm,
       mie_dielectric_rcs_dBsm, mie_dielectric_bistatic_rcs_dBsm,
       dipole_halfwave_Zin_analytic, dipole_resonant_Zin_analytic,
    dipole_halfwave_farfield_analytic, extract_sphere_radius
export AccuracyResult, AntennaAccuracyResult,
       compute_rcs_accuracy, compute_antenna_accuracy,
       print_accuracy_report
export AccuracyCurve, AccuracyCurveSummary, PerformanceBenchmarkResult,
       load_accuracy_curve, summarize_accuracy_curve,
       accuracy_curve_group, accuracy_curve_cut,
       load_performance_results

end # module Accuracy
