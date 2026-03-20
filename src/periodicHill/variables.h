#ifndef PERIODIC_HILL_VARIABLES_H
#define PERIODIC_HILL_VARIABLES_H

// ================================================================
// Periodic Hill Flow — Parameter Control (FluidX3D Edition)
// 對齊 GILBM D3Q19 reference code 的 variables.h 結構
// ================================================================

// ================================================================
// 1. 數學常數
// ================================================================
#define PH_PI   3.14159265358979323846264338327950
#define PH_CS   (1.0/1.732050807568877)  // 1/√3, LBM 聲速

// ================================================================
// 2. 物理域幾何 (Hill height H=1 為參考長度)
// ================================================================
#define PH_LX       (4.5)      // 展向 (spanwise) 長度 / H
#define PH_LY       (9.0)      // 流向 (streamwise) 長度 / H = hill-to-hill 週期長度
#define PH_LZ       (3.036)    // 法向 (wall-normal) 長度 / H
#define PH_H_HILL   (1.0)      // hill 高度 (Re_h 參考長度)

// ================================================================
// 3. 網格設定
// ================================================================
#define PH_NX       32          // 展向格點數 (spanwise, periodic)
#define PH_H        32          // hill 高度的格點數 (resolution parameter)
                                // 增加 PH_H 可提高精度; PH_NY, PH_NZ 自動跟隨
#define PH_NY       (9 * PH_H)                                     // 流向格點數 = LY/H * h
#define PH_NZ       ((uint)(3.036f * (float)PH_H) + 2u)            // 法向格點數 = LZ/H * h + 2(BB層)

// --- 多 GPU 域分解 (Multi-GPU Domain Decomposition) ---
// Dx×Dy×Dz = 總 GPU 數量
// 建議沿最長方向 (Y, 流向) 切分; NY 必須能被 DY 整除
//   單 GPU: DX=1, DY=1, DZ=1  →  1 GPU
//   4  GPU: DX=1, DY=4, DZ=1  →  NY/DY = 288/4 = 72 per GPU
//   8  GPU: DX=1, DY=8, DZ=1  →  NY/DY = 288/8 = 36 per GPU
#define PH_DX       1u
#define PH_DY       4u          // ← 改這裡: 1=單GPU, 4=四GPU, 8=八GPU
#define PH_DZ       1u

// ================================================================
// 4. 物理參數
// ================================================================
#define PH_Re       75          // Reynolds number (基於 H_HILL 和 Uref)
                                // 可選: 75, 1400, 5600, 10595
#define PH_Uref     0.05f       // 參考速度 = bulk velocity (LBM units)
                                // Ma = Uref/cs ≈ 0.087 (足夠低)
                                // GILBM ref: Re75=0.0503, Re1400=0.0776
#define PH_niu      ((float)PH_Uref * (float)PH_H / (float)PH_Re)  // 運動黏度 ν = Uref*H/Re

// ================================================================
// 5. 驅動力 (volume force in streamwise direction)
// ================================================================
// h_eff = (LZ - H_HILL) * h = 有效通道高度 (hill crest 以上)
// Poiseuille analogy: fy = 8*ν*Ub / h_eff²
#define PH_H_EFF    ((PH_LZ - PH_H_HILL) * (float)PH_H)
#define PH_FY       (8.0f * PH_niu * PH_Uref / (PH_H_EFF * PH_H_EFF))

// ================================================================
// 6. 碰撞算子設定 (FluidX3D 由 defines.hpp 控制)
// ================================================================
// FluidX3D: SRT (defines.hpp 中 #define SRT)
// 鬆弛時間: τ = 3ν + 0.5

// ================================================================
// 7. 模擬控制 — 輸出頻率
// ================================================================
#define PH_NDTMIT       50      // 每 N 步輸出 monitor (checkrho, Ustar_Force_record)
#define PH_NDTVTK       1000    // 每 N 步輸出 VTK
#define PH_NDTCKPT      10000   // 每 N 步輸出 binary checkpoint

// ================================================================
// 8. FTT 閾值與統計控制
// ================================================================
// Flow-Through Time: 1 FTT = NY / Uref lattice time steps
//
// Stage 0: FTT < FTT_STATS_START → 只跑瞬時場，不累積統計量
//                                   checkpoint 只寫 rho + u (分佈函數近似)
// Stage 1: FTT ≥ FTT_STATS_START → 累積統計量 (mean, RS)
//                                   VTK 寫入平均統計量
//                                   checkpoint 包含統計量累積和
#define PH_FTT_STATS_START  50.0    // 統計量開始累積
#define PH_FTT_STOP         100.0   // 模擬結束

// ================================================================
// 9. VTK 輸出格式 (binary, 對齊 FluidX3D 原生格式)
// ================================================================
// 0 = 基本: VECTORS velocity + SCALARS u_inst/v_inst/w_inst
//           + (Stage 1) U_mean/W_mean/uu_RS/uw_RS/ww_RS/k_TKE
// 1 = 完整: Level 0 + V_mean + vv_RS
#define PH_VTK_OUTPUT_LEVEL  0

// ================================================================
// 10. 重啟 (Restart) 配置
// ================================================================
// 0 = 冷啟動 (zero velocity → u_bulk 初始化, ρ=1)
// 1 = 從 binary checkpoint 續跑 (讀取 rho + u + 統計量)
#define PH_INIT             0

// INIT=1 用: binary checkpoint 目錄路徑
#define PH_RESTART_DIR      "checkpoint/step_100000"

// ================================================================
// 11. 初始擾動 (觸發 3D 湍流轉捩)
// ================================================================
// Re=75 為層流，不需擾動; Re>1400 建議開啟
#define PH_PERTURB_INIT     0       // 1=注入隨機擾動, 0=不擾動
#define PH_PERTURB_PERCENT  5       // 擾動振幅 (% of Uref), 典型 1-10%

// ================================================================
// 12. 控制台輸出頻率
// ================================================================
#define PH_PRINT_INTERVAL   500     // 每 N 步印出進度 (step, FTT, Ub, Ma)

#endif // PERIODIC_HILL_VARIABLES_H

/*
備註: 座標系對應關係 (與 GILBM reference 完全一致)

  FluidX3D 方向    物理方向        Benchmark 符號 (ERCOFTAC)
  ──────────────────────────────────────────────────────────
  x (i)           展向 spanwise    V (benchmark)
  y (j)           流向 streamwise  U (benchmark)
  z (k)           法向 wall-normal W (benchmark)

  VTK 輸出映射:
    VTK VECTORS velocity  = (u.y, u.x, u.z) / Uref
    VTK u_inst            = u.y / Uref    (streamwise)
    VTK v_inst            = u.x / Uref    (spanwise)
    VTK w_inst            = u.z / Uref    (wall-normal)
    VTK U_mean            = sum_uy / N / Uref
    VTK W_mean            = sum_uz / N / Uref
    VTK uu_RS             = (<uy²> - <uy>²) / Uref²
    VTK uw_RS             = (<uy·uz> - <uy><uz>) / Uref²
    VTK ww_RS             = (<uz²> - <uz>²) / Uref²
    VTK k_TKE             = 0.5(uu + vv + ww) / Uref²
*/
