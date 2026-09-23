"""
daylight_basis_check.py
=======================
驗證「日光光譜只需要 ~3 個基底」

為什麼不能直接拿 CIE D 系列光譜做 PCA？
    CIE 日光模型本身就是 S0 + M1*S1 + M2*S2 三個基底疊出來的，拿它做 PCA 得到 3 是循環論證。
    所以這裡用一個「與 CIE 無關」的簡化物理大氣模型（類 SPCTRAL2）產生日光光譜，
    輸入是 7 個獨立的物理參數（太陽天頂角、氣溶膠濁度與波長指數、單次散射反照率、
    臭氧量、可降水量、直射可見比例、雲量），再看 PCA 到底需要幾個成分。

驗證項目：
    1. 累積變異量 vs 基底數
    2. 光譜重建誤差 (RMSE) 與色度誤差 (delta u'v') vs 基底數
    3. PCA 前兩個變異方向與 CIE S1/S2 的子空間夾角
    4. 每個 PC 主要由哪個物理參數驅動

若有實測日光資料庫（例如 Granada 日光資料庫、CSV：第一列波長、每列一條光譜），
呼叫 analyze(spectra, wl) 即可用同一套流程分析。

依賴：numpy, matplotlib, colour-science (pip install colour-science)
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

WL = np.arange(380, 781, 5.0)          # nm，5 nm 取樣
I560 = int(np.where(WL == 560)[0][0])

# ----------------------------------------------------------------------------
# 1. 物理日光模擬器（簡化版 SPCTRAL2，所有項共用的 cos(zenith) 已因正規化省略）
# ----------------------------------------------------------------------------

def planck(wl_nm, T):
    lam = wl_nm * 1e-9
    h, c, k = 6.626e-34, 2.998e8, 1.381e-23
    return (2 * h * c**2 / lam**5) / (np.exp(h * c / (lam * k * T)) - 1)

E0 = planck(WL, 5778.0)
E0 = E0 / E0[I560]

def air_mass(zenith_deg):
    z = np.radians(zenith_deg)
    return 1.0 / (np.cos(z) + 0.50572 * (96.07995 - zenith_deg) ** -1.6364)   # Kasten & Young 1989

lam_um = WL / 1000.0
TAU_R = 0.008735 * lam_um ** -4.08                               # Rayleigh
K_O3  = 0.13 * np.exp(-0.5 * ((WL - 600) / 55) ** 2)             # 臭氧 Chappuis 帶（近似）
K_H2O = 0.017 * np.exp(-0.5 * ((WL - 720) / 8) ** 2)             # 720 nm 水氣弱吸收
K_O2  = 0.10 * np.exp(-0.5 * ((WL - 760) / 4) ** 2)              # 760 nm 氧氣 A 帶

def simulate_daylight(n, seed=0):
    rng = np.random.default_rng(seed)
    S, P = [], []
    for _ in range(n):
        zen   = rng.uniform(0, 85)                                # 太陽天頂角
        beta  = np.exp(rng.uniform(np.log(0.02), np.log(0.5)))    # Angstrom 濁度
        alpha = rng.uniform(0.5, 2.5)                             # Angstrom 波長指數
        w0    = rng.uniform(0.80, 0.98)                           # 氣溶膠單次散射反照率
        u_o3  = rng.uniform(0.20, 0.45)                           # 臭氧 atm-cm
        pw    = np.exp(rng.uniform(np.log(0.3), np.log(5.0)))     # 可降水 cm
        r = rng.uniform()                                         # 直射可見比例
        sun   = 1.0 if r < 0.6 else (0.0 if r < 0.75 else rng.uniform())
        cloud = 0.0 if rng.uniform() < 0.4 else rng.beta(0.5, 0.5)

        m = air_mass(zen)
        tau_a = beta * lam_um ** -alpha
        T_R  = np.exp(-m * TAU_R)
        T_as = np.exp(-m * w0 * tau_a)
        T_aa = np.exp(-m * (1 - w0) * tau_a)
        T_abs = np.exp(-m * (u_o3 * K_O3 + pw * K_H2O + K_O2))
        Fs = 0.65 + 0.2 * np.cos(np.radians(zen))                 # 氣溶膠前向散射比例

        direct = E0 * T_R * T_as * T_aa * T_abs
        sky_R  = E0 * T_abs * T_aa * (1 - T_R ** 0.95) * 0.5
        sky_a  = E0 * T_abs * T_aa * T_R ** 1.5 * (1 - T_as) * Fs
        clear  = sun * direct + sky_R + sky_a
        overcast = E0 * T_abs * T_R ** 0.4 * 0.3                  # 雲層近乎光譜中性
        s = (1 - cloud) * clear + cloud * overcast
        S.append(s / s[I560])
        P.append([zen, beta, alpha, w0, u_o3, pw, sun, cloud])
    names = ["zenith", "beta", "alpha", "w0", "ozone", "water", "sun", "cloud"]
    return np.array(S), np.array(P), names

# ----------------------------------------------------------------------------
# 2. CIE 資料（S0/S1/S2、色匹配函數）
# ----------------------------------------------------------------------------

def load_cie(wl):
    import colour
    shape = colour.SpectralShape(wl[0], wl[-1], wl[1] - wl[0])
    B = colour.colorimetry.SDS_BASIS_FUNCTIONS_CIE_ILLUMINANT_D_SERIES
    S = np.stack([B[k].copy().align(shape).values for k in ("S0", "S1", "S2")], axis=1)
    cmfs = colour.MSDS_CMFS["CIE 1931 2 Degree Standard Observer"].copy().align(shape).values
    return S, cmfs

def to_uv(spectra, cmfs):
    XYZ = spectra @ cmfs
    d = XYZ[:, 0] + 15 * XYZ[:, 1] + 3 * XYZ[:, 2]
    return np.stack([4 * XYZ[:, 0] / d, 9 * XYZ[:, 1] / d], axis=1)

# ----------------------------------------------------------------------------
# 3. 分析
# ----------------------------------------------------------------------------

def principal_angles(A, B):
    """兩個子空間（欄向量張成）之間的主角度，單位：度"""
    Qa, _ = np.linalg.qr(A)
    Qb, _ = np.linalg.qr(B)
    s = np.linalg.svd(Qa.T @ Qb, compute_uv=False)
    return np.degrees(np.arccos(np.clip(s, -1, 1)))

def analyze(spectra, wl, params=None, param_names=None, kmax=8, tag="sim"):
    n, B = spectra.shape
    mu = spectra.mean(axis=0)
    Xc = spectra - mu
    U, s, Vt = np.linalg.svd(Xc, full_matrices=False)
    evr = s**2 / np.sum(s**2)
    cum = np.cumsum(evr)

    S_cie, cmfs = load_cie(wl)
    uv_true = to_uv(spectra, cmfs)

    print(f"\n=== {tag}: {n} 條光譜, {B} 個波段 ===")
    print(f"{'k':>2} {'單一變異%':>10} {'累積變異%':>10} {'RMSE':>9} {'RMSE最大':>9} {'Δu\'v\'中位':>11} {'Δu\'v\'95%':>10}")
    rmse_k, duv_k = [], []
    for k in range(1, kmax + 1):
        rec = mu + (Xc @ Vt[:k].T) @ Vt[:k]
        err = rec - spectra
        rmse = np.sqrt(np.mean(err**2, axis=1))
        duv = np.linalg.norm(to_uv(rec, cmfs) - uv_true, axis=1)
        rmse_k.append(rmse.mean()); duv_k.append(np.median(duv))
        print(f"{k:>2} {100*evr[k-1]:>10.3f} {100*cum[k-1]:>10.3f} {rmse.mean():>9.4f} {rmse.max():>9.4f} "
              f"{np.median(duv):>11.5f} {np.percentile(duv, 95):>10.5f}")
    print("（光譜已在 560 nm 正規化為 1，RMSE 為相對值；Δu'v' 約 0.002 為色度可辨閾值）")

    # 與 CIE 基底比較：變異方向 span{PC1,PC2} vs span{S1,S2}
    ang2 = principal_angles(Vt[:2].T, S_cie[:, 1:3])
    ang3 = principal_angles(np.column_stack([mu, Vt[:2].T]), S_cie)
    print(f"\nspan{{PC1,PC2}} 與 span{{S1,S2}} 的主角度: {np.round(ang2, 1)} 度")
    print(f"span{{mean,PC1,PC2}} 與 span{{S0,S1,S2}} 的主角度: {np.round(ang3, 1)} 度")

    # 去均值後的變異用 CIE S1/S2 兩個方向擬合，與 PCA k=2 的殘差比較
    # （直接用 S0/S1/S2 擬合原始光譜會被「Planck 太陽 vs 真實太陽」的平均形狀差異主導，不公平）
    coef, *_ = np.linalg.lstsq(S_cie[:, 1:3], Xc.T, rcond=None)
    rec_cie = mu + (S_cie[:, 1:3] @ coef).T
    rmse_cie = np.sqrt(np.mean((rec_cie - spectra)**2, axis=1))
    duv_cie = np.linalg.norm(to_uv(rec_cie, cmfs) - uv_true, axis=1)
    print(f"去均值變異用 CIE S1/S2 擬合: RMSE 平均 {rmse_cie.mean():.4f} (PCA k=2: {rmse_k[1]:.4f}), "
          f"Δu'v' 中位 {np.median(duv_cie):.5f} (PCA k=2: {duv_k[1]:.5f})")

    # 各 PC 分數與物理參數的相關
    if params is not None:
        scores = Xc @ Vt[:4].T
        print("\n各 PC 分數與物理參數的相關係數 (|r| >= 0.3 才顯示):")
        for i in range(4):
            r = [np.corrcoef(scores[:, i], params[:, j])[0, 1] for j in range(params.shape[1])]
            shown = ", ".join(f"{param_names[j]}={r[j]:+.2f}" for j in np.argsort(-np.abs(r)) if abs(r[j]) >= 0.3)
            print(f"  PC{i+1} ({100*evr[i]:.1f}%): {shown}")

    return dict(mu=mu, Vt=Vt, evr=evr, cum=cum, rmse_k=rmse_k, duv_k=duv_k, S_cie=S_cie)

# ----------------------------------------------------------------------------
# 4. 繪圖
# ----------------------------------------------------------------------------

def make_figure(spectra, wl, res, path):
    fig, ax = plt.subplots(2, 2, figsize=(12, 9))

    idx = np.random.default_rng(1).choice(len(spectra), 60, replace=False)
    ax[0, 0].plot(wl, spectra[idx].T, color="gray", alpha=0.35, lw=0.8)
    ax[0, 0].plot(wl, res["mu"], color="k", lw=2, label="mean")
    ax[0, 0].set(title="Simulated daylight spectra (normalized at 560 nm)", xlabel="nm", ylabel="relative SPD")
    ax[0, 0].legend()

    k = np.arange(1, 9)
    ax[0, 1].bar(k, 100 * res["evr"][:8], color="steelblue", label="individual")
    ax[0, 1].plot(k, 100 * res["cum"][:8], "o-", color="darkred", label="cumulative")
    ax[0, 1].axhline(99, ls="--", color="gray", lw=0.8); ax[0, 1].axhline(99.9, ls=":", color="gray", lw=0.8)
    ax[0, 1].set(title="Explained variance", xlabel="number of components", ylabel="%", ylim=(0, 102))
    ax[0, 1].legend()

    S = res["S_cie"]
    for i, c in enumerate(["C0", "C1", "C2"]):
        pc = res["Vt"][i] if i < 2 else res["Vt"][2]
        ax[1, 0].plot(wl, pc / np.abs(pc).max(), color=c, lw=2, label=f"PC{i+1}")
    for i, (c, name) in enumerate(zip(["C0", "C1", "C2"], ["S1", "S2"])):
        sc = S[:, i + 1]
        sign = np.sign(np.dot(sc, res["Vt"][i]))
        ax[1, 0].plot(wl, sign * sc / np.abs(sc).max(), "--", color=c, lw=1.5, label=f"CIE {name}")
    ax[1, 0].axhline(0, color="gray", lw=0.5)
    ax[1, 0].set(title="PCA variation directions vs CIE S1/S2 (scaled, sign-aligned)", xlabel="nm")
    ax[1, 0].legend(ncol=2, fontsize=8)

    ax2 = ax[1, 1]
    ax2.plot(k, res["rmse_k"], "s-", color="steelblue", label="spectral RMSE (mean)")
    ax2.set(title="Reconstruction error vs number of components", xlabel="number of components", ylabel="RMSE")
    ax2.set_yscale("log")
    ax3 = ax2.twinx()
    ax3.plot(k, res["duv_k"], "o-", color="darkred", label="Δu'v' (median)")
    ax3.axhline(0.002, ls="--", color="darkred", lw=0.8)
    ax3.set_ylabel("Δu'v'"); ax3.set_yscale("log")
    h1, l1 = ax2.get_legend_handles_labels(); h2, l2 = ax3.get_legend_handles_labels()
    ax2.legend(h1 + h2, l1 + l2, loc="upper right")

    fig.tight_layout()
    fig.savefig(path, dpi=130)
    print(f"\n圖已存至 {path}")

# ----------------------------------------------------------------------------

if __name__ == "__main__":
    spectra, params, names = simulate_daylight(3000, seed=0)
    res = analyze(spectra, WL, params, names, kmax=8, tag="physics-simulated daylight")
    make_figure(spectra, WL, res, "daylight_basis_check.png")
