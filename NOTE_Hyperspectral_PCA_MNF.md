# 高光譜資料的 SVD / PCA / MNF 筆記

---

## 1. 為什麼要對高光譜做 SVD / PCA

### 資料特性
- 每個像素有數十到數百個波段
- 相鄰波段的反射率高度相關，資訊大量重複
- 有效資訊通常只集中在前幾個到十幾個主成分

### 目的

| 目的 | 說明 |
|---|---|
| 降維與去相關 | 旋轉到互不相關的座標，前幾個 PC 保留 95–99% 變異量 |
| 去雜訊（低秩近似） | $X = U\Sigma V^T$，只保留前 $k$ 個奇異值重建，小奇異值的成分多為雜訊 |
| 避免 Hughes 現象 | 特徵維度高、樣本少會使分類精度下降；先降維再訓練分類器 |
| 估計本質維度 / 端元數 | 奇異值衰減曲線（scree plot）可估訊號子空間維度，是 VCA、N-FINDR 的前置步驟 |
| 視覺化 | 前三個 PC 當 RGB 顯示 |
| 異常 / 目標偵測 | 背景建模為低秩子空間，投影殘差大的像素為異常（RX detector 的基礎） |

### SVD 與 PCA 的關係
- PCA = 對「去均值後的資料矩陣」做 SVD
- 右奇異向量 $V$ = 共變異矩陣的特徵向量（主成分方向）
- $\sigma_i^2 / (n-1)$ = 各主成分的變異量
- 實務上直接做 SVD 比先算共變異矩陣再特徵分解數值更穩定

### 限制
- **只看變異量，不看訊噪比**：雜訊大的波段可能排到前面 → 改用 MNF
- 主成分是線性組合，失去物理光譜意義
- 小面積目標的訊息可能藏在低變異量成分，降維太狠會丟失
- 全域線性轉換

---

## 2. MNF（Minimum Noise Fraction）

### 動機
PCA 按變異量排序，MNF 按 **SNR** 排序：前面的成分最乾淨，後面的成分是純雜訊。（Green et al., 1988）

### 數學定義
資料 = 訊號 + 雜訊：$X = S + N$，$\Sigma_X = \Sigma_S + \Sigma_N$

找投影方向 $w$ 使雜訊比例最小：

$$\text{NF}(w) = \frac{w^T \Sigma_N w}{w^T \Sigma_X w} \quad \rightarrow \quad \min$$

等價於廣義特徵值問題：

$$\Sigma_X w = \lambda \Sigma_N w, \qquad \lambda_i = 1 + \text{SNR}_i$$

- 特徵值由大到小排序
- $\lambda \approx 1$ 的成分幾乎全是雜訊

### 兩步驟實作（Noise-Adjusted PCA, Lee et al., 1990）

**步驟一：雜訊白化**
1. $\Sigma_N = E \Lambda E^T$
2. 白化矩陣 $F = E \Lambda^{-1/2}$
3. $Y = F^T X$ → 雜訊共變異矩陣變成單位矩陣

**步驟二：對白化後資料做 PCA**
- 對 $\Sigma_Y = F^T \Sigma_X F$ 做特徵分解，得 $G$
- 因雜訊已壓平為 1，變異量排序 = SNR 排序

最終轉換矩陣：$W = F\,G$

### 雜訊共變異矩陣的估計（影像資料）

| 方法 | 做法 | 備註 |
|---|---|---|
| 相鄰像素差分（shift difference） | $\Sigma_N \approx \tfrac{1}{2}\text{Cov}(X(i,j) - X(i+1,j))$ | ENVI 預設；水平＋垂直取平均；等價於 MAF |
| 多元線性迴歸殘差 | 用其他波段預測某波段，殘差視為雜訊 | HySime 使用；詳見第 4 節 |
| 均質區域 / 暗電流 | 在均勻區域（水體、校正板）直接算變異量 | 需要已知均勻區域 |

### 保留幾個成分
1. 看特徵值曲線，取 $\lambda > 2$（SNR > 1）
2. 逐張檢視 MNF 成分影像，從有空間結構變成雪花點的地方切（最常用）
3. 看特徵值拐點

### Inverse MNF（去雜訊）
保留前 $k$ 個成分、其餘設零，用 $W^{-1}$ 轉回原始波段空間 → 波段數與物理意義保留。常作為解混（VCA、PPI、N-FINDR）與分類的前處理。

### PCA vs MNF

| | PCA | MNF |
|---|---|---|
| 排序依據 | 變異量 | SNR |
| 需要雜訊估計 | 否 | 是 |
| 波段尺度縮放 | 結果會改變 | 不變（scale-invariant） |
| 雜訊大的波段 | 可能排到前面 | 排到後面 |
| 運算成本 | 一次特徵分解 | 兩次特徵分解＋雜訊估計 |

### 限制
- 假設雜訊是加性、空間不相關、平穩
- 推掃式（pushbroom）感測器的條紋雜訊空間相關 → shift difference 低估，條紋會跑到前面成分
- 場景空間變化劇烈（高解析都市區）→ 訊號差異被誤當雜訊，$\Sigma_N$ 高估
- 做 MNF 前先處理壞波段（水氣吸收帶）與條紋

---

## 3. 樣本沒有空間鄰域時（單像素樣本 / 點光譜）

**重點：MNF 轉換本身不受影響，只有雜訊估計那一步需要換方法。**
兩次特徵分解都在「樣本數 × 波段數」矩陣上做，不需要空間資訊。

### 替代的雜訊估計方法

**(1) 重複量測（最可靠）**
同一樣本量 $r$ 次，每次減掉該樣本平均即為純雜訊：

$$\Sigma_N = \frac{1}{n(r-1)} \sum_{i=1}^{n} \sum_{k=1}^{r} (x_{ik} - \bar{x}_i)(x_{ik} - \bar{x}_i)^T$$

- pooled within-sample covariance，含波段間雜訊相關結構
- 建議每個樣本至少量 3 次
- 若「重複」是同類別的不同樣本 → $\Sigma_N$ 變成類內共變異 → MNF 退化為 Fisher LDA

**(2) 光譜方向差分**
把 shift difference 從空間軸搬到波長軸：

$$\Sigma_N \approx \tfrac{1}{2}\,\text{Cov}(x_b - x_{b+1})$$

或用 Savitzky-Golay / 小波平滑後取殘差。
- 缺點：尖銳吸收特徵（red edge、礦物窄吸收帶）會被誤判為雜訊；相鄰波段雜訊常有相關性 → 低估

**(3) 多元線性迴歸殘差（HySime）**
只需光譜、不需空間鄰域。需 $n \gg B$，否則過擬合。詳見第 4 節。

**(4) 儀器端雜訊模型**
重複量暗電流與白板參考 → 各波段雜訊變異量 $\sigma_b^2$（只有對角線）。
- MNF 退化為「每波段除以雜訊標準差再做 PCA」（雜訊加權 PCA）
- 抓不到波段間雜訊相關；光子雜訊隨訊號強度變，暗電流只量到底噪

**(5) 像素來自影像**
回到整張影像用 shift difference 估 $\Sigma_N$，再把轉換矩陣套用到抽出的像素。雜訊是感測器性質，用整張影像估更準。

### 樣本數 < 波段數的問題
點光譜資料集常是幾十到幾百個樣本對幾百個波段，$\Sigma_X$ 秩 $\leq n-1$：
- 直接對資料矩陣做 SVD，不要先算 $B \times B$ 共變異矩陣
- $\Sigma_N$ 用對角或 Ledoit-Wolf shrinkage 正則化，確保可逆
- 先做波段選擇或合併（binning）降低 $B$

### 實務建議
- 點光譜儀（ASD、實驗室 FTIR）SNR 通常遠高於航空 / 衛星影像，MNF 效益較小
- 化學計量學慣例：SNV、一階微分、MSC 前處理 → PCA
- 有目標變數（濃度、類別）→ PLS / PLS-DA 通常比 MNF + 分類器有效
- MNF 划算的情境：雜訊在不同波段差異大，且能可靠估出

---

## 4. 多元線性迴歸殘差法（Regression-based Noise Estimation）

### 原理
第 $b$ 波段：$x_b = s_b + n_b$

- **訊號可預測**：訊號維度遠低於波段數（物質十幾種 vs 波段兩三百個），任一波段的訊號幾乎是其他波段訊號的線性組合
  $$s_b \approx \sum_{j \neq b} \alpha_j\, s_j$$
- **雜訊不可預測**：若各波段雜訊互相獨立，$n_b$ 與其他波段的任何線性組合都不相關

用其他波段對 $x_b$ 做最小平方迴歸，能解釋的是訊號，殘差就是雜訊：

$$\hat{\beta}_b = \arg\min_\beta \|x_b - Z_{-b}\beta\|^2, \qquad \varepsilon_b = x_b - Z_{-b}\hat{\beta}_b \approx n_b$$

$Z_{-b}$：去掉第 $b$ 欄的資料矩陣（$n \times (B-1)$）。
對每個波段做一次，得雜訊矩陣 $E$（$n \times B$），$\Sigma_N = E^T E / n$。

來源：Roger & Arnold (1996) 針對 AVIRIS 提出；Bioucas-Dias & Nascimento (2008) 在 HySime 中給出高效實作。

### 高效實作：不用真的做 B 次迴歸
天真做法：$B$ 次迴歸，每次解 $(B-1)\times(B-1)$ 系統，成本 $O(B^4)$。
HySime 的關鍵：只算一次 $R = X^T X$ 的逆，所有波段的殘差直接從 $R^{-1}$ 讀出。

**推導**：把第 $b$ 波段排到最後，
$$R = \begin{bmatrix} Z^TZ & Z^Tx_b \\ x_b^TZ & x_b^Tx_b \end{bmatrix}$$
分塊求逆後：
- $R^{-1}$ 右下角元素 = $1/\text{RSS}_b$（殘差平方和倒數）
- 最後一欄上半部 = $-\hat\beta_b / \text{RSS}_b$

整理得：

$$\hat{\beta}_b = -\frac{[R^{-1}]_{-b,\,b}}{[R^{-1}]_{b,b}}, \qquad \varepsilon_b = \frac{X\,[R^{-1}]_{:,\,b}}{[R^{-1}]_{b,b}}$$

整個雜訊矩陣一行完成：

$$E = X\,R^{-1}\,D^{-1}, \qquad D = \text{diag}\big(\text{diag}(R^{-1})\big)$$

每波段雜訊變異量：$\sigma_b^2 = 1 / (n\,[R^{-1}]_{b,b})$（同「逆相關矩陣對角線 = VIF」）

```python
import numpy as np

def regression_noise(X, ridge=1e-6):
    # X: n x B，每列一條光譜；先去均值等於迴歸含截距項
    X = X - X.mean(axis=0)
    n, B = X.shape
    R = X.T @ X / n
    R += ridge * np.trace(R) / B * np.eye(B)   # 避免病態
    Rinv = np.linalg.inv(R)
    E = X @ Rinv / np.diag(Rinv)               # n x B 雜訊估計
    Sigma_N = E.T @ E / n
    return E, Sigma_N
```

- 成本：一次 $O(nB^2)$ 矩陣乘法 + 一次 $O(B^3)$ 求逆
- 影像太大時隨機抽幾萬個像素估 $R$ 即可（$\Sigma_N$ 只有 $B \times B$）
- **Poisson 變體**（訊號相依雜訊）：對 $\sqrt{X}$ 做迴歸，殘差再乘回 $\sqrt{X}$ 還原尺度

### 兩個相反方向的偏差

| 偏差 | 方向 | 原因 | 修正 |
|---|---|---|---|
| 有限樣本過擬合 | 低估 | $B-1$ 個自由參數，$n$ 不夠時殘差被硬壓小；$n = B-1$ 時殘差為零 | RSS 除以 $(n-B+1)$ 而非 $n$；加 ridge；經驗法則 $n \geq 5B \sim 10B$ |
| 預測變數含雜訊 | 高估 | 其他波段的雜訊漏進殘差 | 係數分散在多波段、各自很小，通常可忽略 |

### 使用情境

**適合**
- 沒有空間鄰域的資料（點光譜、單像素樣本）
- 空間紋理豐富的影像（都市、農田邊界）：shift difference 會誤判訊號差異，迴歸法不受影響
- HySime 訊號子空間維度估計：$\Sigma_S = \Sigma_X - \Sigma_N$ → 估端元數量 → VCA 等解混的必要輸入
- 每波段 SNR 剖面：$\sigma_b^2$ 對波長作圖，水氣吸收帶、感測器邊緣波段明顯突出 → 壞波段偵測

**不適合**
- 多光譜（$B$ 只有幾個到十幾個）：波段冗餘不足，殘差含大量訊號
- 雜訊跨波段相關（增益抖動、共模雜訊、跨波段條紋）：會被當訊號解釋掉 → 嚴重低估
- 樣本數少於或接近波段數
- 某波段有其他波段沒有的獨特特徵（極窄發射線、單波段壞像素）

**常見變體**
- 只用相鄰 $k$ 個波段迴歸（局部視窗）：減少參數，適合 $n$ 不大
- SSDC（spectral-spatial decorrelation, Gao et al.）：同時用相鄰波段與相鄰像素當預測變數，影像資料上較穩健

### 驗證方法
把估出的 $E$ 拿去看自相關：真正的雜訊在波長方向與空間方向都應接近白雜訊。若殘差影像仍看得出地物輪廓 → $n$ 不夠或雜訊有結構，估計不可信。

---

## 5. 流程總結

```
高光譜資料 (n × B)
    │
    ├─ 有空間鄰域 ──→ shift difference 估 Σ_N
    ├─ 無空間鄰域 ──→ 重複量測 / 迴歸殘差 / 儀器雜訊模型
    │
    ▼
雜訊白化 (F = E Λ^{-1/2})
    │
    ▼
PCA on 白化資料 ──→ 特徵值 λ = 1 + SNR
    │
    ├─ 保留 λ > 2 或目視判斷
    │
    ▼
Inverse MNF 去雜訊 ──→ 解混 / 分類 / 端元估計
```

---

## 6. 參考文獻
- Green, A. A., Berman, M., Switzer, P., & Craig, M. D. (1988). A transformation for ordering multispectral data in terms of image quality with implications for noise removal. *IEEE TGRS*, 26(1), 65–74.
- Lee, J. B., Woodyatt, A. S., & Berman, M. (1990). Enhancement of high spectral resolution remote-sensing data by a noise-adjusted principal components transform. *IEEE TGRS*, 28(3), 295–304.
- Roger, R. E., & Arnold, J. F. (1996). Reliably estimating the noise in AVIRIS hyperspectral images. *IJRS*, 17(10), 1951–1962.
- Bioucas-Dias, J. M., & Nascimento, J. M. P. (2008). Hyperspectral subspace identification. *IEEE TGRS*, 46(8), 2435–2445.
- Gao, L., Zhang, B., Zhang, X., Zhang, W., & Tong, Q. (2008). A new operational method for estimating noise in hyperspectral images. *IEEE GRSL*, 5(1), 83–87.
