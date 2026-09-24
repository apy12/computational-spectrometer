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
| 去雜訊（低秩近似） | $X = U \Sigma V^T$ ，只保留前 k 個奇異值重建，小奇異值的成分多為雜訊 |
| 避免 Hughes 現象 | 特徵維度高、樣本少會使分類精度下降；先降維再訓練分類器 |
| 估計本質維度 / 端元數 | 奇異值衰減曲線（scree plot）可估訊號子空間維度，是 VCA、N-FINDR 的前置步驟 |
| 視覺化 | 前三個 PC 當 RGB 顯示 |
| 異常 / 目標偵測 | 背景建模為低秩子空間，投影殘差大的像素為異常（RX detector 的基礎） |

### SVD 與 PCA 的關係

- PCA 等於對「去均值後的資料矩陣」做 SVD
- 右奇異向量 $V$ 就是共變異矩陣的特徵向量（主成分方向）
- 各主成分的變異量為 $\sigma_i^2 / (n-1)$
- 實務上直接做 SVD 比先算共變異矩陣再特徵分解數值更穩定

### 限制

- **只看變異量，不看訊噪比**：雜訊大的波段可能排到前面，此時改用 MNF
- 主成分是線性組合，失去物理光譜意義
- 小面積目標的訊息可能藏在低變異量成分，降維太狠會丟失
- 全域線性轉換

---

## 2. MNF（Minimum Noise Fraction）

### 動機

PCA 按變異量排序，MNF 按 **SNR** 排序：前面的成分最乾淨，後面的成分是純雜訊。（Green et al., 1988）

### 數學定義

資料等於訊號加雜訊：

$$
X = S + N , \qquad \Sigma_X = \Sigma_S + \Sigma_N
$$

找投影方向 $w$ 使雜訊比例（noise fraction）最小：

$$
\mathrm{NF}(w) = \frac{w^T \Sigma_N w}{w^T \Sigma_X w}
$$

等價於廣義特徵值問題：

$$
\Sigma_X w = \lambda \Sigma_N w , \qquad \lambda_i = 1 + \mathrm{SNR}_i
$$

- 特徵值由大到小排序
- 特徵值接近 1 的成分幾乎全是雜訊

### 兩步驟實作（Noise-Adjusted PCA, Lee et al., 1990）

**步驟一：雜訊白化**

1. 對雜訊共變異矩陣做特徵分解 $\Sigma_N = E \Lambda E^T$
2. 建立白化矩陣 $F = E \Lambda^{-1/2}$
3. 轉換資料 $Y = F^T X$ ，此時雜訊的共變異矩陣變成單位矩陣

**步驟二：對白化後的資料做 PCA**

- 對 $\Sigma_Y = F^T \Sigma_X F$ 做特徵分解，得特徵向量矩陣 $G$
- 因雜訊已壓平為 1，變異量排序就等於 SNR 排序

最終轉換矩陣：

$$
W = F G
$$

### 雜訊共變異矩陣的估計（影像資料）

| 方法 | 做法 | 備註 |
|---|---|---|
| 相鄰像素差分（shift difference） | $\Sigma_N \approx \frac{1}{2} \mathrm{Cov}(X(i,j) - X(i+1,j))$ | ENVI 預設；水平與垂直取平均；等價於 MAF |
| 多元線性迴歸殘差 | 用其他波段預測某波段，殘差視為雜訊 | HySime 使用；詳見第 4 節 |
| 均質區域 / 暗電流 | 在均勻區域（水體、校正板）直接算變異量 | 需要已知均勻區域 |

### 保留幾個成分

1. 看特徵值曲線，取 $\lambda > 2$ 的成分（即 SNR 大於 1）
2. 逐張檢視 MNF 成分影像，從有空間結構變成雪花點的地方切（最常用）
3. 看特徵值拐點

### Inverse MNF（去雜訊）

保留前 k 個成分、其餘設零，用 $W^{-1}$ 轉回原始波段空間，波段數與物理意義都保留。常作為解混（VCA、PPI、N-FINDR）與分類的前處理。

### PCA vs MNF

| | PCA | MNF |
|---|---|---|
| 排序依據 | 變異量 | SNR |
| 需要雜訊估計 | 否 | 是 |
| 波段尺度縮放 | 結果會改變 | 不變（scale-invariant） |
| 雜訊大的波段 | 可能排到前面 | 排到後面 |
| 運算成本 | 一次特徵分解 | 兩次特徵分解加雜訊估計 |

### 限制

- 假設雜訊是加性、空間不相關、平穩
- 推掃式（pushbroom）感測器的條紋雜訊空間相關，shift difference 會低估，條紋會跑到前面成分
- 場景空間變化劇烈（高解析都市區）時，訊號差異被誤當雜訊， $\Sigma_N$ 高估
- 做 MNF 前先處理壞波段（水氣吸收帶）與條紋

---

## 3. 樣本沒有空間鄰域時（單像素樣本 / 點光譜）

**重點：MNF 轉換本身不受影響，只有雜訊估計那一步需要換方法。**
兩次特徵分解都在「樣本數 × 波段數」矩陣上做，不需要空間資訊。

### 替代的雜訊估計方法

**(1) 重複量測（最可靠）**

同一樣本量 r 次，每次減掉該樣本平均即為純雜訊：

$$
\Sigma_N = \frac{1}{n(r-1)} \sum_{i=1}^{n} \sum_{k=1}^{r} (x_{ik} - \bar{x}_i)(x_{ik} - \bar{x}_i)^T
$$

- 這是 pooled within-sample covariance，含波段間雜訊相關結構
- 建議每個樣本至少量 3 次
- 若「重複」是同類別的不同樣本，則 $\Sigma_N$ 變成類內共變異，MNF 退化為 Fisher LDA

**(2) 光譜方向差分**

把 shift difference 從空間軸搬到波長軸：

$$
\Sigma_N \approx \frac{1}{2} \mathrm{Cov}(x_b - x_{b+1})
$$

或用 Savitzky-Golay / 小波平滑後取殘差。
缺點：尖銳吸收特徵（red edge、礦物窄吸收帶）會被誤判為雜訊；相鄰波段雜訊常有相關性，導致低估。

**(3) 多元線性迴歸殘差（HySime）**

只需光譜、不需空間鄰域。需要樣本數遠大於波段數，否則過擬合。詳見第 4 節。

**(4) 儀器端雜訊模型**

重複量暗電流與白板參考，得各波段雜訊變異量 $\sigma_b^2$ （只有對角線）。

- MNF 退化為「每波段除以雜訊標準差再做 PCA」（雜訊加權 PCA）
- 抓不到波段間雜訊相關；光子雜訊隨訊號強度變，暗電流只量到底噪

**(5) 像素來自影像**

回到整張影像用 shift difference 估 $\Sigma_N$ ，再把轉換矩陣套用到抽出的像素。雜訊是感測器性質，用整張影像估更準。

### 樣本數小於波段數的問題

點光譜資料集常是幾十到幾百個樣本對幾百個波段， $\Sigma_X$ 的秩最多是 n − 1：

- 直接對資料矩陣做 SVD，不要先算 B × B 的共變異矩陣
- $\Sigma_N$ 用對角或 Ledoit-Wolf shrinkage 正則化，確保可逆
- 先做波段選擇或合併（binning）降低 B

### 實務建議

- 點光譜儀（ASD、實驗室 FTIR）的 SNR 通常遠高於航空 / 衛星影像，MNF 效益較小
- 化學計量學慣例：SNV、一階微分、MSC 前處理，再做 PCA
- 有目標變數（濃度、類別）時，PLS / PLS-DA 通常比 MNF 加分類器有效
- MNF 划算的情境：雜訊在不同波段差異大，且能可靠估出

---

## 4. 多元線性迴歸殘差法（Regression-based Noise Estimation）

### 原理

第 b 波段的觀測值是訊號加雜訊：

$$
x_b = s_b + n_b
$$

**訊號可預測**：訊號維度遠低於波段數（物質十幾種 vs 波段兩三百個），任一波段的訊號幾乎是其他波段訊號的線性組合：

$$
s_b \approx \sum_{j \neq b} \alpha_j s_j
$$

**雜訊不可預測**：若各波段雜訊互相獨立，則 $n_b$ 與其他波段的任何線性組合都不相關。

用其他波段對 $x_b$ 做最小平方迴歸，能解釋的是訊號，殘差就是雜訊：

$$
\hat{\beta}_b = \arg\min_{\beta} \lVert x_b - Z_{-b} \beta \rVert^2 ,
\qquad
\varepsilon_b = x_b - Z_{-b} \hat{\beta}_b \approx n_b
$$

其中 $Z_{-b}$ 是去掉第 b 欄的資料矩陣，大小為 n × (B − 1)。
對每個波段做一次，得到 n × B 的雜訊矩陣 $E$ ，雜訊共變異矩陣為：

$$
\Sigma_N = \frac{E^T E}{n}
$$

來源：Roger & Arnold (1996) 針對 AVIRIS 提出；Bioucas-Dias & Nascimento (2008) 在 HySime 中給出高效實作。

### 高效實作：不用真的做 B 次迴歸

天真做法：做 B 次迴歸，每次解一個 (B − 1) × (B − 1) 的線性系統，成本 $O(B^4)$ 。
HySime 的關鍵：只算一次 $R = X^T X$ 的逆矩陣，所有波段的殘差直接從 $R^{-1}$ 讀出。

**推導**：把第 b 波段排到最後，

$$
R = \begin{bmatrix} Z^T Z & Z^T x_b \cr x_b^T Z & x_b^T x_b \end{bmatrix}
$$

分塊求逆後：

- $R^{-1}$ 右下角的元素等於 $1 / \mathrm{RSS}_b$ （殘差平方和的倒數）
- 最後一欄的上半部等於 $-\hat{\beta}_b / \mathrm{RSS}_b$

整理得：

$$
\hat{\beta}_b = - \frac{[R^{-1}]_{-b,b}}{[R^{-1}]_{b,b}} ,
\qquad
\varepsilon_b = \frac{X [R^{-1}]_{:,b}}{[R^{-1}]_{b,b}}
$$

整個雜訊矩陣一行完成：

$$
E = X R^{-1} D^{-1} , \qquad D = \mathrm{diag}\left( \mathrm{diag}(R^{-1}) \right)
$$

每個波段的雜訊變異量：

$$
\sigma_b^2 = \frac{1}{n [R^{-1}]_{b,b}}
$$

這跟統計學裡「逆相關矩陣的對角線是 VIF」是同一件事。

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

- 成本：一次 $O(nB^2)$ 的矩陣乘法加一次 $O(B^3)$ 的求逆
- 影像太大時隨機抽幾萬個像素估 R 即可（雜訊共變異矩陣只有 B × B）
- **Poisson 變體**（訊號相依雜訊）：對 $\sqrt{X}$ 做迴歸，殘差再乘回 $\sqrt{X}$ 還原尺度

### 兩個相反方向的偏差

| 偏差 | 方向 | 原因 | 修正 |
|---|---|---|---|
| 有限樣本過擬合 | 低估 | 有 B − 1 個自由參數，n 不夠時殘差被硬壓小；n = B − 1 時殘差為零 | RSS 除以 (n − B + 1) 而非 n；加 ridge；經驗法則 n ≥ 5B 到 10B |
| 預測變數含雜訊 | 高估 | 其他波段的雜訊漏進殘差 | 係數分散在多波段、各自很小，通常可忽略 |

### 使用情境

**適合**

- 沒有空間鄰域的資料（點光譜、單像素樣本）
- 空間紋理豐富的影像（都市、農田邊界）：shift difference 會誤判訊號差異，迴歸法不受影響
- HySime 訊號子空間維度估計：由 $\Sigma_S = \Sigma_X - \Sigma_N$ 估端元數量，是 VCA 等解混的必要輸入
- 每波段 SNR 剖面：把 $\sigma_b^2$ 對波長作圖，水氣吸收帶、感測器邊緣波段明顯突出，可當壞波段偵測

**不適合**

- 多光譜（B 只有幾個到十幾個）：波段冗餘不足，殘差含大量訊號
- 雜訊跨波段相關（增益抖動、共模雜訊、跨波段條紋）：會被當訊號解釋掉，嚴重低估
- 樣本數少於或接近波段數
- 某波段有其他波段沒有的獨特特徵（極窄發射線、單波段壞像素）

**常見變體**

- 只用相鄰 k 個波段迴歸（局部視窗）：減少參數，適合 n 不大
- SSDC（spectral-spatial decorrelation, Gao et al.）：同時用相鄰波段與相鄰像素當預測變數，影像資料上較穩健

### 驗證方法

把估出的 $E$ 拿去看自相關：真正的雜訊在波長方向與空間方向都應接近白雜訊。若殘差影像仍看得出地物輪廓，代表 n 不夠或雜訊有結構，估計不可信。

---

## 5. 實際設定：diffuser 均勻化光場、非同調光源

設定：樣本經 diffuser 後成像，sensor 看到近乎均勻的光場，取 ROI 平均當該樣本的光譜。光源為鹵素燈 / LED（非同調，無 speckle）。

### 拍攝協定

- 每個樣本連拍 5–10 張影格，不要只存平均
- 暗電流：遮光後拍 20–50 張，同曝光時間
- 白板參考：透過同一片 diffuser 拍 20–50 張取平均；每 10–20 個樣本或 15–30 分鐘重拍，追蹤燈源漂移
- 全程固定曝光與增益，雜訊才近似平穩
- 反射率：

$$
R = \frac{S - D}{W - D}
$$

### 兩個層級的雜訊

**像素層級**：ROI 內像素間的差異。訊號完全相同，shift difference 就是純雜訊，n 有上萬個像素，估得非常穩。這是讀出雜訊加光子雜訊，各像素獨立。

**影格層級**：同一樣本不同影格的 ROI 平均之間的差異。這是最後分析用的光譜真正的雜訊，等於像素雜訊除以 N，再加上所有像素共享的時間性雜訊（燈源閃爍、色溫漂移）：

$$
\Sigma_{N,\mathrm{mean}} = \frac{\Sigma_{N,\mathrm{pixel}}}{N} + \Sigma_{\mathrm{common}}
$$

實務上直接用影格間的 pooled within-sample covariance 估影格層級雜訊，自由度是「樣本數 × (影格數 − 1)」，通常仍小於 B，需加 shrinkage。共模雜訊通常是低秩的（整體增益抖動加平滑的光譜形狀變化），shrinkage 對它很友善。

**診斷**：比較兩者的對角線。若影格層級雜訊遠大於像素層級除以 N，代表共模雜訊主導，多平均像素沒用，要穩定燈源或多拍影格。

### 白板參考的貢獻

W 是多張平均，雜訊本身小，但它的誤差以相同形狀進到每個樣本，在 PCA 裡會變成一個共有成分。鹵素燈 / LED 通常不嚴重；自然光會隨雲層大幅漂移，需頻繁重拍參考，或接受第一個 PC 可能是照明變化並在分析時剔除。

### 波段相依的雜訊：MNF 的用武之地

反射率的雜訊變異量約為：

$$
\mathrm{Var}(R_b) \propto \frac{R_b}{W_b - D_b}
$$

鹵素燈在藍光端輸出弱、LED 有光譜空隙、感測器邊緣波段量子效率低，這些波段的 W − D 很小，反射率雜訊會大好幾倍。PCA 會把這些高雜訊波段排到前面，MNF 會正確地把它們壓到後面。做 MNF 前先看每波段 SNR 剖面，把明顯壞的波段砍掉。

樣本間的亮度依賴只要反射率沒差到一個數量級可以忽略；差很多時對原始計數做平方根或 Anscombe 轉換再算。

### MNF 實作

```python
import numpy as np
from sklearn.covariance import LedoitWolf

# frames[i]: (r, B) 第 i 個樣本每張影格的 ROI 平均反射率
r = frames[0].shape[0]
X = np.stack([f.mean(0) for f in frames])                            # n x B 樣本光譜
E = np.concatenate([f - f.mean(0, keepdims=True) for f in frames])   # 影格離差
Sigma_N = LedoitWolf(assume_centered=True).fit(E).covariance_ * r/(r-1)

w, V = np.linalg.eigh(Sigma_N)
F = V / np.sqrt(w)                            # 白化矩陣
Xc = X - X.mean(0)
U, s, Gt = np.linalg.svd(Xc @ F, full_matrices=False)   # 直接 SVD，不建 B×B 矩陣
lam = s**2 / (len(X) - 1)                     # ≈ 1 + SNR
W = F @ Gt.T                                  # MNF 轉換矩陣
scores = Xc @ W
```

n 小於 B 時直接對白化資料做 SVD，不要先組 $\Sigma_Y$ 。

### 自我一致性檢查

- 尾端的 $\lambda$ 應該接近 1。整體明顯大於 1 表示 $\Sigma_N$ 低估；小於 1 表示高估。這是 MNF 特有的免費檢查
- ROI 內像素的標準差與 shift difference 估出的雜訊應接近；差很多代表光場不夠均勻（diffuser 不夠、灰塵、暗角），縮小 ROI 或換 diffuser
- 影格層級的雜訊在波長軸上應接近白雜訊；若看得出平滑的光譜形狀，是燈源漂移，多拍影格或加穩壓

### 此設定下不需要的東西

- 迴歸殘差法：已有更直接的雜訊來源，可留作交叉驗證
- 空間去條紋：diffuser 已把空間結構抹掉

---

## 6. 應用：感測器響應曲線（sensor response profiles）的選擇

問題：設計多通道光譜感測器時，要根據應用情境（環境光、Raman）從候選的響應曲線池裡選出一組通道。

**核心觀念：分解的對象是目標訊號空間，不是感測器響應本身。**
感測器響應的分解只能看出通道有多冗餘，回答不了「這組通道適不適合這個應用」。

### 6.1 框架

$$
y = A s + n
$$

- $A$ ： m × B 感測矩陣，每一列是一個通道的響應曲線
- $s$ ： 待測光譜（B × 1）
- $y$ ： m 個通道的讀值
- 目標：從 $y$ 推回 $s$ 或推回真正要的量（CCT、濃度、峰強度）的誤差最小

### 6.2 目標光譜的基底

先建目標光譜庫（環境光：各種 CCT 的日光、LED、螢光燈、白熾燈；Raman：分析物光譜加螢光背景），做 PCA：

$$
s \approx P c , \qquad P \in \mathbb{R}^{B \times k}
$$

- k 是這類光譜的自由度，也是通道數的下限
- 日光只需要 3 個基底（CIE 的 S0、S1、S2 就是日光 PCA 的結果，見第 7 節）
- 加入 LED 與螢光燈後約需 6–10 個
- Raman 取決於分析物數量加背景的平滑成分

### 6.3 評估感測器組合：把 A 投影到基底上

$$
M = A P \quad (m \times k)
$$

$M$ 的奇異值代表每個基底方向被這組感測器「看到」的強度。某個奇異值接近零，表示該方向的光譜變化在輸出裡幾乎沒反應，後處理救不回來。

考慮通道雜訊後：

$$
M_w = \Sigma_n^{-1/2} A P
$$

$M_w$ 的奇異值直接是各基底方向的 SNR。這與 MNF 是同一個邏輯（雜訊白化後看訊號子空間），只是白化的是感測器通道的雜訊，而且在設計階段就能算出來。

### 6.4 設計指標

貝氏框架下，基底係數的先驗 $\Sigma_c$ 取 PCA 特徵值的對角矩陣，後驗共變異矩陣：

$$
\Sigma_{\mathrm{post}} = \left( \Sigma_c^{-1} + M^T \Sigma_n^{-1} M \right)^{-1}
$$

| 指標 | 定義 | 意義 |
|---|---|---|
| A-optimal | 最小化 $\mathrm{tr}(\Sigma_{\mathrm{post}})$ | 期望重建誤差 |
| D-optimal | 最大化 $\log\det(\Sigma_c^{-1} + M^T \Sigma_n^{-1} M)$ | 讀值與光譜的互資訊 |
| 任務導向 | 解 $M^T w \approx P^T f$ | 線性泛函 $f^T s$ 能否被量到 |

任務導向指標：若只需要一個線性泛函（三刺激值、Raman 峰積分、濃度迴歸向量），要求 $f^T P$ 落在 $M$ 的列空間裡。最小平方殘差是「原理上量不到的部分」， $w^T \Sigma_n w$ 是雜訊放大量。這是 Luther 條件的推廣。

```python
def design_score(A_sub, P, Sigma_c, Sigma_n):
    # A_sub: m x B 候選通道響應；P: B x k 目標光譜基底
    M = A_sub @ P
    info = np.linalg.inv(Sigma_c) + M.T @ np.linalg.solve(Sigma_n, M)
    post = np.linalg.inv(info)
    return np.trace(post), np.linalg.slogdet(info)[1]   # A-opt 越小越好，D-opt 越大越好
```

選擇流程：從候選池做 greedy forward selection，每次加入讓 $\mathrm{tr}(\Sigma_{\mathrm{post}})$ 下降最多的通道，直到邊際改善小於門檻。log det 是次模函數，greedy 有接近最佳的理論保證。

### 6.5 感測器響應本身的分解能看什麼

對候選池的 $A$ 做 SVD，奇異值衰減曲線顯示這些響應曲線實際張成幾維。兩條幾乎共線的濾光片會貢獻一個很小的奇異值，代表冗餘。但冗餘不一定沒用（重複通道可平均雜訊），最終仍以 6.4 的雜訊感知指標決定。

### 6.6 兩種應用的差異

**環境光**
- 目標光譜平滑、低維，寬帶重疊濾光片就有效（類似人眼三錐體）
- Luther 條件：要估 XYZ，CIE 色匹配函數須近似落在感測器響應的線性張成裡
- 需要紅外與紫外截止，日光和白熾燈在可見光外能量很大，會污染寬帶通道

**Raman**
- 訊號是稀疏窄峰疊在平滑螢光背景上，光子極少
- PCA 基底不適合（峰位固定但強度獨立變化），改用 NMF 或「已知峰位 + 低階多項式背景」字典
- 核心矛盾：窄帶濾光片特異性好但丟光子；寬帶重疊（Hadamard 型）保留通量但需後端解混
- SNR 是綁死的約束，通道數增加代表每通道光子減少，指標要在固定總曝光下比較
- MOE（multivariate optical element）：把 PLS 迴歸向量做成濾光片穿透率，一個通道直接讀出濃度

### 6.7 容易踩的坑

- PCA 基底按變異量排序，不是按任務重要性；有目標變數時用 PLS 或 LDA 方向，或至少用任務導向指標驗證
- 光譜庫不夠多樣，k 會低估；留一部分光譜庫做驗證
- 實體濾光片響應非負，做不出有負瓣的曲線，只能靠通道相減實現，會放大雜訊
- 製程公差：對 $A$ 加擾動做蒙地卡羅，確認指標對公差不敏感

---

## 7. 驗證：日光只需要 3 個基底

腳本：`daylight_basis_check.py`（需 `pip install colour-science`）

### 7.1 為什麼不能拿 CIE D 系列做 PCA

CIE 日光模型本身就是三個基底疊出來的：

$$
S_D(\lambda) = S_0(\lambda) + M_1 S_1(\lambda) + M_2 S_2(\lambda)
$$

對它做 PCA 必然得到 3，是循環論證。

### 7.2 方法

用與 CIE 無關的簡化物理大氣模型（類 SPCTRAL2）產生 3000 條日光光譜：

- 太陽：5778 K Planck
- 直射透過率：Rayleigh 散射、Ångström 氣溶膠（濁度 β、指數 α、單次散射反照率）、臭氧 Chappuis 帶、水氣 720 nm、氧氣 760 nm
- 天空光：Rayleigh 散射項（偏藍）加氣溶膠前向散射項（偏白）
- 雲：近乎光譜中性的漫射項，以雲量混合
- 直射可見比例 0–1（陽光 vs 陰影）
- 8 個物理參數獨立隨機取樣，光譜在 560 nm 正規化為 1

然後做 PCA，看累積變異、光譜 RMSE、色度誤差 Δu'v' 隨基底數的變化，並與 CIE S1/S2 比較子空間夾角。

### 7.3 結果

| k | 累積變異 | 光譜 RMSE | Δu'v' 中位 | Δu'v' 95% |
|---|---|---|---|---|
| 1 | 84.9% | 0.034 | 0.0014 | 0.0056 |
| 2 | 98.4% | 0.012 | 0.0006 | 0.0030 |
| **3** | **99.4%** | **0.007** | **0.0003** | **0.0013** |
| 4 | 99.8% | 0.004 | 0.0002 | 0.0009 |
| 5 | 99.95% | 0.002 | 0.0001 | 0.0005 |

Δu'v' 可辨閾值約 0.002。

- **3 個基底時 95% 的光譜色度誤差都在閾值以下**：「日光只需要 3 個基底」是色度學上的結論
- 光譜保真（RMSE < 0.2%）需要約 5 個
- 與 Hernández-Andrés 等人 2001 年用 2600 條實測 Granada 日光的結論一致：色度 3 個、光譜 5–7 個

### 7.4 與 CIE 基底的比較

- 第一變異方向與 CIE S1 夾角 9°（幾乎重合，即 CCT 軸）
- 第二方向與 S2 夾角 26°（部分重合）
- PC1 主要由直射可見比例驅動（陽光 vs 陰影），PC2、PC3 由太陽天頂角驅動，PC4 由水氣驅動

### 7.5 限制

- 8 個輸入參數只產生 3 個主要成分：低維性來自物理（大氣效應在光譜上平滑且高度共線），不是取樣巧合
- 簡化模型的高階基底形狀不會與實測完全相同；太陽光譜用 Planck 近似，沒有 Fraunhofer 結構
- 嚴格驗證：把實測資料庫（Granada、Judd 1964 的 622 條）餵給腳本的 `analyze(spectra, wl)`，流程相同

---

## 8. HySime 用在像素 responses：多通道感測器與重建式光譜儀

### 8.1 HySime 的輸入本來就是像素

HySime 把影像攤成「n 個像素 × B 個波段」矩陣，場景像素當樣本正是它的用法。它分兩段：

1. 用迴歸殘差估 $\Sigma_n$（第 4 節）
2. 在訊號共變異矩陣的特徵向量上逐一判斷「該方向訊號功率是否大於雜訊功率」，數出 k

第 2 段只需要 $\Sigma_y$ 和 $\Sigma_n$ ，**不在乎 $\Sigma_n$ 從哪來**。第 1 段不適用時，第 2 段仍可搭配外部提供的 $\Sigma_n$ 使用。

需要的三個條件：n ≫ B（5–10 倍）；通道冗餘 k ≪ B；雜訊在通道間不相關。

### 8.2 對照兩種像素

| 情況 | 迴歸估雜訊 | 補救 |
|---|---|---|
| 高光譜像素（B 約 100–300），n = 300 | n 與 B 同量級，過擬合、殘差趨近零 | 局部 ±10 波段迴歸、ridge、像素連續時用 shift difference；第 2 段不受影響 |
| 設計過的多通道感測器（m 約 4–32） | 通道不冗餘，殘差混入訊號， $\Sigma_n$ 高估、k 低估 | $\Sigma_n$ 從別處來（8.3），只用第 2 段 |

### 8.3 多通道情況下 Σ_n 的來源（依可靠度）

1. **靜態場景多拍幾張**：影格間變異就是雜訊，在使用的資料層級上估
2. **暗電流加平場的雜訊模型**： $\sigma^2 = \sigma_{\mathrm{read}}^2 + y / g$ ，g 為每 DN 的電子數；通道獨立，對角但隨訊號變
3. **樣本共變異的尾端特徵值**：各通道先除以雜訊標準差，尾端 m − k 個特徵值應平坦且接近 1；用 Marchenko–Pastur 分布檢查
4. **均質區域的像素間差分**

### 8.4 去噪：冗餘從像素維度或時間維度來

通道少時光譜方向沒有冗餘，但 300 個像素通常只含少數材質，訊號集中在 k 維子空間，像素間的冗餘就是去噪來源。

```python
def hysime_subspace(Sigma_y, Sigma_n):
    # Sigma_y: 觀測共變異 (m x m)；Sigma_n: 雜訊共變異，任何來源
    Rx = Sigma_y - Sigma_n
    Rx = (Rx + Rx.T) / 2
    _, E = np.linalg.eigh(Rx)
    E = E[:, ::-1]
    p = np.einsum('ji,jk,ki->i', E, Sigma_y, E)    # 各方向觀測功率
    s = np.einsum('ji,jk,ki->i', E, Sigma_n, E)    # 各方向雜訊功率
    keep = (-p + 2 * s) < 0                          # 訊號功率 > 雜訊功率
    return int(keep.sum()), E[:, keep]

def wiener_denoise(Y, Sigma_y, Sigma_n):
    # Y: n x m 像素讀值；子空間投影的軟版本
    mu = Y.mean(0)
    Sigma_s = Sigma_y - Sigma_n
    w, V = np.linalg.eigh(Sigma_s)
    Sigma_s = (V * np.clip(w, 0, None)) @ V.T       # 裁成半正定
    G = Sigma_s @ np.linalg.inv(Sigma_s + Sigma_n)
    return mu + (Y - mu) @ G.T
```

HySime 判斷式 $p_i > 2 \sigma_i^2$ 等價於「該方向 SNR 大於 1」，與 MNF 的 $\lambda > 2$ 是同一件事；差別在 HySime 用訊號共變異的特徵向量，MNF 用廣義特徵向量。

**目的是重建光譜時不要先去噪再重建**：貝氏重建本身含去噪，

$$
\hat{s} = \Sigma_s A^T \left( A \Sigma_s A^T + \Sigma_n \right)^{-1} y
$$

先去噪再重建等於收縮兩次，會壓掉弱訊號。雜訊只處理一次，放在重建裡。

### 8.5 重建式光譜儀：300 個像素各有任意形狀的響應

設定：diffuser 均勻化光場，300 個像素看到同一條光譜 $s$ ，每個像素有自己的響應曲線 $a_i(\lambda)$ ，可能多峰、寬帶、峰位不按順序：

$$
y_i = \int a_i(\lambda) s(\lambda) d\lambda + n_i , \qquad y = A s + n , \quad A \in \mathbb{R}^{300 \times B}
$$

**「順序亂掉」對子空間方法完全沒影響**

PCA、SVD、MNF、HySime 全迴歸都是通道置換不變的：重排 $y$ 的元素，共變異矩陣只是行列跟著重排，特徵值、子空間、殘差不變。

會壞掉的是假設相鄰通道相似的方法：光譜方向 shift difference、局部視窗迴歸、Savitzky-Golay、「逐張看成分影像判斷 k」。這些全部不能用，k 只能用特徵值或 HySime 判斷式決定。

**「多峰、寬帶」對迴歸估計反而是好事**

冗餘來自 $s$ 的低維性，不是響應曲線形狀。 $s \approx P c$ 只有 k 個自由度，所以 $y = A P c$ 不管 $A$ 長什麼樣都落在 k 維子空間，冗餘度 300 / k，任一像素的訊號都能被其他像素完美線性預測。300 個像素是不同物理像素，雜訊獨立，比光柵光譜儀相鄰波段有 crosstalk 更符合假設。前提仍是 n ≫ 300：要有幾千條不同入射光譜堆成 Y。

**單次量測的雜訊估計：用校正過的 A**

一次量測給 300 個數字，訊號只有 k 個自由度，剩下 300 − k 個自由度全是雜訊加模型誤差。只要 $A$ 有校正（單色儀掃描每個像素響應），單次量測就能自己估雜訊：

```python
def reconstruct_and_check(y, A, P, Sigma_c, sigma2):
    # y: (m,) 一次量測；A: m x B 校正響應；P: B x k 目標光譜基底
    # Sigma_c: k x k 基底係數先驗；sigma2: (m,) 各像素雜訊變異量（暗電流 + 光子雜訊模型）
    M = A @ P
    Wn = 1.0 / sigma2                                   # 對角雜訊白化
    info = np.linalg.inv(Sigma_c) + (M.T * Wn) @ M
    c_hat = np.linalg.solve(info, (M.T * Wn) @ y)       # MAP 係數
    s_hat = P @ c_hat
    r = y - M @ c_hat
    chi2 = np.sum(r**2 * Wn) / (len(y) - M.shape[1])    # 應接近 1
    return s_hat, chi2, np.linalg.inv(info)             # 後驗共變異可轉成光譜誤差條
```

`chi2` 是即時自我檢查：接近 1 表示雜訊模型與校正都對；明顯大於 1 表示 $A$ 漂移（濾光片峰位隨溫度移動）或 $s$ 跑出光譜庫子空間；小於 1 表示雜訊模型高估。每一發都能做，比 HySime 更直接。

**去噪就是重建**

MAP 解已是最佳線性去噪，300 個像素對 k 個未知數，平均增益約 $\sqrt{300 / k}$ 。先驗依應用換：環境光用光譜庫 PCA 基底加高斯先驗；Raman 改成峰位字典加稀疏先驗（NNLS / LASSO），因為 Raman 光譜是稀疏而非低秩，高斯先驗會抹平窄峰。

**多峰寬帶響應真正要擔心的兩件事**

- *雜訊異質性*：寬帶像素光子多、SNR 高；窄峰像素光子少、SNR 低，可差一個數量級。 $\Sigma_n$ 不能當成 $\sigma^2 I$ ，要逐像素建模並白化；設計指標看 $\Sigma_n^{-1/2} A P$ 的奇異值
- *校正誤差比光子雜訊更危險*：尖銳多峰響應對峰位漂移極敏感， $\Delta A \cdot s$ 在像素間是相關的（同一溫度漂移影響所有像素），不符合雜訊獨立假設，也不會被 300 對 k 的平均消掉。對策：chi2 監控、溫度補償校正表、或把 $A$ 的不確定性當成額外共變異項加進 $\Sigma_n$ 。寬帶平滑響應對漂移較穩健，這是多峰與寬帶的取捨

**什麼時候還是用 HySime**

累積幾千次量測後，把 Y（n × 300）丟給 HySime 當交叉驗證：估出的 $\Sigma_n$ 應與暗電流模型一致，估出的 k 應與光譜庫 PCA 維度一致。對不上通常是校正漂移或光譜庫不夠多樣。

---

## 9. 流程總結

```
高光譜資料 (n × B)
    │
    ├─ 有空間鄰域 ──→ shift difference 估 Σ_N
    ├─ 無空間鄰域 ──→ 重複量測 / 迴歸殘差 / 儀器雜訊模型
    ├─ diffuser 均勻化 ──→ ROI 內像素差分（像素層級）+ 多影格離差（影格層級）
    │
    ▼
雜訊白化 (F = E Λ^{-1/2})
    │
    ▼
PCA on 白化資料 ──→ 特徵值 λ = 1 + SNR
    │
    ├─ 保留 λ > 2 或目視判斷
    ├─ 檢查尾端 λ ≈ 1
    │
    ▼
Inverse MNF 去雜訊 ──→ 解混 / 分類 / 端元估計


感測器設計 (選通道)
    │
    ├─ 目標光譜庫 ──→ PCA 基底 P (k 維)
    ├─ 候選通道 A ──→ M_w = Σ_n^{-1/2} A P
    │
    ▼
A-opt / D-opt / 任務導向指標 ──→ greedy 選通道 ──→ 公差蒙地卡羅


重建式光譜儀 (300 像素, 任意響應, 單次量測 y)
    │
    ├─ 校正 A ──→ M = A P
    ├─ 逐像素雜訊模型 σ_i² ──→ 白化
    │
    ▼
MAP 重建 ŝ = P ĉ ──→ chi2 = rᵀΣ_n⁻¹r / (m − k) ≈ 1 ?
    │
    ├─ 累積 n ≫ 300 次量測後 ──→ HySime 交叉驗證 Σ_n 與 k
```

---

## 10. 參考文獻

- Green, A. A., Berman, M., Switzer, P., & Craig, M. D. (1988). A transformation for ordering multispectral data in terms of image quality with implications for noise removal. *IEEE TGRS*, 26(1), 65–74.
- Lee, J. B., Woodyatt, A. S., & Berman, M. (1990). Enhancement of high spectral resolution remote-sensing data by a noise-adjusted principal components transform. *IEEE TGRS*, 28(3), 295–304.
- Roger, R. E., & Arnold, J. F. (1996). Reliably estimating the noise in AVIRIS hyperspectral images. *IJRS*, 17(10), 1951–1962.
- Bioucas-Dias, J. M., & Nascimento, J. M. P. (2008). Hyperspectral subspace identification. *IEEE TGRS*, 46(8), 2435–2445.
- Gao, L., Zhang, B., Zhang, X., Zhang, W., & Tong, Q. (2008). A new operational method for estimating noise in hyperspectral images. *IEEE GRSL*, 5(1), 83–87.
- Judd, D. B., MacAdam, D. L., & Wyszecki, G. (1964). Spectral distribution of typical daylight as a function of correlated color temperature. *JOSA*, 54(8), 1031–1040.
- Hernández-Andrés, J., Romero, J., Nieves, J. L., & Lee, R. L. (2001). Color and spectral analysis of daylight in southern Europe. *JOSA A*, 18(6), 1325–1335.
- Nelson, M. P., Aust, J. F., Dobrowolski, J. A., Verly, P. G., & Myrick, M. L. (1998). Multivariate optical computation for predictive spectroscopy. *Analytical Chemistry*, 70(1), 73–82.
- Bird, R. E., & Riordan, C. (1986). Simple solar spectral model for direct and diffuse irradiance on horizontal and tilted planes at the Earth's surface for cloudless atmospheres. *J. Climate Appl. Meteor.*, 25(1), 87–97. (SPCTRAL2)
