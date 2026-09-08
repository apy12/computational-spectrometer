# 計算光譜儀與字典學習(Dictionary Learning)技術筆記

## 1. 計算光譜儀(Computational Spectrometer)是什麼

計算光譜儀的核心想法,是用比波長種類少很多的量測通道,搭配運算方法把完整光譜「猜」回來,而不是像傳統光柵光譜儀那樣,一個感測器對應一個波長。要量 N 個波長,傳統作法需要 N 個獨立通道;計算光譜儀改用 M 個(M ≪ N)經過特殊設計的濾波通道,每個通道對多個波長都有反應、反應曲線各不相同,再靠演算法從少量量測值反推出完整光譜。

Dictionary learning(字典學習)正是讓這個「反推」變得可靠的關鍵技術之一。

## 2. Dictionary Learning:ELI5

想像你要形容很多很多種顏色配方,每一種配方就是一條光譜。如果每種顏色都要準備一支專屬蠟筆,你需要無限多支蠟筆,太麻煩了。但你會發現:顏色雖然千變萬化,只要有一小盒「常用色蠟筆」,大部分顏色只要混合其中兩三支就能調出來,不必每種顏色都各自準備一支。

字典學習做的事,就是自動「看過」很多真實世界的光譜範例後,幫你挑出這盒最好用的「常用色蠟筆」——這盒蠟筆就叫做字典(D),每一支蠟筆叫做一個 atom。挑得好的話,幾乎任何一條新光譜,都能只用其中少數幾支的組合就精確描述出來。這種「只需要少數幾支就夠」的性質,叫做**稀疏表示法(sparse representation)**。

## 3. 數學模型:為什麼需要字典學習

### 3.1 量測模型與欠定問題

計算光譜儀的量測過程可以寫成:

```
y = Φx
```

- `x`:完整光譜(N 維向量,未知)
- `Φ`:感測矩陣,每一列是一個濾波通道的光譜響應
- `y`:量測結果(只有 M 個數字,M ≪ N)

因為 M < N,這是一個「欠定」方程組:理論上有無限多組 x 都能滿足這個等式,單靠 y 解不出唯一答案。

### 3.2 稀疏假設

字典學習提供一個額外假設來「收斂」答案:假設真實光譜 x 可以寫成

```
x ≈ Dα
```

其中 D 是 N×K 的字典矩陣(K 個 atom),α 是係數向量,且絕大多數為 0(稀疏)。代入量測方程:

```
y = Φx ≈ ΦDα = Aα
```

只要組合矩陣 A = ΦD 滿足某些數學性質(如 restricted isometry property),即使 M 遠小於 K,也能透過凸優化(如 LASSO / Basis Pursuit)或貪婪演算法(如 OMP)穩定地把稀疏的 α 找出來——這正是壓縮感知(compressive sensing)理論的核心結果。

### 3.3 字典 D 怎麼學出來:K-SVD

字典學習用一批真實範例光譜當訓練資料,自動找出最適合的 D。以 K-SVD 為例,是交替最佳化的兩步驟循環:

1. **稀疏編碼步驟**:固定字典 D,對每條訓練光譜 xᵢ,用 OMP 等演算法求出最稀疏的係數 αᵢ,使 xᵢ ≈ Dαᵢ 誤差最小。
2. **字典更新步驟**:固定所有 αᵢ,一次更新一個 atom(K-SVD 用 SVD 分解殘差矩陣來更新),讓整體重建誤差進一步下降。

兩步不斷交替直到收斂。最後得到的 D,每個 atom 代表訓練資料中反覆出現的「典型光譜形狀」——例如某材料特有的吸收峰、某光源的發光輪廓等。

### 3.4 對計算光譜儀的意義

字典是針對應用場景(特定材料庫、拉曼光譜、螢光光譜、遙測影像等)學出來的,比通用基底(傅立葉、小波)更「懂」該領域光譜的樣子,能用更少非零係數(更稀疏)描述真實光譜。稀疏度越高,理論上需要的量測通道數 M 就能越少——直接對應到光譜儀可以做得更小、更便宜,只用少量物理元件(濾波片、奈米結構、量子點濾波器等)就能達到不錯的重建品質。

## 4. 重建穩定性:Atom 數量與 Coherence

### 4.1 Mutual coherence 的定義

字典裡每個 atom 可視為高維空間中的一個向量。把所有向量正規化成單位長度後,mutual coherence 定義為任兩個不同 atom 內積絕對值的最大值:

```
μ(A) = max_{i≠j} |⟨a_i, a_j⟩|
```

- 兩向量夾角越小 → 內積越接近 1 → coherence 越高 → 越難分辨是哪個 atom 在起作用
- 夾角越接近 90° → 內積越接近 0 → coherence 越低 → 越容易分辨

### 4.2 稀疏恢復的理論界線

稀疏恢復理論(Donoho & Elad、Tropp 等)給出:若真正的稀疏係數 α 只有 s 個非零值,且

```
s < (1 + 1/μ(A)) / 2
```

則用 OMP 或 L1 最小化都能保證唯一、正確地找回 α。反過來解 μ:

```
μ(A) < 1 / (2s − 1)
```

例:若光譜通常只需 3 個 atom 組合(s=3),coherence 必須小於 1/5=0.2 才能理論上保證唯一解。s 越大,對 coherence 要求越嚴格。此界線是最壞情況下的理論保證,實務上有雜訊時通常需要更低的 coherence 才穩定。

### 4.3 Atom 數量 K 的取捨

| 情況 | 影響 |
|---|---|
| K 太少 | 字典表達能力不足,許多真實光譜無法用少數 atom 精確逼近,重建誤差(approximation error)上升 |
| K 太多(過完備) | 表達能力更強,但 atom 被塞進同一個 N 維空間,彼此更容易靠近 → coherence 上升;也容易 overfitting、稀疏編碼計算成本增加 |

實務上通常用交叉驗證:保留部分光譜不參與訓練,只用來測試不同 K 值下的重建誤差,找出誤差最小、coherence 又不會太高的 K。部分 K-SVD 變形會在字典更新步驟加入「鼓勵 atom 互不相關」的懲罰項,主動抑制 coherence 上升。

### 4.4 Φ 與 D 的耦合:不只是 D 自己的問題

真正決定重建穩定性的,是組合矩陣 A = ΦD 的 coherence,而非 D 單獨在原本 N 維空間裡的 coherence。

Φ 是降維映射(N 維壓到 M 維,M ≪ N),降維必然丟失資訊。其中一種後果是:原本在 N 維空間裡夾角很大、容易分辨的兩個 atom,經 Φ 投影到 M 維量測空間後,夾角可能被壓縮得更小——即使 D 自己設計得很好(atom 彼此低相關),Φ 選不好一樣可能讓 A=ΦD 的 coherence 大幅升高,讓理論上可穩定重建的稀疏度 s 變得不可靠。

### 4.5 實務策略

- **感測矩陣最佳化(optimized projections)**:D 學好後,反過來調整 Φ,目標讓 A=ΦD 的 Gram 矩陣(所有欄位兩兩內積組成的矩陣)盡量接近單位矩陣,主動壓低 coherence。
- **聯合設計(joint design of Φ and D)**:若 Φ 的物理實現有彈性(如可設計的薄膜濾波片、metasurface),可把 Φ 也當成可學習參數,跟 D 一起針對「投影後」的 coherence 做最佳化。
- **限制 K 並搭配上述兩者**:適度控制字典大小,搭配感測矩陣優化,通常比單獨加大 K 更能兼顧表達能力與重建穩定性。

## 5. 小結

- Atom 數量與 Φ-D 相容性其實是同一問題的兩面:都是在問「經過實際量測後,系統還能不能分辨是哪幾個 atom 在起作用」。
- 字典學習的價值,在於針對特定應用領域自動找出「最省、最能代表真實光譜」的基本形狀組合,讓計算光譜儀能用遠少於波長數的量測通道,穩定重建出完整光譜。
- 實務設計時需同時考慮:字典表達能力(K)、字典內部相關性、感測矩陣 Φ 的設計,以及三者交互作用下的重建穩定性。

## 6. 參考文獻

**字典學習基礎演算法**

- Aharon, M., Elad, M., & Bruckstein, A. (2006). K-SVD: An algorithm for designing overcomplete dictionaries for sparse representation. *IEEE Transactions on Signal Processing*, 54(11), 4311–4322. https://doi.org/10.1109/TSP.2006.881199
- Engan, K., Aase, S. O., & Husøy, J. H. (1999). Method of optimal directions for frame design. In *Proc. IEEE ICASSP*, Vol. 5, pp. 2443–2446.

**壓縮感知與稀疏恢復理論**

- Donoho, D. L. (2006). Compressed sensing. *IEEE Transactions on Information Theory*, 52(4), 1289–1306.
- Donoho, D. L., & Elad, M. (2003). Optimally sparse representation in general (nonorthogonal) dictionaries via ℓ1 minimization. *Proceedings of the National Academy of Sciences*, 100(5), 2197–2202.
- Tropp, J. A. (2004). Greed is good: Algorithmic results for sparse approximation. *IEEE Transactions on Information Theory*, 50(10), 2231–2242.

**感測矩陣與字典的聯合最佳化(Coherence 最佳化)**

- Elad, M. (2007). Optimized projections for compressed sensing. *IEEE Transactions on Signal Processing*, 55(12), 5695–5702. https://doi.org/10.1109/TSP.2007.900760
- Duarte-Carvajalino, J. M., & Sapiro, G. (2009). Learning to sense sparse signals: Simultaneous sensing matrix and sparsifying dictionary optimization. *IEEE Transactions on Image Processing*, 18(7), 1395–1408.

**計算光譜儀應用與回顧文獻**

- Kim, C., Park, D., & Lee, H.-N. (2020). Compressive sensing spectroscopy using a residual convolutional neural network. *Sensors*, 20(3), 594. https://doi.org/10.3390/s20030594
- Zhang, G., Xu, T., Sun, B., et al. (2025). Parallel dictionary reconstruction and fusion for spectral recovery in computational imaging spectrometers. *Sensors*, 25(15), 4556. https://doi.org/10.3390/s25154556
- Xue, Q., Yang, Y., Ma, W., Zhang, H., Zhang, D., Lan, X., Gao, L., Zhang, J., & Tang, J. (2024). Advances in miniaturized computational spectrometers. *Advanced Science*, 11, 2404448. https://doi.org/10.1002/advs.202404448
- Guan, Q., Lim, Z. H., Sun, H., Chew, J. X. Y., & Zhou, G. (2023). Review of miniaturized computational spectrometers. *Sensors*, 23(21), 8768. https://doi.org/10.3390/s23218768
