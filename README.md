# 🧬 GPCR-Nanodisc Integration Predictor

GPCR-Nanodisc Integration Predictor は、Gタンパク質共役受容体（GPCR）をナノディスク（Nanodisc）へ再構成する際の安定性を、アミノ酸配列のみから予測・評価するWebアプリケーションです。

ESMFoldによる構造予測、回転半径（$R_g$）に基づくMSPサイズの推奨、および疎水性ベルト解析による脂質との相性（疎水性ミスマッチ）のスコアリングを統合しています。

---

## 🌟 Key Features

* **Structure Prediction:** ESMFold APIを用いて配列から3D構造（PDB）を高速に生成。
* **Scaffold Recommendation:** タンパク質の広がり（$R_g$）から、最適なMSP（MSP1D1, MSP1E3D1, etc.）を自動提案。
* **Hydrophobic Belt Analysis:** Kyte-Doolittleスケールを用い、Z軸方向の疎水性プロファイルを算出。
* **Lipid Compatibility:** 選択した脂質（DMPC, POPC等）の厚さとタンパク質の膜貫通領域の整合性をスコアリング。
* **Interactive 3D View:** `py3Dmol` を活用し、予測された膜境界（Box）と構造をブラウザ上で視覚化。

---

## 🧪 Scientific Background

GPCRの構造機能解析において、膜環境を模倣したナノディスクへの再構成は極めて重要ですが、適切な**Scaffold Protein (MSP)**や**Lipid**の選択はしばしば試行錯誤を要します。

本ツールは以下の物理化学的指標に基づき、実験の成功率を高めるための意思決定を支援します。

1. **Radius of Gyration ($R_g$):** タンパク質の物理的なサイズを定義し、ナノディスク内に収まるかを確認します。
2. **Hydrophobic Mismatch:** 脂質二重層の疎水性コアの厚さとタンパク質の疎水性ベルトの幅が一致しない場合、タンパク質の変性や凝集のリスクが高まります。

---

## 🚀 Getting Started

### Installation
```bash
git clone [https://github.com/your-username/gpcr-nanodisc-predictor.git](https://github.com/your-username/gpcr-nanodisc-predictor.git)
cd gpcr-nanodisc-predictor
pip install -r requirements.txt
Usage
Bash
streamlit run app.py
🧑‍🔬 Author
[Your Name / ユーザー名]

Ph.D. in Agriculture

Specializing in Protein Engineering & Antibody Drug Discovery