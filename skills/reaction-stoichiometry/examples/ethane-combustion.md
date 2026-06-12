# 示例：乙烷燃烧 —— 配平 / 化学计量 / 产率 / 绿色指标 / 图表

本示例演示 `reaction-stoichiometry` Skill 的 `stoich` 子命令完整能力：从未配平方程式出发，
配平、判断限制试剂、计算理论产量与剩余量、根据实际产量算百分产率、给出原子经济性，
并**生成一张产量柱状图（PNG 文件）**——这是普通 LLM 做不到的。

## 用户提示（飞书）

```text
反应 C2H6 + O2 -> CO2 + H2O。我有 10 g 乙烷和 40 g 氧气，实际只得到了 25 g CO2。
请配平方程式，判断谁是限制试剂，算出 CO2、H2O 的理论产量、氧气剩余量、CO2 的百分产率
和原子经济性，并画一张产量对比图。
```

## Agent 调用脚本

```bash
python "<SKILL_DIR>/scripts/stoichiometry_helper.py" stoich \
  --equation "C2H6 + O2 -> CO2 + H2O" \
  --given "C2H6=10g, O2=40g" \
  --actual "CO2=25g" \
  --target CO2 \
  --plot "compound_outputs/ethane/yield.png"
```

## 脚本输出（节选）

```text
[RESULT] balanced=2 C2H6 + 7 O2 -> 4 CO2 + 6 H2O
[RESULT] percent_yield[CO2]=85.41%
[RESULT] limiting_reagent=C2H6
[RESULT] extent_mol=0.166279
[RESULT] atom_economy=61.96%
[RESULT] chart_file=.../compound_outputs/ethane/yield.png
[JSON] {"mode": "stoich", "balanced": "2 C2H6 + 7 O2 -> 4 CO2 + 6 H2O",
        "limiting_reagent": ["C2H6"], "atom_economy_percent": 61.96,
        "e_factor": 0.614, "percent_yield": {"CO2": {"actual_g": 25.0,
        "theoretical_g": 29.271, "percent_yield": 85.41}}, ...}
```

## 最终回复（飞书）

```text
已完成配平与化学计量计算。

配平方程式：2 C2H6 + 7 O2 -> 4 CO2 + 6 H2O
限制试剂：C2H6（乙烷）　　反应进度 ξ：0.1663 mol
原子经济性：61.96%　　理论 E-factor：0.614
CO2 百分产率：85.41%（实际 25.00 g / 理论 29.27 g）

反应物
| 物质 | 系数 | 投料 (g/mol) | 消耗 (g/mol) | 剩余 (g/mol) | 限制? |
|---|---:|---|---|---|:--:|
| C2H6 | 2 | 10.00 / 0.3326 | 10.00 / 0.3326 | 0.00 / 0.0000 | ✅ |
| O2 | 7 | 40.00 / 1.2501 | 37.24 / 1.1640 | 2.76 / 0.0861 | — |

产物（理论产量）
| 物质 | 系数 | 理论产量 (mol) | 理论产量 (g) |
|---|---:|---:|---:|
| CO2 | 4 | 0.6651 | 29.27 |
| H2O | 6 | 0.9977 | 17.97 |

产量对比图：compound_outputs/ethane/yield.png（理论 vs 实际）

说明：乙烷为限制试剂，氧气过量 2.76 g。理论产量为上限值，实际产率 85.41%。
```

（回复时把 `yield.png` 作为图片附件发回飞书。）

## 其他子命令快速示例

```bash
# 配平复杂氧化还原（LLM 易算错，脚本保证正确）
python "<SKILL_DIR>/scripts/stoichiometry_helper.py" balance \
  --equation "KMnO4 + HCl -> KCl + MnCl2 + H2O + Cl2"
#  -> 2 KMnO4 + 16 HCl -> 2 KCl + 2 MnCl2 + 8 H2O + 5 Cl2

# 燃烧分析反推分子式：0.240 g 样品 -> CO2 0.352 g, H2O 0.144 g, M=180.16
python "<SKILL_DIR>/scripts/stoichiometry_helper.py" empirical \
  --combustion "CO2=0.352g, H2O=0.144g" --sample 0.240g --molar-mass 180.16
#  -> 实验式 CH2O，分子式 C6H12O6（葡萄糖）

# 反应焓：甲烷燃烧
python "<SKILL_DIR>/scripts/stoichiometry_helper.py" thermo \
  --equation "CH4 + O2 -> CO2 + H2O"
#  -> ΔH = -890.5 kJ，放热
```
