# 自编 Agent Skill：reaction-stoichiometry（反应分析）

本目录是 EXP3 作业中**自行编写**的 Agent Skill，用于在 OpenClaw / 飞书中完成
一个与课堂示例不同的化学任务集：**方程式配平 + 化学计量 + 实验式反推 + 热化学**。

## 为什么不是"普通 LLM 也能做"

裸 LLM 在这些任务上**不可靠甚至做不到**，而本 Skill 解决了：

1. **精确配平**：用有理数线性代数求解系数，避免 LLM 在复杂氧化还原中"幻觉"出错误系数。
2. **确定性计算**：限制试剂、理论产量、百分产率、原子经济性 / E-factor、反应焓 ΔH 全部由脚本精确算出，可复现、可审计。
3. **生成真实文件**：可输出产量柱状图（PNG，没有 matplotlib 时退回 SVG）和 JSON 结果文件——文本模型无法产出文件。

核心功能仅用 Python 标准库（`fractions`），不联网、不依赖 RDKit，演示时不会因网络或依赖缺失而失败。

## 四个子命令

| 子命令 | 功能 |
|---|---|
| `balance`   | 配平方程式 + 各物质摩尔质量 |
| `stoich`    | 限制试剂、理论产量、剩余量、百分产率、原子经济性 / E-factor，可选生成产量柱状图 |
| `empirical` | 由燃烧分析或元素质量百分比反推实验式 / 分子式 |
| `thermo`    | 用内置标准生成焓按盖斯定律算反应焓 ΔH |

## 目录结构

```text
skills/
└── reaction-stoichiometry/
    ├── SKILL.md                         # Skill 定义（触发条件、工作流、输出格式）
    ├── scripts/
    │   └── stoichiometry_helper.py      # 配平与化学计量核心脚本（无第三方依赖）
    └── examples/
        └── ethane-combustion.md         # 乙烷燃烧完整示例
```

## 让 OpenClaw 加载

在本目录（`skills/`）或仓库根目录执行其一：

```bash
# 方式 A：Skills CLI
npx skills add ./reaction-stoichiometry -y
npx skills list | grep reaction-stoichiometry

# 方式 B：OpenClaw 旧工作流
openclaw skills install ./reaction-stoichiometry --as reaction-stoichiometry
openclaw skills list | grep reaction-stoichiometry
```

加载后建议新开一个会话，让 skill 快照重新生成。

## 脚本自检

```bash
S=./reaction-stoichiometry/scripts/stoichiometry_helper.py
python "$S" balance   --equation "KMnO4 + HCl -> KCl + MnCl2 + H2O + Cl2"
python "$S" stoich    --equation "C2H6 + O2 -> CO2 + H2O" --given "C2H6=10g, O2=40g" --actual "CO2=25g" --plot out/yield.png
python "$S" empirical --combustion "CO2=0.352g, H2O=0.144g" --sample 0.240g --molar-mass 180.16
python "$S" thermo    --equation "CH4 + O2 -> CO2 + H2O"
```

任一打印出 `[RESULT] ...` 与 `[JSON] ...` 即表示脚本可用
（`python` 可替换为 `python3` 或 `uv run`）。

## 在飞书中的用法

直接用自然语言描述任务即可，无需说出 Skill 名称，例如：

```text
帮我配平 KMnO4 + HCl -> KCl + MnCl2 + H2O + Cl2，并给出各物质摩尔质量。
```

```text
合成氨 N2 + H2 -> NH3，我有 14 g 氮气和 5 g 氢气，谁是限制试剂？最多生成多少克氨？
```

```text
乙烷燃烧 C2H6 + O2 -> CO2 + H2O，投 10 g 乙烷和 40 g 氧气，实际得 25 g CO2，
求百分产率和原子经济性，并画产量对比图。
```

```text
某有机物 0.240 g 完全燃烧得 0.352 g CO2 和 0.144 g H2O，摩尔质量 180，求分子式。
```

```text
甲烷燃烧 CH4 + O2 -> CO2 + H2O 是放热还是吸热？反应焓是多少？
```
