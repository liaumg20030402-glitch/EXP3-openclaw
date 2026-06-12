---
name: reaction-stoichiometry
description: >
  反应分析工具：当用户给出化学反应方程式（可未配平）或化学组成数据，并要求以下任一任务时使用 ——
  (1) 配平方程式、计算各物质摩尔质量；
  (2) 根据反应物用量判断限制试剂、计算理论产量、过量剩余量，给定实际产物质量时计算百分产率，
      并计算原子经济性 / E-factor 等绿色化学指标，可选地生成产量柱状图（PNG/SVG 文件）；
  (3) 由燃烧分析数据（CO2/H2O 质量 + 样品质量）或元素质量百分比反推实验式与分子式；
  (4) 用内置标准生成焓数据按盖斯定律计算反应焓变 ΔH。
  应根据用户描述自动推断任务，用户无需说出 Skill 名称。
  不要用于安全/SDS/GHS（用 chemical-safety-report）、分子3D结构（用 compound-3d-profile）、
  量子化学几何优化（用 xtb-demo-workflow），也不要用于反应机理推断或合成路线设计。
compatibility: >
  核心功能仅依赖 Python 标准库（fractions 精确有理数配平），无需联网、无需 RDKit。
  生成图表时优先用 matplotlib 输出 PNG；若不可用则自动退回纯标准库 SVG。
  可用 `python`、`python3` 或 `uv run` 任一方式执行随附脚本。
license: MIT
skill_type: primary
teaching_role: teaching-example
metadata:
  author: user
  version: '0.2'
---

# 反应分析（配平 / 化学计量 / 实验式 / 热化学）

这是一个**自包含的计算型 Skill**，把化学反应与组成数据转化为可靠、可复现的定量结果。
它做的是**普通 LLM 单独做不可靠、甚至做不到**的事：

1. **精确配平**：用有理数线性代数求解配平系数，避免 LLM 在复杂氧化还原中"幻觉"系数。
2. **确定性计算**：限制试剂、理论产量、产率、原子经济性、ΔH 全部由脚本精确算出。
3. **生成真实文件**：可输出产量柱状图（PNG/SVG）与 JSON 结果文件——文本模型无法生成文件。

核心逻辑封装在随附脚本 `scripts/stoichiometry_helper.py`。Agent 只负责解析用户意图、
调用脚本、把 `[JSON]` 结果整理成中文报告。用户用自然语言描述任务即可。

## 运行环境与脚本调用

脚本仅核心用标准库（matplotlib 可选）。按优先级选可用解释器：

```bash
uv run "<SKILL_DIR>/scripts/stoichiometry_helper.py" <subcommand> ...
python "<SKILL_DIR>/scripts/stoichiometry_helper.py" <subcommand> ...
python3 "<SKILL_DIR>/scripts/stoichiometry_helper.py" <subcommand> ...
```

`<SKILL_DIR>` 是本 `SKILL.md` 所在目录。本机（Windows + Anaconda）绝对路径示例：

```bash
python "E:/硕士课程/EXP3-openclaw/EXP3-openclaw/skills/reaction-stoichiometry/scripts/stoichiometry_helper.py" balance --equation "C2H6 + O2 -> CO2 + H2O"
```

脚本成功时打印若干 `[RESULT] key=value` 行和**一行以 `[JSON]` 开头的完整结果**；
解析 `[JSON]` 即可拿到全部数据。出错时打印 `[ERROR] ...` 并返回非零退出码。

## 四个子命令与触发判断

| 用户意图 | 子命令 |
|---|---|
| 只要配平 / 摩尔质量 | `balance` |
| 给了反应物用量，问限制试剂 / 产量 / 产率 / 绿色指标 / 要图 | `stoich` |
| 给了燃烧数据或元素百分比，求化学式 | `empirical` |
| 问反应放热/吸热、反应焓 ΔH | `thermo` |

### 1. balance —— 配平与摩尔质量

```bash
python "<SKILL_DIR>/scripts/stoichiometry_helper.py" balance \
  --equation "KMnO4 + HCl -> KCl + MnCl2 + H2O + Cl2"
```

读取 `balanced`、`reactant_coefficients`、`product_coefficients`、`species[*].molar_mass`。

### 2. stoich —— 化学计量 + 绿色化学指标 + 可选图表

```bash
python "<SKILL_DIR>/scripts/stoichiometry_helper.py" stoich \
  --equation "C2H6 + O2 -> CO2 + H2O" \
  --given "C2H6=10g, O2=40g" \
  --actual "CO2=25g" \
  --target CO2 \
  --plot "<output_dir>/yield.png"
```

- `--given`（必需）：反应物用量，单位 `g/mg/kg/mol/mmol`。
- `--actual`（可选）：某产物实际获得量 → 计算**百分产率** `percent_yield`。
- `--target`（可选）：目标产物，用于**原子经济性** `atom_economy_percent`（默认取第一个产物）。
- `--plot`（可选）：生成产量柱状图，返回 `chart_file` 绝对路径；有 matplotlib 出 PNG，否则 SVG。

脚本以「反应进度 ξ = min(给定反应物 mol / 系数)」确定限制试剂，并算出理论产量、剩余量、
原子经济性与理论 E-factor。读取 `limiting_reagent`、`products[*].theoretical_g`、
`reactants[*].leftover_g`、`atom_economy_percent`、`e_factor`、`percent_yield`、`chart_file`。

### 3. empirical —— 实验式 / 分子式反推

由燃烧分析（含 C、H、可能含 O 的有机物）：

```bash
python "<SKILL_DIR>/scripts/stoichiometry_helper.py" empirical \
  --combustion "CO2=0.352g, H2O=0.144g" --sample 0.240g --molar-mass 180.16
```

或由元素质量百分比：

```bash
python "<SKILL_DIR>/scripts/stoichiometry_helper.py" empirical \
  --percent "C=40.0, H=6.7, O=53.3" --molar-mass 180.16
```

`--molar-mass`（可选）：给定摩尔质量则进一步给出**分子式** `molecular_formula`。
读取 `empirical_formula`、`empirical_formula_mass`、`molecular_formula`。
注意：燃烧法中 O 由"样品质量 − C、H 质量"差减得出，需提供 `--sample`。

### 4. thermo —— 反应焓 ΔH（盖斯定律）

```bash
python "<SKILL_DIR>/scripts/stoichiometry_helper.py" thermo \
  --equation "CH4 + O2 -> CO2 + H2O"
```

脚本先配平，再用内置标准生成焓表算 ΔH_rxn = Σν·ΔHf(产物) − Σν·ΔHf(反应物)（kJ）。
读取 `dH_rxn_kJ_per_mol_rxn`、`nature`（放热/吸热）、`terms`。
若某物种不在内置表中，脚本会报错并提示用 `--dhf "物种=数值"`（kJ/mol）补充：

```bash
python "<SKILL_DIR>/scripts/stoichiometry_helper.py" thermo \
  --equation "C2H4 + H2 -> C2H6" --dhf "C2H4=52.4, C2H6=-84.0"
```

⚠️ 内置 ΔHf° 假定常规标准态（H2O 取液态；金属/氧化物/盐取固态；气体按写法）。
如用户的状态不同，需用 `--dhf` 覆盖，并在报告中注明假设。

## 化学式书写规则（传给脚本时）

- 元素符号区分大小写：`Co`(钴) ≠ `CO`(一氧化碳)。
- 括号/方括号与下标：`Ca(OH)2`、`Fe2(SO4)3` 支持。
- 水合物点号：`CuSO4·5H2O`（`·`、`•`、`*`、`.` 均可）。
- 状态符 `(s)/(l)/(g)/(aq)` 与离子电荷会被忽略，不影响计算。
- 用户用中文名（如「乙烷」「葡萄糖」）时，先在回复中确认化学式再传脚本。

## 最终回复格式

根据子命令选用对应模板，**数值全部取自脚本输出，不要手算改动**。

### balance

```text
已完成配平。
配平结果：<balanced>

| 物质 | 系数 | 摩尔质量 (g/mol) |
|---|---:|---:|
```

### stoich

```text
已完成配平与化学计量计算。

配平方程式：<balanced>
限制试剂：<limiting_reagent>　　反应进度 ξ：<extent> mol
原子经济性：<atom_economy_percent>%　　理论 E-factor：<e_factor>
（如有实际产量）百分产率：<species> <percent_yield>%

反应物
| 物质 | 系数 | 投料 (g/mol) | 消耗 (g/mol) | 剩余 (g/mol) | 限制? |
|---|---:|---|---|---|:--:|

产物（理论产量）
| 物质 | 系数 | 理论产量 (mol) | 理论产量 (g) |
|---|---:|---:|---:|

（如生成图表）产量柱状图：<chart_file>
说明：理论产量为完全反应、无损失的上限值。
```

若生成了 `chart_file`，应把该图片文件作为附件/图片发回飞书。

### empirical

```text
实验式：<empirical_formula>（式量 <empirical_formula_mass> g/mol）
（如给摩尔质量）分子式：<molecular_formula>（倍数 ×<multiplier>）
各元素物质的量：...
```

### thermo

```text
配平方程式：<balanced>
反应焓 ΔH = <dH_rxn_kJ_per_mol_rxn> kJ（每摩尔反应），<nature>
各物种 ΔHf° 贡献：...
假设：<assumptions>
```

## 错误处理与恢复

| 现象 | 原因 | 恢复操作 |
|---|---|---|
| `Unknown element symbol` | 化学式拼写/大小写错误 | 提示核对，如 `co`→`Co`/`CO` |
| `cannot be balanced` | 两边元素无法守恒 | 提示反应式可能写错或缺物种 |
| `underdetermined (N independent reactions)` | 含多个独立反应 | 请用户拆成单个净反应 |
| `Given species ... not found` | 用量化学式与方程式写法不一致 | 用与方程式完全相同的写法 |
| `Bad amount value` | 用量格式错误 | 用 `物种=数值单位`，如 `O2=40g` |
| `No standard formation enthalpy for ...` | 物种不在内置 ΔHf 表 | 用 `--dhf "物种=数值"` 补充 |
| `[WARN] matplotlib unavailable` | 无 matplotlib | 脚本自动改出 SVG，照常返回 `chart_file` |

恢复原则：脚本报错时**先修正输入再重试，绝不自行手算结果**；可交付的部分结果照常交付。

## 智能体检查清单

- [ ] 从提示推断子命令（balance / stoich / empirical / thermo），不要求用户说 Skill 名。
- [ ] 中文物质名先转化学式并在回复中确认。
- [ ] 配平、产量、产率、ΔH 等数值全部取自脚本 `[JSON]`，不手算。
- [ ] stoich 任务至少给一个反应物用量；给了实际产量才报产率。
- [ ] 生成了 `chart_file` 时把图片发回飞书。
- [ ] thermo 报告注明 ΔHf° 标准态假设。
- [ ] 报告中标注理论产量是上限值。
- [ ] 脚本报错先修正输入再重试。

## 完整示例

见 `examples/` 目录：

- [ethane-combustion.md](examples/ethane-combustion.md)：乙烷燃烧的端到端示例 ——
  配平、限制试剂、理论产量、百分产率、原子经济性，并生成产量柱状图（PNG）。
