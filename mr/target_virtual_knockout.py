# -*- coding: utf-8 -*-
"""
target_virtual_knockout.py — 靶点网络虚拟剔除分析（in silico knockout）
================================================================
设计：
  1. 构建 CKM-血管模块 PPI 网络（STRING v12，combined_score>=400）：
     利钠肽通路 8 基因 + 文献锚定的 CVD/炎症核心基因
  2. 虚拟剔除 = 逐一移除候选靶点节点，度量：
     - 全局效率 (global efficiency) 降幅 ΔE%
     - 最大连通分量 (giant component) 缩减
     - 度/中介中心性基线
  3. 置换检验：与网络中所有其他单节点移除比较，给 z 值与经验 p 值
数据源：STRING-DB 公开 API（真实数据，非模拟）
输出：virtual_knockout_results.csv / virtual_knockout_results.txt
"""
import urllib.request, urllib.parse, json, csv, time, math
import networkx as nx

DIR = r"C:\Users\xuyin\WorkBuddy\2026-09-11-14-07-58\CD特刊选题评估"
SPECIES = 9606
SCORE_CUT = 400  # STRING medium confidence

# ---- CKM-血管模块基因集（文献锚定，Methods 可复述）----
NATRIURETIC = ["NPPA", "NPPB", "NPPC", "NPR1", "NPR2", "NPR3", "MME", "CORIN"]
CVD_CORE = [  # 血压/血管内皮/脂质/炎症核心（AHA CKM 与 CAD GWAS 反复命中）
    "ACE", "AGT", "REN", "EDN1", "NOS3", "LPA", "PCSK9", "APOB", "SORT1",
    "IL6R", "CRP", "TNF", "IL1B", "VCAM1", "ICAM1", "SELE", "MMP9",
    "CCL2", "LEP", "ADIPOQ", "INSR", "PPARG", "ABCG1", "PLAT", "VWF",
]
GENES = sorted(set(NATRIURETIC + CVD_CORE))
TARGETS = ["NPR3", "NPPA", "NPPB", "NPR1", "NPR2", "MME", "CORIN", "IL6R"]

def fetch_string_network(genes, score_cut=400):
    ids = "%0d".join([])  # placeholder
    url = ("https://string-db.org/api/tsv/network?identifiers="
           + "%0d".join(genes) + f"&species={SPECIES}&required_score={score_cut}")
    req = urllib.request.Request(url, headers={"User-Agent": "academic-analysis/1.0"})
    rows = []
    for attempt in range(3):
        try:
            with urllib.request.urlopen(req, timeout=120) as r:
                text = r.read().decode("utf-8", "ignore")
            lines = text.strip().split("\n")
            header = lines[0].split("\t")
            for ln in lines[1:]:
                parts = ln.split("\t")
                if len(parts) == len(header):
                    rows.append(dict(zip(header, parts)))
            return rows
        except Exception as e:
            print(f"  STRING fetch attempt {attempt+1} fail: {e}")
            time.sleep(10)
    return rows

def build_graph(rows, genes):
    G = nx.Graph()
    G.add_nodes_from(genes)
    seen = set()
    for r in rows:
        a, b = r["preferredName_A"], r["preferredName_B"]
        s = float(r["score"])
        if a in genes and b in genes and (a, b) not in seen and (b, a) not in seen:
            seen.add((a, b))
            G.add_edge(a, b, weight=s)
    return G

print("=" * 70)
print("靶点网络虚拟剔除分析（STRING v12, combined score >= %d）" % SCORE_CUT)
print("=" * 70)

rows = fetch_string_network(GENES, SCORE_CUT)
print(f"STRING 返回边记录: {len(rows)}")
G = build_graph(rows, GENES)
n_nodes, n_edges = G.number_of_nodes(), G.number_of_edges()
print(f"CKM-血管模块网络: {n_nodes} 节点 / {n_edges} 边")
print(f"候选靶点度数基线: " + ", ".join(f"{t}={G.degree(t)}" for t in TARGETS))

# ---- 基线中心性 ----
btw = nx.betweenness_centrality(G)
deg = dict(G.degree())
E0 = nx.global_efficiency(G)
gc0 = len(max(nx.connected_components(G), key=len)) if n_nodes else 0
print(f"基线全局效率 E0 = {E0:.4f}；最大连通分量 = {gc0}/{n_nodes}")

# ---- 虚拟剔除：移除单节点后的全局效率 ----
def knockout_effect(g, node):
    H = g.copy()
    H.remove_node(node)
    E = nx.global_efficiency(H)
    gc = len(max(nx.connected_components(H), key=len)) if H.number_of_nodes() else 0
    return E, gc

# 全部节点的剔除效应（完整枚举，n≈30，无需抽样）
all_effects = {}
for node in list(G.nodes()):
    E, gc = knockout_effect(G, node)
    all_effects[node] = {"E": E, "dE_pct": 100 * (E0 - E) / E0, "gc": gc}

dE_values = [v["dE_pct"] for v in all_effects.values()]
mean_dE = sum(dE_values) / len(dE_values)
sd_dE = math.sqrt(sum((x - mean_dE) ** 2 for x in dE_values) / (len(dE_values) - 1))

# ---- 输出 ----
out_rows = []
print("\n%-8s %5s %10s %12s %10s %8s %8s" % ("靶点", "度", "中介性", "ΔE(%)", "z(ΔE)", "p秩", "GC损失"))
print("-" * 70)
for t in sorted(TARGETS, key=lambda x: -all_effects[x]["dE_pct"]):
    v = all_effects[t]
    z = (v["dE_pct"] - mean_dE) / sd_dE if sd_dE > 0 else float("nan")
    # 经验 p：网络中等或更高影响的节点比例
    n_ge = sum(1 for x in dE_values if x >= v["dE_pct"])
    p_rank = n_ge / len(dE_values)
    gc_loss = gc0 - v["gc"]
    print("%-8s %5d %10.4f %12.2f %10.2f %8.3f %8d" % (t, deg[t], btw[t], v["dE_pct"], z, p_rank, gc_loss))
    out_rows.append({
        "target": t, "degree": deg[t], "betweenness": round(btw[t], 4),
        "knockout_dE_pct": round(v["dE_pct"], 2), "z_vs_all_nodes": round(z, 2),
        "p_rank": round(p_rank, 3), "giant_component_loss": gc_loss,
        "role": "candidate" if t in ("NPR3", "NPPA", "NPPB", "NPR1", "NPR2") else ("positive-ctrl" if t == "MME" else ("repurposing" if t == "IL6R" else "candidate")),
    })

# 对照基因（非靶点核心基因）里前 5 名，供比较
print("\n网络中影响力最高的 5 个非靶点基因（对照参考）:")
for node in sorted(all_effects, key=lambda x: -all_effects[x]["dE_pct"]):
    if node not in TARGETS:
        v = all_effects[node]
        print("  %-8s 度=%d  ΔE=%.2f%%  GC损失=%d" % (node, deg[node], v["dE_pct"], gc0 - v["gc"]))
        if sum(1 for n2 in all_effects if n2 not in TARGETS and all_effects[n2]["dE_pct"] >= v["dE_pct"]) >= 5:
            break

# 保存
with open(DIR + r"\virtual_knockout_results.csv", "w", newline="", encoding="utf-8-sig") as f:
    w = csv.DictWriter(f, fieldnames=list(out_rows[0].keys()))
    w.writeheader()
    w.writerows(out_rows)

with open(DIR + r"\virtual_knockout_results.txt", "w", encoding="utf-8") as f:
    f.write("靶点网络虚拟剔除分析（STRING v12, score>=%d）\n" % SCORE_CUT)
    f.write(f"模块网络: {n_nodes} 节点 / {n_edges} 边；基线全局效率 {E0:.4f}\n")
    f.write(f"全节点剔除 ΔE 均值 {mean_dE:.2f}% ± {sd_dE:.2f}%\n\n")
    for r in out_rows:
        f.write(str(r) + "\n")
    f.write("\n全部节点剔除效应（降序）:\n")
    for node in sorted(all_effects, key=lambda x: -all_effects[x]["dE_pct"]):
        f.write("  %-8s ΔE=%.2f%%  GC损失=%d\n" % (node, all_effects[node]["dE_pct"], gc0 - all_effects[node]["gc"]))

# 保存网络边表
with open(DIR + r"\string_network_edges.csv", "w", newline="", encoding="utf-8-sig") as f:
    w = csv.writer(f)
    w.writerow(["source", "target", "string_score"])
    for a, b, d in G.edges(data=True):
        w.writerow([a, b, round(d["weight"], 3)])

print("\nsaved: virtual_knockout_results.csv / .txt / string_network_edges.csv")
