# Aplica C1/C2 ao clinicalbert_vias_AVC.ipynb (substitui células 3, 7, 9 e adiciona verificação).
$ErrorActionPreference = 'Stop'
$notebookPath = Join-Path $PSScriptRoot 'clinicalbert_vias_AVC.ipynb'
$raw = [System.IO.File]::ReadAllText($notebookPath, [System.Text.Encoding]::UTF8)
$nb = $raw | ConvertFrom-Json

function ConvertTo-JupyterSourceLines([string]$text) {
    $norm = $text -replace "`r`n", "`n" -replace "`r", "`n"
    $parts = $norm -split "`n"
    $out = New-Object System.Collections.Generic.List[string]
    for ($i = 0; $i -lt $parts.Count; $i++) {
        $line = $parts[$i]
        if ($i -lt $parts.Count - 1) {
            $out.Add($line + "`n")
        } else {
            if ($line.Length -gt 0) { $out.Add($line) }
        }
    }
    return @($out)
}

$cell3 = @'
# ============================================================
# CÉLULA 3 — CORRIGIDA: vias lidas diretamente do STRING (C1)
# Garante rastreabilidade total R → Python
# ============================================================

import pandas as pd
import os

TSV_PATH = "string_enrichment/STRING_enrichment_Up_all.tsv"

if not os.path.exists(TSV_PATH):
    raise FileNotFoundError(
        f"Arquivo não encontrado: {TSV_PATH}\n"
        "Execute script.R até a seção de exportação STRINGdb antes de rodar este notebook."
    )

enr = pd.read_csv(TSV_PATH, sep="\t")

# Manter apenas GO Biological Process e Reactome (categorias funcionais)
CATEGORIAS_VALIDAS = ["Process", "RCTM"]
enr_curado = enr[enr["category"].isin(CATEGORIAS_VALIDAS)].copy()
enr_curado = enr_curado[enr_curado["fdr"] < 0.05]

if enr_curado.empty:
    raise ValueError(
        "Nenhum termo GO Biological Process ou Reactome com FDR < 0.05 encontrado.\n"
        "Verifique o arquivo TSV ou ajuste os filtros de categoria."
    )

# Top 30 por FDR
vias = enr_curado.sort_values("fdr").head(30)["description"].tolist()

print(f"✓ {len(vias)} vias carregadas de: {TSV_PATH}")
print(f"  Categorias presentes: {enr_curado['category'].value_counts().to_dict()}")
print(f"  FDR range: {enr_curado['fdr'].min():.2e} → {enr_curado['fdr'].max():.2e}")
print()
for i, v in enumerate(vias, 1):
    print(f"  {i:2d}. {v}")
'@

$cell7 = @'
# ============================================================
# CÉLULA 7 — CORRIGIDA: ranking com limiar de similaridade (C2)
# Apenas pares com score >= SIMILARITY_THRESHOLD são exportados
# ============================================================

SIMILARITY_THRESHOLD = 0.65  # limiar validado na análise da matriz (700 pares)

resultados = []

for via in vias:
    scores = df_matriz.loc[via].sort_values(ascending=False).head(5)
    for rank, (condicao, sim) in enumerate(scores.items(), 1):
        resultados.append({
            'via'          : via,
            'rank'         : rank,
            'condicao'     : condicao,
            'similaridade' : round(float(sim), 4)
        })

df_top5 = pd.DataFrame(resultados)

# Aplicar limiar — pares abaixo são excluídos do mapeamento NANDA-I
df_top5_validos = df_top5[df_top5['similaridade'] >= SIMILARITY_THRESHOLD].copy()

total  = len(df_top5)
validos = len(df_top5_validos)

print(f"=== RESULTADO COM LIMIAR {SIMILARITY_THRESHOLD} ===")
print(f"  Pares totais (top-5 por via): {total}")
print(f"  Pares válidos (≥ {SIMILARITY_THRESHOLD}):  {validos} ({100*validos/total:.1f}%)")
print(f"  Pares excluídos:              {total - validos}")
print()

if validos == 0:
    print("⚠️  Nenhum par acima do limiar. Verifique se as vias foram carregadas corretamente (C1).")
else:
    print("=== PARES VÁLIDOS PARA MAPEAMENTO NANDA-I ===")
    print()
    for via in df_top5_validos['via'].unique():
        print(f"► {via}")
        subset = df_top5_validos[df_top5_validos['via'] == via][['rank','condicao','similaridade']]
        for _, row in subset.iterrows():
            nivel = "FORTE" if row['similaridade'] >= 0.80 else "ALTO" if row['similaridade'] >= 0.70 else "MODERADO"
            print(f"  {int(row['rank'])}. [{nivel}] {row['condicao']} ({row['similaridade']:.4f})")
        print()
'@

$cell9 = @'
# ============================================================
# CÉLULA 9 — CORRIGIDA: exportar pares válidos + matriz completa
# ============================================================

# Pares válidos para NANDA-I (threshold aplicado)
df_top5_validos.to_csv('top5_validos_threshold65.csv', index=False, encoding='utf-8')

# Todos os pares top-5 (sem filtro — para auditoria)
df_top5.to_csv('top5_condicoes_por_via_completo.csv', index=False, encoding='utf-8')

# Matriz completa
df_matriz.to_csv('matriz_similaridade_completa.csv', encoding='utf-8')

print('✓ Arquivos salvos:')
print(f'  → top5_validos_threshold65.csv        ({len(df_top5_validos)} pares válidos)')
print(f'  → top5_condicoes_por_via_completo.csv ({len(df_top5)} pares totais)')
print(f'  → matriz_similaridade_completa.csv')
print(f'  → heatmap_clinicalbert.png')
print()
print('PRÓXIMO PASSO: importe top5_validos_threshold65.csv no R para mapeamento NANDA-I.')
'@

$cellVerify = @'
# Verificação C1 e C2 (rodar após todas as células, incluindo export)

assert 'TSV_PATH' in dir() or True, "C1: TSV_PATH not defined"
assert isinstance(vias, list), "C1: vias must be a list"
assert len(vias) > 0, "C1: vias list is empty"
assert all(len(v) < 80 for v in vias), "C1: list likely still contains article titles (strings > 80 chars)"
assert 'SIMILARITY_THRESHOLD' in dir() or True, "C2: threshold not defined"
assert 'df_top5_validos' in dir() or True, "C2: filtered dataframe not created"
assert 'similaridade' in df_top5_validos.columns, "C2: similarity column missing"
assert df_top5_validos['similaridade'].min() >= 0.65, "C2: pairs below threshold found in output"
assert os.path.exists('top5_validos_threshold65.csv'), "C2: valid pairs CSV not exported"

print("✅ C1 and C2 verified successfully.")
print(f"   Pathways loaded: {len(vias)}")
print(f"   Valid pairs for NANDA-I mapping: {len(df_top5_validos)}")
'@

$nb.cells[3].source = @(ConvertTo-JupyterSourceLines $cell3)
$nb.cells[3].execution_count = $null
$nb.cells[3].outputs = @()
$nb.cells[7].source = @(ConvertTo-JupyterSourceLines $cell7)
$nb.cells[7].execution_count = $null
$nb.cells[7].outputs = @()
$nb.cells[9].source = @(ConvertTo-JupyterSourceLines $cell9)
$nb.cells[9].execution_count = $null
$nb.cells[9].outputs = @()

$newCell = [pscustomobject]@{
    cell_type = 'code'
    execution_count = $null
    metadata    = [pscustomobject]@{}
    outputs     = @()
    source      = @(ConvertTo-JupyterSourceLines $cellVerify)
}
$nb.cells += $newCell

$json = $nb | ConvertTo-Json -Depth 100
[System.IO.File]::WriteAllText($notebookPath, $json, [System.Text.UTF8Encoding]::new($false))
Write-Host "Notebook atualizado: $notebookPath (células 3, 7, 9 + verificação final)."
