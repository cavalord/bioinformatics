#!/usr/bin/env bash
set -euo pipefail
 
# ==============================================================================
# Script de preparação e envio de VCF de exoma (WES) para imputação genômica
# via TOPMed Imputation Server (GRCh38)
#
# Pré-requisitos:
#   - bcftools, tabix, bgzip (pacote htslib) e curl instalados e no PATH
#   - Um FASTA de referência GRCh38 (ex.: GRCh38_full_analysis_set_plus_decoy_hla.fa)
#   - Uma conta na TOPMed Imputation Server + token de API
#     (https://imputation.biodatacatalyst.nhlbi.nih.gov -> username -> Profile -> API Tokens)
#     O token expira em 30 dias.
#
# Uso:
#   1. Edite as variáveis na seção CONFIGURAÇÕES abaixo
#   2. export TOPMED_API_TOKEN="seu-token-aqui"
#   3. bash imputacao_topmed.sh prepare   # normaliza, corrige build/nomenclatura, separa por cromossomo
#   4. bash imputacao_topmed.sh submit    # envia os arquivos via API
#
# Observação importante (dados de exoma/WES):
#   A imputação preenche posições ausentes usando como "âncoras" apenas as
#   variantes presentes no seu VCF. Como o exoma só cobre regiões exônicas,
#   a qualidade da imputação (R2) tende a cair MUITO à medida que a posição-alvo
#   se afasta de uma região capturada. Espere R2 baixo para SNPs intrônicos/
#   intergênicos distantes de éxons — sempre filtre a saída por R2 antes de
#   usar os dados no pgsc_calc.
# ==============================================================================
 
# ---------------------- CONFIGURAÇÕES (edite antes de rodar) ----------------------
INPUT_VCF="/caminho/para/exoma_coorte.vcf.gz"     # VCF original da coorte (WES), multi-cromossomo
REF_FASTA="/caminho/para/GRCh38_full_analysis_set_plus_decoy_hla.fa"
OUTDIR="/caminho/para/saida_imputacao"
THREADS=8
 
# Parâmetros do job na TOPMed Imputation Server
API_URL="https://imputation.biodatacatalyst.nhlbi.nih.gov/api/v2"
API_TOKEN="${TOPMED_API_TOKEN:-}"          # defina via: export TOPMED_API_TOKEN=... (não deixe o token hardcoded aqui)
JOB_NAME="minha_coorte_exoma_$(date +%Y%m%d)"
REFPANEL="apps@topmed-r3"
BUILD="hg38"
PHASING="eagle"
POPULATION="all"                            # "all" é o mais seguro se a coorte tem ancestralidade mista/desconhecida
R2_FILTER="0.3"                             # descarta variantes imputadas de baixa qualidade já na saída do servidor
MODE="imputation"                           # use "qconly" para rodar só o QC sem gastar a cota de imputação
 
mkdir -p "$OUTDIR"/{norm,split,logs}
 
# ---------------------- Funções ----------------------
 
check_deps() {
  for tool in bcftools tabix bgzip curl; do
    command -v "$tool" >/dev/null 2>&1 || { echo "ERRO: '$tool' não encontrado no PATH."; exit 1; }
  done
}
 
prepare() {
  check_deps
 
  echo "[1/4] Normalizando VCF (split de multialélicos, left-align de indels, checagem contra o FASTA de referência)..."
  bcftools norm \
    -m -any \
    -f "$REF_FASTA" \
    --threads "$THREADS" \
    -Oz -o "$OUTDIR/norm/coorte.norm.vcf.gz" \
    "$INPUT_VCF" 2> "$OUTDIR/logs/01_norm.log"
  tabix -f -p vcf "$OUTDIR/norm/coorte.norm.vcf.gz"
 
  local FINAL_VCF="$OUTDIR/norm/coorte.norm.vcf.gz"
 
  echo "[2/4] Ajustando nomenclatura de cromossomos para o padrão 'chr' (exigido para build hg38)..."
  if bcftools view -h "$FINAL_VCF" | grep -q '^##contig=<ID=1,'; then
    for i in $(seq 1 22) X Y MT; do printf '%s\tchr%s\n' "$i" "$i"; done > "$OUTDIR/chr_rename_map.txt"
    bcftools annotate \
      --rename-chrs "$OUTDIR/chr_rename_map.txt" \
      --threads "$THREADS" \
      -Oz -o "$OUTDIR/norm/coorte.chr.vcf.gz" \
      "$FINAL_VCF" 2>> "$OUTDIR/logs/02_rename.log"
    tabix -f -p vcf "$OUTDIR/norm/coorte.chr.vcf.gz"
    FINAL_VCF="$OUTDIR/norm/coorte.chr.vcf.gz"
  else
    echo "  -> Prefixo 'chr' já presente, nenhuma renomeação necessária."
  fi
 
  echo "[3/4] Garantindo header em VCF v4.2 (a TOPMed Imputation Server rejeita 4.3)..."
  bcftools view -h "$FINAL_VCF" | sed 's/##fileformat=VCFv4\.3/##fileformat=VCFv4.2/' > "$OUTDIR/norm/new_header.txt"
  bcftools reheader -h "$OUTDIR/norm/new_header.txt" -o "$OUTDIR/norm/coorte.final.vcf.gz" "$FINAL_VCF"
  tabix -f -p vcf "$OUTDIR/norm/coorte.final.vcf.gz"
  FINAL_VCF="$OUTDIR/norm/coorte.final.vcf.gz"
 
  echo "[4/4] Separando por cromossomo (autossomos 1-22; adicione chrX manualmente se for usar)..."
  for chr in $(seq 1 22); do
    bcftools view -r "chr${chr}" "$FINAL_VCF" \
      -Oz -o "$OUTDIR/split/coorte.chr${chr}.vcf.gz" \
      2>> "$OUTDIR/logs/03_split.log"
    tabix -f -p vcf "$OUTDIR/split/coorte.chr${chr}.vcf.gz"
  done
 
  bcftools stats "$FINAL_VCF" > "$OUTDIR/logs/04_bcftools_stats.txt"
 
  echo ""
  echo "Pronto! Arquivos por cromossomo em: $OUTDIR/split/"
  echo "Estatísticas de QC em: $OUTDIR/logs/04_bcftools_stats.txt"
}
 
submit() {
  check_deps
  if [[ -z "$API_TOKEN" ]]; then
    echo "ERRO: defina seu token com  export TOPMED_API_TOKEN='seu-token-aqui'  antes de rodar 'submit'."
    exit 1
  fi
 
  local files_args=()
  for f in "$OUTDIR"/split/coorte.chr*.vcf.gz; do
    files_args+=(-F "files=@${f}")
  done
 
  if [[ ${#files_args[@]} -eq 0 ]]; then
    echo "ERRO: nenhum arquivo encontrado em $OUTDIR/split/ — rode '$0 prepare' primeiro."
    exit 1
  fi
 
  echo "Enviando ${#files_args[@]} arquivos para a TOPMed Imputation Server..."
  curl "$API_URL/jobs/submit/imputationserver" \
    -X "POST" \
    -H "X-Auth-Token: ${API_TOKEN}" \
    -F "job-name=${JOB_NAME}" \
    "${files_args[@]}" \
    -F "refpanel=${REFPANEL}" \
    -F "build=${BUILD}" \
    -F "phasing=${PHASING}" \
    -F "population=${POPULATION}" \
    -F "r2Filter=${R2_FILTER}" \
    -F "mode=${MODE}" \
    -F "meta=yes" \
    | tee "$OUTDIR/logs/05_submit_response.json"
 
  echo ""
  echo "Job enviado. Acompanhe o status em: https://imputation.biodatacatalyst.nhlbi.nih.gov"
  echo "A senha para descriptografar os resultados chega por e-mail — guarde-a, ela NÃO é reenviada em caso de perda."
}
 
# ---------------------- Execução ----------------------
case "${1:-}" in
  prepare) prepare ;;
  submit)  submit ;;
  *)
    echo "Uso: $0 {prepare|submit}"
    exit 1
    ;;
esac
 
