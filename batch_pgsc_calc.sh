#!/usr/bin/env bash
#
# Roda o pgsc_calc para as listas de PGS IDs restantes, criando um
# diretório de saída "loteN" para cada uma.
#
# Mapeamento assumido a partir do script original:
#   lista_prs_0.txt  -> lote1  (já rodado)
#   lista_prs_1.txt  -> lote2  (já rodado)
#   lista_prs_2.txt  -> lote3
#   ...
#   lista_prs_69.txt -> lote70
#
# Se o seu mapeamento for diferente (ex: lista_prs_N.txt -> loteN),
# troque a linha "lote=$((i + 1))" abaixo por "lote=${i}".

set -uo pipefail

BASE="/home/venus/mar/vtarghetta/wxs_analysis/pgscalc"
NEXTFLOW="/home/venus/mar/vtarghetta/nextflow"
LOGDIR="${BASE}/logs"
mkdir -p "${LOGDIR}"

# já rodamos lista_prs_0.txt e lista_prs_1.txt -> começa em 2
for i in $(seq 2 69); do
    lote=$((i + 1))
    outdir="${BASE}/phs000980/lote${lote}"
    listfile="${BASE}/lista_prs_${i}.txt"
    log="${LOGDIR}/lote${lote}.log"

    echo "=== Rodando lista_prs_${i}.txt -> lote${lote} ==="

    if [[ ! -f "${listfile}" ]]; then
        echo "AVISO: ${listfile} não encontrado, pulando lote${lote}."
        continue
    fi

    mkdir -p "${outdir}"

    "${NEXTFLOW}" run pgscatalog/pgsc_calc \
        -profile docker \
        --input "${BASE}/samplesheet-pgscalc.csv" \
        --target_build GRCh38 \
        --pgs_id $(cat "${listfile}") \
        --run_ancestry "${BASE}/pgsc_HGDP+1kGP_v1.tar.zst" \
        --genotypes_cache ./cache \
        --outdir "${outdir}/" \
        --max_cpus 16 --max_memory 120.GB --max_time 240.h \
        --min_overlap 0.0 \
        -resume 2>&1 | tee "${log}"

    status=${PIPESTATUS[0]}
    if [[ ${status} -ne 0 ]]; then
        echo "ERRO: lote${lote} falhou (status ${status}). Veja ${log}."
        echo "Continuando para o próximo lote..."
    else
        echo "OK: lote${lote} concluído."
    fi
done

echo "Loop finalizado."
