#!/usr/bin/env bash
set -u
cd "$(dirname "$0")/.."
for i in $(seq 1 480); do
  n=$(find results/drem_bronchial/models -name tf_ranking.csv | wc -l | tr -d ' ')
  if [ "$n" -ge 32 ]; then break; fi
  sleep 60
done
echo "rankings present: $n"
BTC_DREM_SET=drem_bronchial NUMBA_CACHE_DIR="$PWD/.numba_cache" python scripts/collect_drem_ali.py 2>&1 | tail -4
python scripts/compare_drem_viper.py --drem-set drem_bronchial --label "bronchial GSE233145" 2>&1 | tail -40
