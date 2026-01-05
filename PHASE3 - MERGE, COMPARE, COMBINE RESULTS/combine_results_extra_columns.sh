#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob

# 1) Per-file processing: keep DCv2 lines, reorder cols, compute percent (1 decimal), write per-TYPE result
for f in *_bcftools_first_*_comparison_results_GTvsGT.txt; do
    GENO="${f%%_bcftools_first_*}"
    TYPE=$(grep -oE 'het|hom' <<<"$f" | head -n1)

    awk -v TYPE="$TYPE" '
        BEGIN {
            OFS = "\t"
            print "mother","embryo",TYPE"_sites_total",TYPE"_sites_mismatch",TYPE"_match_percentage"
        }
        /^DCv2/ {
            perc = ($6 != 0) ? sprintf("%.1f", ($7/$6)*100) : "NA"
            print $2, $3, $6, $4, perc
        }
    ' "$f" > "${GENO}_${TYPE}_results.tsv"
done

# 2) Join het + hom per GENO (strip headers BEFORE sort; full outer join on mother+embryo)
for het in *_het_results.tsv; do
    GENO="${het%_het_results.tsv}"
    hom="${GENO}_hom_results.tsv"
    [[ -f "$hom" ]] || continue

    tail -n +2 "$het" | sort -t$'\t' -k2,2 > het.tmp
    tail -n +2 "$hom" | sort -t$'\t' -k2,2 > hom.tmp

    awk -F'\t' -v OFS='\t' '
        NR==FNR { key=$1 FS $2; data[key]=$0; next }            
        { key=$1 FS $2; seen[key]=1
          if (key in data) { print data[key], $3, $4, $5 }      
          else              { print $1, $2, "NA","NA","NA", $3, $4, $5 }
        }
        END {
          for (k in data) if (!(k in seen)) {
            split(data[k], a, FS)
            print a[1], a[2], a[3], a[4], a[5], "NA","NA","NA"
          }
        }
    ' het.tmp hom.tmp > combined.tmp

    {
      printf "mother\tembryo\thet_sites_total\thet_sites_mismatch\thet_match_percentage\thom_sites_total\thom_sites_mismatch\thom_match_percentage\n"
      cat combined.tmp
    } > "${GENO}_results_combined.tsv"
done

rm -f het.tmp hom.tmp combined.tmp 2>/dev/null || true

# 3) Build the super file with a single header
files=( *_results_combined.tsv )
out="all_results_combined.tsv"
files=( "${files[@]/$out}" )

[[ ${#files[@]} -gt 0 ]] || { echo "No per-GENO combined files found."; exit 0; }

head -n 1 "${files[0]}" > "$out"
for f in "${files[@]}"; do
  [[ -z "$f" ]] && continue
  tail -n +2 "$f" >> "$out"
done

# 4) Add total_sites, total_mismatch, total_match_percentage
awk -F'\t' -v OFS='\t' '
    NR==1 {
        print $0,"total_sites","total_mismatch","total_match_percentage"
        next
    }
    {
        het_tot = ($3 == "NA" ? 0 : $3)
        hom_tot = ($6 == "NA" ? 0 : $6)
        het_mis = ($4 == "NA" ? 0 : $4)
        hom_mis = ($7 == "NA" ? 0 : $7)

        total_sites = het_tot + hom_tot
        total_mismatch = het_mis + hom_mis

        if (total_sites > 0) {
            perc = sprintf("%.1f", ((total_sites - total_mismatch) / total_sites) * 100)
        } else {
            perc = "NA"
        }

        print $0, total_sites, total_mismatch, perc
    }
' "$out" > tmp && mv tmp "$out"
