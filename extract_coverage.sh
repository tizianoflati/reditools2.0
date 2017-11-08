samtools depth $1 | grep -vP "\t0$" | awk '{print $0 > "cov/"$1}'
