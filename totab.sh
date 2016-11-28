for i in salida*txt; do awk '/^K=/ { for (i = 2; i < NF; i++) {printf "%s ",$i; }; print ""}' $i > ${i/salida/tabla}; done
