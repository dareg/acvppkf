#!/usr/bin/awk -f
/HOOK/ {name = $3 "@" $2;  total[name]+=$4}
END{
	PROCINFO["sorted_in"] = "@ind_str_asc"
	for(k in total){
		print k, total[k]
	}
}
