#!/bin/bash

file=~/issm/share_setup/Profiles_bash.txt

IFS=$'\r\n' GLOBIGNORE='*' command eval  'points=($(cat $file))'

#for i in ${!points[@]}"
#do
for i in $(seq 0 4 $((${#points[@]}-1)));
do
	x1=${points[i]}
	x2=${points[i+1]}
	y1=${points[i+2]}
	y2=${points[i+3]}

	r.profile -g input=BedMachine_clip coordinates=$x1,$y1,$x2,$y2 >> all_profiles.txt
	echo new_profile >> all_profiles.txt
done	


