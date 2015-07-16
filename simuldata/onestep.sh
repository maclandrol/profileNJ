#!/bin/bash

data=$1
output=data/"$data"_trees.nw
final_output=data/$data/$data.rtrees
./reconcile.py run -s tests/testinput/plantes.nw -g data/$data/$data.nwk --show_branch_tag  -S genetospecie.smap --display_losses --output data/$data/reroot/$data --reroot &&
./profileNJ -s tests/testinput/plantes.nw  -g data/$data/$data.nwk -gl 1 --slimit=-1  --plimit=-1 -r 'best' -S genetospecie.smap -o data/$data/"$data"_20.trees --seuil=20 -d alignements/$data.dist &&
./profileNJ -s tests/testinput/plantes.nw  -g data/$data/$data.nwk -gl 1 --slimit=-1  --plimit=-1 -r 'best' -S genetospecie.smap -o data/$data/"$data"_30.trees --seuil=30 -d alignements/$data.dist &&
./profileNJ -s tests/testinput/plantes.nw  -g data/$data/$data.nwk -gl 1 --slimit=-1  --plimit=-1 -r 'best' -S genetospecie.smap -o data/$data/"$data"_40.trees --seuil=40 -d alignements/$data.dist &&
./profileNJ -s tests/testinput/plantes.nw  -g data/$data/$data.nwk -gl 1 --slimit=-1  --plimit=-1 -r 'best' -S genetospecie.smap -o data/$data/"$data"_50.trees --seuil=50 -d alignements/$data.dist &&
./profileNJ -s tests/testinput/plantes.nw  -g data/$data/$data.nwk -gl 1 --slimit=-1  --plimit=-1 -r 'best' -S genetospecie.smap -o data/$data/"$data"_60.trees --seuil=60 -d alignements/$data.dist &&
./profileNJ -s tests/testinput/plantes.nw  -g data/$data/$data.nwk -gl 1 --slimit=-1  --plimit=-1 -r 'best' -S genetospecie.smap -o data/$data/"$data"_70.trees --seuil=70 -d alignements/$data.dist &&
./profileNJ -s tests/testinput/plantes.nw  -g data/$data/$data.nwk -gl 1 --slimit=-1  --plimit=-1 -r 'best' -S genetospecie.smap -o data/$data/"$data"_80.trees --seuil=80 -d alignements/$data.dist &&
./profileNJ -s tests/testinput/plantes.nw  -g data/$data/$data.nwk -gl 1 --slimit=-1  --plimit=-1 -r 'best' -S genetospecie.smap -o data/$data/"$data"_90.trees --seuil=90 -d alignements/$data.dist &&
./treecat.py data/$data/"$data"_*.trees data/$data/reroot/$data.trees data/"$data"_root.nwk > $output && cat $output | grep -v '>' > "$final_output"