
q=Kuenenia_stuttgartiensis_fromNCBI_chromosome.faa
db=all_nitrogen_methane_genes.pep
/90days/s4506266/softwares/diamond-v0.9.29.130/diamond blastp -q Kuenenia_stuttgartiensis_fromNCBI_chromosome.faa --db $db --outfmt 6 --strand both --max-target-seqs 10 -o Kuenenia_check_nitrogen_metabolism.diamond_fmt6
echo done

