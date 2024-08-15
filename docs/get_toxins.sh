python ../bin/grouping.py \
    -i /home/shannc/Bio_SDD/MUIC_senior_project/workflow/results/C_indra/2-Second_pass/C_indra_all_wcoverage.tsv \
    -o /home/shannc/Bio_SDD/MUIC_senior_project/workflow/results/C_indra/Analysis/toxin_groups.tsv \
    -t /home/shannc/Bio_SDD/MUIC_senior_project/workflow/data/reference/toxin_groups.tsv \
    -m toxin \
    -p /home/shannc/Bio_SDD/MUIC_senior_project/workflow/data/reference/Pfam-A.clans.tsv \
    -a /home/shannc/Bio_SDD/MUIC_senior_project/workflow/config/protein_groups.toml \
    -e /home/shannc/Bio_SDD/MUIC_senior_project/workflow/results/C_indra/Analysis/toxin_groups_evidence.tsv

python ../bin/grouping.py \
    -i /home/shannc/Bio_SDD/MUIC_senior_project/workflow/results/C_indra.calibrated/2-Second_pass/C_indra.calibrated_all_wcoverage.tsv \
    -o /home/shannc/Bio_SDD/MUIC_senior_project/workflow/results/C_indra.calibrated/Analysis/toxin_groups.tsv \
    -t /home/shannc/Bio_SDD/MUIC_senior_project/workflow/data/reference/toxin_groups.tsv \
    -m toxin \
    -p /home/shannc/Bio_SDD/MUIC_senior_project/workflow/data/reference/Pfam-A.clans.tsv \
    -a /home/shannc/Bio_SDD/MUIC_senior_project/workflow/config/protein_groups.toml \
    -e /home/shannc/Bio_SDD/MUIC_senior_project/workflow/results/C_indra.calibrated/Analysis/toxin_groups_evidence.tsv
