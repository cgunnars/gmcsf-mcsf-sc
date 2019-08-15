M_D1 = ./data/reads.Sample1_MCSF_R1_001.fastq_bq10_star_corrected.umi.dge.txt
M_D2 = ./data/mcsf_day6_1.txt ./data/mcsf_day6_2.txt
M_FILES = $(M_D1) $(M_D2)
M_DONORS = 1 2 2

G_D1 = ./data/reads.Sample2_gMCSF_R1_001.fastq_bq10_star_corrected.umi.dge.txt
G_D2 = .data/gmscf_day6_2.txt ./data/gmcsf_day6_1.txt
G_FILES = $(G_D1) $(G_D2)
G_DONORS = 1 2 2

load: ./data/raw-mergeruns.rds
.PHONY: ./src/loadData.R ./src/performQC.R $(M_FILES) $(G_FILES)

./data/raw-mergeruns.rds: ./src/loadData.R $(M_FILES) $(G_FILES)
	Rscript $< --mfiles $(M_FILES) --gfiles $(G_FILES) \
	--mdonors $(M_DONORS) --gdonors $(G_DONORS) -o $@

qc: ./data/qc-mergeruns.rds
./data/qc-mergeruns.rds: ./src/performQC.R ./data/raw-mergeruns.rds
	Rscript $< -i ./data/raw-mergeruns.rds -o $@
