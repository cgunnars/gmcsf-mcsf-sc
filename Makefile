M_D1 = data/reads.Sample1_MCSF_R1_001.fastq_bq10_star_corrected.umi.dge.txt
M_D2 = data/mcsf_day6_1.txt data/mcsf_day6_2.txt
M_FILES = $(M_D1) $(M_D2)
M_DONORS = 1 2 2

G_D1 = data/reads.Sample2_gMCSF_R1_001.fastq_bq10_star_corrected.umi.dge.txt
G_D2 = data/gmscf_day6_2.txt data/gmcsf_day6_1.txt
G_FILES = $(G_D1) $(G_D2)
G_DONORS = 1 2 2

RAW_DATA: $(M_FILES) $(G_FILES)

load: data/raw-mergeruns.Rds
#.PHONY: ./src/loadData.R ./src/performQC.R $(M_FILES) $(G_FILES)

data/raw-mergeruns.Rds: src/loadData.R $(RAW_DATA)
	Rscript $< --mfiles $(M_FILES) --gfiles $(G_FILES) \
	--mdonors $(M_DONORS) --gdonors $(G_DONORS) -o $@

qc: data/qc-mergeruns.Rds
data/qc-mergeruns.Rds: src/performQC.R data/raw-mergeruns.Rds
	Rscript $< -i data/raw-mergeruns.rds -o $@

markers: data/markers.txt data/markers.Rds
data/markers.txt: src/markers.R data/qc-mergeruns.Rds
	Rscript $< -i data/qc-mergeruns.Rds -o data/markers --quick
data/markers.Rds: 
	@if test -f $@; then :; else \
		rm -f $<; \
		make $<; \
	fi	

fig/markers-stim.svg: src/plotMarkers.R data/qc-mergeruns.Rds
	Rscript $< -i data/qc-mergeruns.Rds -m data/markers.txt -o fig/markers
fig/markers-donor.svg: src/plotMarkers.R data/qc-mergeruns.Rds
	@if test -f $@; then :; else \
		rm -f $<; \
		make $<; \
	fi
