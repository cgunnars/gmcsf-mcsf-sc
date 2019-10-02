M_D1 = data/raw/reads.Sample1_MCSF_R1_001.fastq_bq10_star_corrected.umi.dge.txt
M_D2 = data/raw/mcsf_day6_1.txt data/raw/mcsf_day6_2.txt
M_FILES = $(M_D1) $(M_D2)
M_DONORS = 1 2 2

G_D1 = data/raw/reads.Sample2_gMCSF_R1_001.fastq_bq10_star_corrected.umi.dge.txt
G_D2 = data/raw/gmscf_day6_2.txt data/raw/gmcsf_day6_1.txt
G_FILES = $(G_D1) $(G_D2)
G_DONORS = 1 2 2

$(RAW_DATA): $(M_FILES) $(G_FILES)

load: data/raw/raw-mergeruns.Rds
#.PHONY: ./src/loadData.R ./src/performQC.R $(M_FILES) $(G_FILES)

data/raw/raw-mergeruns.Rds: src/loadData.R $(RAW_DATA)
	Rscript $< --mfiles $(M_FILES) --gfiles $(G_FILES) \
	--mdonors $(M_DONORS) --gdonors $(G_DONORS) -o $@

qc: data/qc-mergeruns.Rds
data/qc-mergeruns.Rds: src/performQC.R #data/raw-mergeruns.Rds
	Rscript $< -i data/raw/raw-mergeruns.rds -o $@

markers: data/markers-m.txt data/markers-g.txt data/markers.rds
data/markers.rds: src/markers.R #data/qc-mergeruns.Rds
	Rscript $< -i data/qc-mergeruns.Rds -o data/markers #--quick
data/markers-m.txt: data/markers.rds
	@if test -f $@; then :; else \
		rm -f $<; \
		make $<; \
	fi
data/markers-g.txt: data/markers.rds
	@if test -f $@; then :; else \
		rm -f $<; \
		make $<; \
	fi	

markers-split: data/markers-split-1-g.txt data/markers-split-2-g.txt \
	       data/markers-split-1-m.txt data/markers-split-2-m.txt \
	       data/markers-split.rds
data/markers-split.rds: src/markers.R
	Rscript $< -i data/qc-mergeruns.Rds -o data/markers-split --quick --split
data/markers-split-%.txt: data/markers-split.rds
	@if test -f $@; then :; else \
		rm -f $<; \
		make $<; \
	fi

fig-new-markers: fig/markers-g-stim.svg fig/markers-g-donor.svg fig/markers-m-stim.svg fig/markers-m-donor.svg
fig-can-markers: fig/markers-canonical-stim.svg fig/markers-canonical-donor.svg
fig-markers: fig-new-markers fig-can-markers
fig-qc: fig/qc-metrics-stim.svg fig/qc-metrics-donor.svg
fig/markers-%-stim.svg: src/plotMarkers.R src/plotutils.R #data/qc-mergeruns.Rds
	Rscript $< -i data/qc-mergeruns.Rds -m data/markers-$*.txt -o fig/markers-$*
fig/markers-%-donor.svg: 
	@if test -f $@; then :; else \
		rm -f $<; \
		make $<; \
	fi
fig/qc-metrics-stim.svg: src/plotQC.R src/plotutils.R
	Rscript $< -i data/qc-mergeruns.Rds -m data/qc-metrics.txt -o fig/qc-metrics \
	-ncol 3

fig2: fig/fig2a-markers-canonical-subset.svg \
      fig/fig2c-markers-all-subset.svg \
      fig/fig2d-mf-simp.svg

fig/fig2a-%.svg: src/plotMarkers.R src/plotutils.R
	Rscript $< -i data/qc-mergeruns.Rds -m data/$*.txt -o fig/fig2a-$* \
	-ncol 2

fig/fig2c-%.svg: src/plotMarkers.R src/plotutils.R
	Rscript $< -i data/qc-mergeruns.Rds -m data/$*.txt -o fig/fig2c-$* \
	-ncol 4 

ENRICH_FILES = data/enrich/g-mf.csv data/enrich/m-mf.csv \
	       data/enrich/g-bp.csv data/enrich/m-bp.csv \
	       data/enrich/m-kegg.csv data/enrich/g-kegg.csv \
	       data/enrich/m-reactome.csv data/enrich/g-reactome.csv
enrich: $(ENRICH_FILES)

ENRICH_ALL = data/enrich/all-g-mf.csv data/enrich/all-m-mf.csv \
	     data/enrich/all-g-bp.csv data/enrich/all-m-bp.csv \
	     data/enrich/all-g-kegg.csv data/enrich/all-m-kegg.csv \
	     data/enrich/all-g-reactome.csv data/enrich/all-g-reactome.csv
enrichall: $(ENRICH_ALL)
data/enrich/%-mf.csv: src/pathwayAnalysis.R 
	Rscript $< -i data/background-genes.txt -m data/markers-$*.txt -o data/enrich/$*
data/enrich/%.csv: data/enrich/m-mf.csv
	@if test -f $@; then :; else \
		rm -f $<; \
		make $<; \
	fi

heatmap: fig/heatmap-markers-all.svg fig/heatmap-markers.svg fig/heatmap-markers-canonical.svg
data/markers-canonical.txt: 
data/markers-canonical-subset.txt:
data/markers-%.txt: data/markers-%-g.txt data/markers-%-m.txt
	cat data/$*-g.txt data/$*-m.txt > data/$*.txt ;
fig/heatmap-%.svg: src/plotHeatmap.R data/%.txt src/plotutils.R
	Rscript $< -i data/qc-mergeruns.Rds -m data/$*.txt -o $@

volcano: fig/volcano-markers.svg
fig/volcano-markers.svg: src/plotVolcano.R src/plotutils.R
	Rscript $< -i data/markers.rds -o $@


goplot: fig/g-mf-simp.svg fig/m-mf-simp.svg \
	fig/g-bp-simp.svg fig/m-bp-simp.svg \
	fig/all-g-mf-simp.svg fig/all-m-mf-simp.svg \
	fig/all-g-bp-simp.svg fig/all-m-bp-simp.svg 

fig/fig2d-%-simp.svg: src/plotGO.R src/plotutils.R
	Rscript $< -i data/enrich/g-$*-simp.csv data/enrich/m-$*-simp.csv \
	-o $@
fig/%-simp.svg: src/plotGO.R src/plotutils.R
	Rscript $< -i data/enrich/$*-simp.csv -o $@

