## user defined  parameters are located here
################################################################################
include PARAMETERS


## general behaviour
################################################################################
SHELL:=/bin/bash
.DELETE_ON_ERROR:
ifeq ($(SECONDARY),YES)
# don't delete intermediate files
.SECONDARY:
endif


## paths
################################################################################
# root for fastas and thresholds
PROJDIR:=/home/maticzkd/projects/RBPaffinity
FA_DIR:=$(PROJDIR)/data/fasta
THR_DIR:=$(PROJDIR)/data/thresholds/
# expect binaries to reside in pwd/bin, otherwise this variable must be overwritten
PWD:=$(shell pwd)
BINDIR:=$(PWD)/bin
DATADIR:=$(PWD)/data


## global binaries
################################################################################
PERL:=/usr/local/perl/bin/perl
RBIN:=/usr/local/R/2.15.1-lx/bin/R --vanilla
SVRTRAIN:=/home/maticzkd/src/libsvm-3.12/svm-train -s 3 -t 0 -m $(SVR_CACHE)
SVRPREDICT:=/home/maticzkd/src/libsvm-3.12/svm-predict
PERF:=/home/maticzkd/src/stat/perf
SHUF:=/home/maticzkd/src/coreutils-8.15/src/shuf
#FASTAPL:=$(PERL) /usr/local/user/RNAtools/fastapl
FASTAPL:=/home/maticzkd/co/RNAtools/fastapl
#FASTA2GSPAN:=$(PERL) /usr/local/user/RNAtools/fasta2shrep_gspan.pl
FASTA2GSPAN:=/home/maticzkd/co/RNAtools/fasta2shrep_gspan.pl
# set parralel environment for omp: -v OMP_NUM_THREADS=24 -pe '*' 24
SVMSGDNSPDK:=/home/maticzkd/local/svmsgdnspdk_130606/EDeN
CAT_TABLES:=$(PERL) /home/maticzkd/co/MiscScripts/catTables.pl
BEDGRAPH2BIGWIG:=/usr/local/ucsctools/2012-02/bin/bedGraphToBigWig
BASH:=/bin/bash
BEDTOOLS:=/usr/local/user/BEDTools-Version-2.17.0/bin/bedtools


## project internal tools
################################################################################
LINESEARCH:=$(PERL) $(BINDIR)/lineSearch.pl
COMBINEFEATURES:=$(PERL) $(BINDIR)/combineFeatures.pl
CREATE_EXTENDED_ACC_GRAPH:=$(PERL) $(BINDIR)/createExtendedGraph.pl
MERGE_GSPAN:=$(PERL) $(BINDIR)/merge_gspan.pl
FILTER_FEATURES:=$(PERL) $(BINDIR)/filter_features.pl
SUMMARIZE_MARGINS:=$(PERL) $(BINDIR)/summarize_margins.pl
MARGINS2BG:=$(PERL) $(BINDIR)/margins2bg.pl
VERTEX2NTMARGINS:=$(PERL) $(BINDIR)/vertex2ntmargins.pl
PLOTLC:=$(BASH) $(BINDIR)/plotlc.sh
CHECK_SYNC_GSPAN_CLASS:=$(BASH) $(BINDIR)/check_sync_gspan_class.sh


## set appropriate id (used to determine which parameter sets to use)
################################################################################
ifeq ($(SVM),SVR)
METHOD_ID=svr
endif
ifeq ($(SVM),TOPSVR)
METHOD_ID=svr
endif
ifeq ($(SVM),SGD)
METHOD_ID=sgd
endif


## set targets for RNAcompete evaluation
################################################################################
ifeq ($(EVAL_TYPE),RNACOMPETE)
# filenames for full data sets
FULL_BASENAMES:=$(patsubst %,%_data_full_A,$(PROTEINS)) \
			$(patsubst %,%_data_full_B,$(PROTEINS))

# filenames of data sets containing only weakly structured sequences
BRUIJN_BASENAMES:=$(patsubst %,%_data_bruijn_A,$(PROTEINS)) \
			$(patsubst %,%_data_bruijn_B,$(PROTEINS))

# extract prefixes for further assembling of target filenames
ifeq ($(TRAINING_SETS),FULL)
BASENAMES:=$(FULL_BASENAMES)
else
ifeq ($(TRAINING_SETS),WEAK)
BASENAMES:=$(BRUIJN_BASENAMES)
else
BASENAMES:=$(FULL_BASENAMES) $(BRUIJN_BASENAMES)
endif
endif

# general class statistics
CSTAT_FILES:=$(patsubst %,%.cstats,$(FULL_BASENAMES))

# generate staticstics on positive/negative composition
classstats : summary.cstats $(CSTAT_FILES)

endif


## set targets for CLIP-seq evaluation
################################################################################
ifeq ($(EVAL_TYPE),CLIP)
BASENAMES=$(PROTEINS)
endif


## set targets common to all evaluations
################################################################################
# parameter files (from linesearch or default values)
PARAM_FILES:=$(patsubst %,%.param,$(BASENAMES))
# results of crossvalidations
CV_FILES:=$(patsubst %,%.train.cv,$(BASENAMES))
# models
MODEL_FILES:=$(patsubst %,%.train.model,$(BASENAMES))
# final results spearman correlation
CORRELATION_FILES:=$(patsubst %,%.test.correlation,$(BASENAMES))
# final results from perf
PERF_FILES:=$(patsubst %,%.test.perf,$(BASENAMES))
# nucleotide-wise margins
TESTPART_FILES:=$(patsubst %,%.test.nt_margins.summarized,$(BASENAMES))
# nucleotide-wise margins as bigWig
TESTPART_BIGWIG:=$(patsubst %,%.test.nt_margins.summarized.bw,$(BASENAMES))
# files for learningcurve
LC_FILES:=$(patsubst %,%.lc.png,$(BASENAMES))

## general feature and affinity creation (overridden where apropriate)
################################################################################
%.feature : RADIUS=$(shell grep '^R ' $*.param | cut -f 2 -d' ')
%.feature : DISTANCE=$(shell grep '^D ' $*.param | cut -f 2 -d' ')
%.feature : BITSIZE=$(shell grep '^b ' $*.param | cut -f 2 -d' ')
%.feature : DIRECTED=$(shell grep '^DIRECTED ' $*.param | cut -f 2 -d' ')
%.feature : %.gspan.gz %.affy | %.param
	$(CHECK_SYNC_GSPAN_CLASS) $*.gspan.gz $*.affy
	$(SVMSGDNSPDK) -a FEATURE -i $< -r $(RADIUS) -d $(DISTANCE) -g $(DIRECTED)
	cat $<.feature | grep -v \"^\$\" | paste -d' ' $*.affy - > $@
	-rm -f $<.feature

# extract affinities from fasta
# expected to reside in last field of fasta header
%.affy : %.fa
	$(FASTAPL) -e '@ry = split(/\s/,$$head); print $$ry[-1], "\n"' < $< > $@
#	$(FASTAPL) -e 'print $$head[-1], "\n"' < $< > $@

## receipes specific to graph type
################################################################################
ifeq ($(GRAPH_TYPE),ONLYSEQ)
# line search parameters
LSPAR:=$(DATADIR)/ls.$(METHOD_ID).onlyseq.parameters

%.gspan.gz : VIEWPOINT=$(subst nil,,$(shell grep '^VIEWPOINT ' $*.param | cut -f 2 -d' '))
%.gspan.gz : %.fa | %.param
	$(FASTA2GSPAN) $(VIEWPOINT) --seq-graph-t -nostr -stdout -fasta $< | gzip > $@; exit $${PIPESTATUS[0]}
endif

################################################################################
ifeq ($(GRAPH_TYPE),STRUCTACC)
# line search parameters
LSPAR:=$(DATADIR)/ls.$(METHOD_ID).structacc.parameters

%.gspan.gz : VIEWPOINT=$(subst nil,,$(shell grep '^VIEWPOINT ' $*.param | cut -f 2 -d' '))
%.gspan.gz : %.fa | %.param
	$(CREATE_EXTENDED_ACC_GRAPH) $(VIEWPOINT) -fa $< -W $(W_PRIMARY) -L $(L_PRIMARY) | gzip > $@; exit $${PIPESTATUS[0]}
endif

################################################################################
ifeq ($(GRAPH_TYPE),SHREP)
# line search parameters
LSPAR:=$(DATADIR)/ls.$(METHOD_ID).shrep.parameters

%.gspan.gz : ABSTRACTION=$(shell grep '^ABSTRACTION ' $*.param | cut -f 2 -d' ')
%.gspan.gz : STACK=$(subst nil,,$(shell grep '^STACK ' $*.param | cut -f 2 -d' '))
%.gspan.gz : CUE=$(subst nil,,$(shell grep '^CUE ' $*.param | cut -f 2 -d' '))
%.gspan.gz : VIEWPOINT=$(subst nil,,$(shell grep '^VIEWPOINT ' $*.param | cut -f 2 -d' '))
%.gspan.gz : %.fa | %.param
	$(FASTA2GSPAN) --seq-graph-t --seq-graph-alph $(STACK) $(CUE) $(VIEWPOINT) -stdout -t $(ABSTRACTION) -M 3 -wins '$(SHAPES_WINS)' -shift '$(SHAPES_SHIFT)' -fasta $< | gzip > $@; exit $${PIPESTATUS[0]}
endif

################################################################################
ifeq ($(GRAPH_TYPE),PROBSHREP)
# line search parameters
LSPAR:=$(DATADIR)/ls.$(METHOD_ID).shrep.parameters

%.gspan.gz : ABSTRACTION=$(shell grep '^ABSTRACTION ' $*.param | cut -f 2 -d' ')
%.gspan.gz : STACK=$(subst nil,,$(shell grep '^STACK ' $*.param | cut -f 2 -d' '))
%.gspan.gz : CUE=$(subst nil,,$(shell grep '^CUE ' $*.param | cut -f 2 -d' '))
%.gspan.gz : VIEWPOINT=$(subst nil,,$(shell grep '^VIEWPOINT ' $*.param | cut -f 2 -d' '))
%.gspan.gz : %.fa | %.param
	$(FASTA2GSPAN) $(STACK) $(CUE) $(VIEWPOINT) -stdout -q -Tp 0.05 -t $(ABSTRACTION) -M 3 -wins '$(SHAPES_WINS)' -shift '$(SHAPES_SHIFT)' -fasta $< | gzip > $@; exit $${PIPESTATUS[0]}

%.feature : RADIUS=$(shell grep '^R ' $*.param | cut -f 2 -d' ')
%.feature : DISTANCE=$(shell grep '^D ' $*.param | cut -f 2 -d' ')
%.feature : BITSIZE=$(shell grep '^b ' $*.param | cut -f 2 -d' ')
%.feature : DIRECTED=$(shell grep '^DIRECTED ' $*.param | cut -f 2 -d' ')
%.feature : %.gspan.gz %.affy | %.param
	# remove t and w, convert s to t
	cat $< | grep -v -e '^t' -e '^w' | sed 's/^s/t/' > $*_singleshreps
	# write out probabilities
	cat $< | grep '^s' | $(PERL) -ne '($$prob) = /SHAPEPROB ([0-9.]+)/; print $$prob, "\n"' > $*_probs
	# write out shrep membership
	cat $< | awk '/^t/{i++}/^s/{print i}' > $*_groups
	# compute features
	$(SVMSGDNSPDK) -a FEATURE -i $*_singleshreps -r $(RADIUS) -d $(DISTANCE) -g $(DIRECTED)
	# compute probability-weighted features
	$(COMBINEFEATURES) $*_singleshreps.feature $*_probs $*_groups > $*
	# add affinities to features
	cat $* | grep -v \"^\$\" | paste -d' ' $*.affy - > $@
	-rm -f $*_singleshreps.feature
endif

################################################################################
ifeq ($(GRAPH_TYPE),CONTEXTSHREP)
# line search parameters
LSPAR:=$(DATADIR)/ls.$(METHOD_ID).shrep_context.parameters

%.gspan.gz : ABSTRACTION=$(shell grep '^ABSTRACTION ' $*.param | cut -f 2 -d' ')
%.gspan.gz : STACK=$(subst nil,,$(shell grep '^STACK ' $*.param | cut -f 2 -d' '))
%.gspan.gz : CUE=$(subst nil,,$(shell grep '^CUE ' $*.param | cut -f 2 -d' '))
%.gspan.gz : VIEWPOINT=$(subst nil,,$(shell grep '^VIEWPOINT ' $*.param | cut -f 2 -d' '))
%.gspan.gz : %.fa | %.param
	$(FASTA2GSPAN) $(STACK) $(CUE) $(VIEWPOINT) --seq-graph-t --seq-graph-alph -abstr -stdout -t $(ABSTRACTION) -M 3 -wins '$(SHAPES_WINS)' -shift '$(SHAPES_SHIFT)' -fasta $< | gzip > $@; exit $${PIPESTATUS[0]}

%.sequence : %.gspan.gz
	zcat $< | awk '/^t/{print $$NF}' > $@

%.top_wins : %.nt_margins.summarized selectTopWinShreps.R
	/usr/local/R/2.15.1-lx/bin/R --slave --no-save --args $< < selectTopWinShreps.R | sort -k3,3nr | head -n $(TOP_WINDOWS) | sort -k1,1n > $@

%.sequence_top_wins : %.sequence %.top_wins
	$(PERL) subTopWins.pl --input $< --locations $*.top_wins --win_size $(MARGINS_WINDOW) > $@

%.struct_annot_top_wins : %.struct_annot %.top_wins
	$(PERL) subTopWins.pl --input $< --locations $*.top_wins --win_size $(MARGINS_WINDOW) > $@

%.truncated : %
	cat $< | awk 'length($$0)==$(MARGINS_WINDOW)' > $@

%.pup : %
	cat $< | tr 'HBIEM' 'UUUUU' | tr 'S' 'P' > $@

# %.gspan.gz : ABSTRACTION=$(shell grep '^ABSTRACTION ' $*.param | cut -f 2 -d' ')
# %.gspan.gz : STACK=$(subst nil,,$(shell grep '^STACK ' $*.param | cut -f 2 -d' '))
# %.gspan.gz : CUE=$(subst nil,,$(shell grep '^CUE ' $*.param | cut -f 2 -d' '))
# %.gspan.gz : VIEWPOINT=$(subst nil,,$(shell grep '^VIEWPOINT ' $*.param | cut -f 2 -d' '))
# %.gspan.gz : %.fa | %.param
# 	$(FASTA2GSPAN) $(STACK) $(CUE) $(VIEWPOINT) --seq-graph-t --seq-graph-alph -abstr -stdout -t $(ABSTRACTION) -M 3 -wins '$(SHAPES_WINS)' -shift '$(SHAPES_SHIFT)' -fasta $< | gzip > $@; exit $${PIPESTATUS[0]}

endif

################################################################################
ifeq ($(GRAPH_TYPE),MEGA)
# line search parameters
LSPAR:=$(DATADIR)/ls.$(METHOD_ID).mega.parameters

# accessibility graphs
%.acc.gspan.gz : VIEWPOINT=$(subst nil,,$(shell grep '^VIEWPOINT ' $*.param | cut -f 2 -d' '))
%.acc.gspan.gz : %.fa
	$(CREATE_EXTENDED_ACC_GRAPH) $(VIEWPOINT) -fa $< -W $(W_PRIMARY) -L $(L_PRIMARY) | gzip > $@; exit $${PIPESTATUS[0]}

# shrep graphs
%.shrep.gspan.gz : ABSTRACTION=$(shell grep '^ABSTRACTION ' $*.param | cut -f 2 -d' ')
%.shrep.gspan.gz : STACK=$(subst nil,,$(shell grep '^STACK ' $*.param | cut -f 2 -d' '))
%.shrep.gspan.gz : CUE=$(subst nil,,$(shell grep '^CUE ' $*.param | cut -f 2 -d' '))
%.shrep.gspan.gz : VIEWPOINT=$(subst nil,,$(shell grep '^VIEWPOINT ' $*.param | cut -f 2 -d' '))
%.shrep.gspan.gz : %.fa | %.param
	$(FASTA2GSPAN) --seq-graph-t --seq-graph-alph $(STACK) $(CUE) $(VIEWPOINT) -stdout -t $(ABSTRACTION) -M 3 -wins '$(SHAPES_WINS)' -shift '$(SHAPES_SHIFT)' -fasta $< | gzip > $@; exit $${PIPESTATUS[0]}

# merge gspans
%.gspan.gz : %.shrep.gspan %.acc.gspan
	$(MERGE_GSPAN) -shrep $*.shrep.gspan -acc $*.acc.gspan | gzip > $@; exit $${PIPESTATUS[0]}
endif


## receipes specific to SVM type
################################################################################
# support vector regression
ifeq ($(SVM),SVR)
# results from cross validation
%.cv_svr : C=$(shell grep '^c ' $*.param | cut -f 2 -d' ')
%.cv_svr : EPSILON=$(shell grep '^e ' $*.param | cut -f 2 -d' ')
%.cv_svr : %.feature | %.param
	time $(SVRTRAIN) -c $(C) -p $(EPSILON) -h 0 -v $(CV_FOLD) $< > $@

# final result of cross validation: squared correlation coefficient
%.cv : %.cv_svr
	cat $< | grep 'Cross Validation Squared correlation coefficient' | perl -ne 'print /(\d+.\d+)/' > $@

# SVR model
%.model : C=$(shell grep '^c' $*.param | cut -f 2 -d' ')
%.model : EPSILON=$(shell grep '^e' $*.param | cut -f 2 -d' ')
%.model : %.feature | %.param
	time $(SVRTRAIN) -c $(C) -p $(EPSILON) $< $@

# SVR predictions
%.test.predictions_svr : %.train.model %.test.feature
	time $(SVRPREDICT) $*.test.feature $< $@

# affinities and predictions default format
%.predictions_affy : %.predictions_svr %.affy
	paste $*.affy $< > $@

# class membership and predictions default format
%.predictions_class : %.predictions_svr %.class
	paste $*.class $< > $@

endif


## support vector regression using sgd-derived subset of top features
################################################################################
ifeq ($(SVM),TOPSVR)
# results from cross validation
%.cv_svr : C=$(shell grep '^c ' $*.param | cut -f 2 -d' ')
%.cv_svr : EPSILON=$(shell grep '^e ' $*.param | cut -f 2 -d' ')
%.cv_svr : %.feature | %.param
	time $(SVRTRAIN) -c $(C) -p $(EPSILON) -h 0 -v $(CV_FOLD) $< > $@

# final result of cross validation: squared correlation coefficient
%.cv : %.cv_svr
	cat $< | grep 'Cross Validation Squared correlation coefficient' | perl -ne 'print /(\d+.\d+)/' > $@

# train model; this one directly works on gspans
%.sgd_model : EPOCHS=$(shell grep '^EPOCHS ' $*.param | cut -f 2 -d' ')
%.sgd_model : LAMBDA=$(shell grep '^LAMBDA ' $*.param | cut -f 2 -d' ')
%.sgd_model : RADIUS=$(shell grep '^R ' $*.param | cut -f 2 -d' ')
%.sgd_model : DISTANCE=$(shell grep '^D ' $*.param | cut -f 2 -d' ')
%.sgd_model : BITSIZE=$(shell grep '^b ' $*.param | cut -f 2 -d' ')
%.sgd_model : DIRECTED=$(shell grep '^DIRECTED ' $*.param | cut -f 2 -d' ')
%.sgd_model : %.gspan.gz %.class | %.param
	$(CHECK_SYNC_GSPAN_CLASS) $*.gspan.gz $*.class
	$(SVMSGDNSPDK) -g $(DIRECTED) -b $(BITSIZE) -a TRAIN -i $*.gspan.gz -t $*.class -m $@ -r $(RADIUS) -d $(DISTANCE) -e $(EPOCHS) -l $(LAMBDA)

%.test.filter : %.train.filter
	ln -s $< $@

%.filter : NFEAT=$(shell cat $< | grep '^w ' | sed 's/^w //' | tr ' :' "\n\t" | wc -l)
%.filter : TENP=$(shell echo "$(NFEAT) / 5" | bc)
%.filter : %.sgd_model
	cat $< | grep '^w ' | sed 's/^w //' | tr ' :' "\n\t" | sort -k2,2gr | head -n $(TENP) | cut -f 1 | sort -n > $@

%.feature_filtered : %.feature %.filter
	$(FILTER_FEATURES) --feature $< --filter $*.filter > $@

# SVR model
%.model : C=$(shell grep '^c' $*.param | cut -f 2 -d' ')
%.model : EPSILON=$(shell grep '^e' $*.param | cut -f 2 -d' ')
%.model : %.feature_filtered | %.param
	time $(SVRTRAIN) -c $(C) -p $(EPSILON) $< $@

# SVR predictions
%.test.predictions_svr : %.train.model %.test.feature_filtered
	time $(SVRPREDICT) $*.test.feature_filtered $< $@

# affinities and predictions default format
%.predictions_affy : %.predictions_svr %.affy
	# combine affinities and predictions
	paste $*.affy $< > $@

# class membership and predictions default format
%.predictions_class : %.predictions_svr %.class
	# combine affinities and predictions
	paste $*.class $< > $@
endif


## stochastic gradient descent
################################################################################
ifeq ($(SVM),SGD)
# extract single performance measure, used for linesearch decisions
%.cv : %.cv.perf
	cat $< | grep 'APR' | awk '{print $$NF}' > $@

# train model; this one directly works on gspans
%.model : EPOCHS=$(shell grep '^EPOCHS ' $*.param | cut -f 2 -d' ')
%.model : LAMBDA=$(shell grep '^LAMBDA ' $*.param | cut -f 2 -d' ')
%.model : RADIUS=$(shell grep '^R ' $*.param | cut -f 2 -d' ')
%.model : DISTANCE=$(shell grep '^D ' $*.param | cut -f 2 -d' ')
%.model : BITSIZE=$(shell grep '^b ' $*.param | cut -f 2 -d' ')
%.model : DIRECTED=$(shell grep '^DIRECTED ' $*.param | cut -f 2 -d' ')
%.model : %.gspan.gz %.class | %.param
	$(CHECK_SYNC_GSPAN_CLASS) $*.gspan.gz $*.class
	$(SVMSGDNSPDK) -g $(DIRECTED) -e $(EPOCHS) -l $(LAMBDA) -b $(BITSIZE) -a TRAIN -i $*.gspan.gz -t $*.class -m $@ -r $(RADIUS) -d $(DISTANCE)

# evaluate model
%.test.predictions_sgd : EPOCHS=$(shell grep '^EPOCHS ' $*.param | cut -f 2 -d' ')
%.test.predictions_sgd : LAMBDA=$(shell grep '^LAMBDA ' $*.param | cut -f 2 -d' ')
%.test.predictions_sgd : RADIUS=$(shell grep '^R ' $*.param | cut -f 2 -d' ')
%.test.predictions_sgd : DISTANCE=$(shell grep '^D ' $*.param | cut -f 2 -d' ')
%.test.predictions_sgd : BITSIZE=$(shell grep '^b ' $*.param | cut -f 2 -d' ')
%.test.predictions_sgd : DIRECTED=$(shell grep '^DIRECTED ' $*.param | cut -f 2 -d' ')
%.test.predictions_sgd : %.train.model %.test.gspan.gz %.test.class | %.param
	$(SVMSGDNSPDK) -g $(DIRECTED) -r $(RADIUS) -d $(DISTANCE) -e $(EPOCHS) -l $(LAMBDA) -b $(BITSIZE) -a TEST -m $< -i $*.test.gspan.gz -t $*.test.class
	mv $*.test.gspan.gz.prediction $*.test.predictions_sgd

# affinities and predictions default format
%.predictions_affy : %.predictions_sgd %.affy
	cat $< | awk '{print $$2}' | paste $*.affy - > $@

# class membership and predictions default format: class{-1,1}, prediction
%.predictions_class : %.predictions_sgd %.class
	cat $< | awk '{print $$2}' | paste $*.class - > $@

# results from crossvalidation cast into default format: class{-1,1}, prediction
%.cv.predictions_class : EPOCHS=$(shell grep '^EPOCHS ' $*.param | cut -f 2 -d' ')
%.cv.predictions_class : LAMBDA=$(shell grep '^LAMBDA ' $*.param | cut -f 2 -d' ')
%.cv.predictions_class : RADIUS=$(shell grep '^R ' $*.param | cut -f 2 -d' ')
%.cv.predictions_class : DISTANCE=$(shell grep '^D ' $*.param | cut -f 2 -d' ')
%.cv.predictions_class : BITSIZE=$(shell grep '^b ' $*.param | cut -f 2 -d' ')
%.cv.predictions_class : DIRECTED=$(shell grep '^DIRECTED ' $*.param | cut -f 2 -d' ')
%.cv.predictions_class : %.gspan.gz %.class | %.param
	$(CHECK_SYNC_GSPAN_CLASS) $*.gspan.gz $*.class
	$(SVMSGDNSPDK) -g $(DIRECTED) -b $(BITSIZE) -a CROSS_VALIDATION -c $(CV_FOLD) -m $*.model -i $< -t $*.class -r $(RADIUS) -d $(DISTANCE) -e $(EPOCHS) -l $(LAMBDA)
	cat $<.cv_predictions | awk '{print $$2==1?1:-1, $$4}' > $@
	-rm  -f $<.cv_predictions$* $*.model_*

# compute margins of graph vertices
# vertex_margins format: seqid verticeid margin
%.test.vertex_margins : EPOCHS=$(shell grep '^EPOCHS ' $*.param | cut -f 2 -d' ')
%.test.vertex_margins : LAMBDA=$(shell grep '^LAMBDA ' $*.param | cut -f 2 -d' ')
%.test.vertex_margins : RADIUS=$(shell grep '^R ' $*.param | cut -f 2 -d' ')
%.test.vertex_margins : DISTANCE=$(shell grep '^D ' $*.param | cut -f 2 -d' ')
%.test.vertex_margins : BITSIZE=$(shell grep '^b ' $*.param | cut -f 2 -d' ')
%.test.vertex_margins : DIRECTED=$(shell grep '^DIRECTED ' $*.param | cut -f 2 -d' ')
%.test.vertex_margins : %.test.gspan.gz %.test.class %.train.model | %.param
	$(SVMSGDNSPDK) -g $(DIRECTED) -r $(RADIUS) -d $(DISTANCE) -b $(BITSIZE) -e $(EPOCHS) -l $(LAMBDA) -a TEST_PART -m $*.train.model -i $*.test.gspan.gz -t $*.test.class
	mv $<.prediction_part $*.test.vertex_margins

# dictionary of all graph vertices
# dict file format: seqid verticeid nt pos
%.vertex_dict : %.gspan.gz
	zcat $< | awk 'BEGIN{seqid=-1}/^t/{seqid++; vertex_id=0; nt_pos=0}/^v/&&!/\^/&&!/#/{print seqid, vertex_id++, $$3, nt_pos++}/^V/&&!/\^/&&!/#/{nt_pos++}' > $@

# compute nucleotide-wise margins from vertice margins
%.nt_margins : %.vertex_margins %.vertex_dict
	cat $< | $(VERTEX2NTMARGINS) -dict $*.vertex_dict | awk '$$2!=0' > $@

# format (tab separated): sequence id, sequence position, margin,
#                         min, max, mean, median, sum
%.nt_margins.summarized : %.nt_margins
	@echo ""
	@echo "summarizing nucleotide-wise margins:"
	$(SUMMARIZE_MARGINS) -W $(MARGINS_WINDOW) < $< > $@

%.nt_margins.summarized.bedGraph : %.nt_margins.summarized %.bed
	@echo ""
	@echo "converting margins to bedGraph"
	$(MARGINS2BG) -bed $*.bed --aggregate $(MARGINS_MEASURE) < $< | \
	$(BEDTOOLS) sort > $@

# compute learningcurve
# svmsgdnspdk creates LEARNINGCURVE_SPLITS many files of the format output.lc_predictions_{test,train}_fold{1..LEARNINGCURVE_SPLITS.ID}
# format: id, class, prediction, margin
# we evaluate each one and summarize the probabilities in the following format:
# SPLIT train_performance test_performance
# this is done using svmsgdnspdk default parameters
%.lc.perf : EPOCHS=$(shell grep '^EPOCHS ' $*.train.param | cut -f 2 -d' ')
%.lc.perf : LAMBDA=$(shell grep '^LAMBDA ' $*.train.param | cut -f 2 -d' ')
%.lc.perf : RADIUS=$(shell grep '^R ' $*.train.param | cut -f 2 -d' ')
%.lc.perf : DISTANCE=$(shell grep '^D ' $*.train.param | cut -f 2 -d' ')
%.lc.perf : BITSIZE=$(shell grep '^b ' $*.train.param | cut -f 2 -d' ')
%.lc.perf : DIRECTED=$(shell grep '^DIRECTED ' $*.train.param | cut -f 2 -d' ')
%.lc.perf : %.train.class %.train.gspan.gz
	-rm -f $*.train.dat_lc;
	LC=10; \
	NUM_REP=10; \
	lcn=$$((LC+1)); \
	for r in $$(seq 1 $$NUM_REP); \
	do $(SVMSGDNSPDK) -a LEARNING_CURVE -g $(DIRECTED) -r $(RADIUS) -d $(DISTANCE) -b $(BITSIZE) -e $(EPOCHS) -l $(LAMBDA) -i $*.train.gspan.gz -t $*.train.class -p $$lcn -R $$r > /dev/null; \
	for i in $$(seq 1 $$LC); \
	do \
	dim=$$(cat  $*.train.gspan.gz.lc_predictions_train_fold_$$i | wc -l); \
	echo -n "calculating learningcurve iteration $$i"; \
	echo -n "$$dim " >> $*.train.dat_lc; \
	cat $*.train.gspan.gz.lc_predictions_train_fold_$$i | \
	awk '{print $$2,$$4}' | $(PERF) -APR -ROC -ACC -t 0 -PRF 2> /dev/null | \
	awk '{printf("%s %s ",$$1,$$2)}END{printf("\n")}' >> $*.train.dat_lc; \
	cat $*.train.gspan.gz.lc_predictions_test_fold_$$i | awk '{print $$2,$$4}' | \
	$(PERF) -APR -ROC -ACC -t 0 -PRF 2> /dev/null | \
	awk '{printf("%s %s ",$$1,$$2)}END{printf("\n")}' >> $*.train.dat_lc; \
	done; \
	done; \
	cat $*.train.dat_lc | awk 'NR%2==1{printf("%s ",$$0)}NR%2==0{print $$0}' | column -t > $@
	-rm -f $*.train.gspan.gz.lc_predictions_t*_fold_* $*.train.dat_lc

%.lc.png : %.lc.perf
	$(PLOTLC) $< $@
	cat $@.fit_log | \
	grep 'ats' | \
	grep '=' | \
	tail -n 1 | \
	awk '{print $$3}' > $@.train_limit

endif


## evaluations specific to RNAcompete analysis
################################################################################
ifeq ($(EVAL_TYPE),RNACOMPETE)

# class memberships {-1,0,1}
%.class : BASENAME=$(firstword $(subst _, ,$<))
%.class : HT=$(shell grep $(BASENAME) $(THR_DIR)/positive.txt | cut -f 2 -d' ')
%.class : LT=$(shell grep $(BASENAME) $(THR_DIR)/negative.txt | cut -f 2 -d' ')
%.class : %.affy
	cat $< | awk '{ if ($$1 > $(HT)) {print 1} else { if ($$1 < $(LT)) {print -1} else {print 0} } }' > $@

# some statistics about class distribution
%.cstats : BASENAME=$(firstword $(subst _, ,$<))
%.cstats : TYPE=$(word 3,$(subst _, ,$<))
%.cstats : SET=$(word 4,$(subst ., ,$(subst _, ,$<)))
%.cstats : HT=$(shell grep $(BASENAME) $(THR_DIR)/positive.txt | cut -f 2 -d' ')
%.cstats : LT=$(shell grep $(BASENAME) $(THR_DIR)/negative.txt | cut -f 2 -d' ')
%.cstats : HN=$(shell cat $< | grep '^>' | awk '$$NF > $(HT)' | wc -l)
%.cstats : LN=$(shell cat $< | grep '^>' | awk '$$NF < $(LT)' | wc -l)
%.cstats : %.fa
	$(PERL) -e 'print join("\t", "$(BASENAME)", "$(SET)", "$(LT)", "$(LN)", "$(HT)", "$(HN)"),"\n"' > $@

# final class summary
summary.cstats : $(CSTAT_FILES)
	( $(PERL) -e 'print join("\t", "protein", "set", "negative threshold", "negative instances", "positive threshold", "positive instances"),"\n"'; \
	cat $^ | sort -k1,2 ) > $@
endif


## evaluations specific to CLIP analysis
################################################################################
ifeq ($(EVAL_TYPE),CLIP)
# combine input sequences
%.fa : %.positives.fa %.negatives.fa %.unknowns.fa
	( $(FASTAPL) -p -1 -e '$$head .= " 1";' < $*.positives.fa; \
	  $(FASTAPL) -p -1 -e '$$head .= " -1";' < $*.negatives.fa; \
	  $(FASTAPL) -p -1 -e '$$head .= " 0";' < $*.unknowns.fa ) > $@

%.bed : %.positives.bed %.negatives.bed %.unknowns.bed
	cat $^ > $@

# for clip data, affinities are actually the class
%.class : %.affy
	ln -sf $< $@
endif


## misc helper receipes
################################################################################

# if needed, unzip gspans (for svr)
%.gspan : %.gspan.gz
	zcat $< > $@

ifeq ($(DO_LINESEARCH),NO)
# just use defaults instead of doing line search
%.param : $(LSPAR)
	cut -f 1,2 -d' ' < $< > $@
else
ifeq ($(DO_SGDOPT),YES)
# do parameter optimization by line search but also use sgd-internal optimization
%.param : %.ls.fa $(LSPAR)
	$(LINESEARCH) -fa $< -param $(LSPAR) -mf Makefile -of $@ -bindir $(PWD) -sgdopt 2> >(tee $@.log >&2)
# call sgdsvmnspdk optimization and write file containing optimized parameters

%.ls.param : EPOCHS=$(shell grep '^EPOCHS ' $*.ls_sgdopt.param | cut -f 2 -d' ')
%.ls.param : LAMBDA=$(shell grep '^LAMBDA ' $*.ls_sgdopt.param | cut -f 2 -d' ')
%.ls.param : RADIUS=$(shell grep '^R ' $*.ls_sgdopt.param | cut -f 2 -d' ')
%.ls.param : DISTANCE=$(shell grep '^D ' $*.ls_sgdopt.param | cut -f 2 -d' ')
%.ls.param : BITSIZE=$(shell grep '^b ' $*.ls_sgdopt.param | cut -f 2 -d' ')
%.ls.param : DIRECTED=$(shell grep '^DIRECTED ' $*.ls_sgdopt.param | cut -f 2 -d' ')
%.ls.param : %.ls_sgdopt.param %.ls_sgdopt.gspan.gz %.ls_sgdopt.class
	$(SVMSGDNSPDK) -g $(DIRECTED) -r $(RADIUS) -d $(DISTANCE) -e $(EPOCHS) \
	-l $(LAMBDA) -b $(BITSIZE) -a PARAMETERS_OPTIMIZATION \
	-i $*.ls_sgdopt.gspan.gz -t $*.ls_sgdopt.class -m $@ -p $(SGDOPT_STEPS) > /dev/null
	( cat $< | grep -v -e '^D ' -e '^R ' -e '^EPOCHS ' -e '^LAMBDA '; \
	cat $*.ls_sgdopt.gspan.gz.opt_param | awk '{print "R",$$2,"\nD",$$4,"\nEPOCHS",$$6,"\nLAMBDA",$$8}' \
	) > $@
	rm $*.ls_sgdopt.gspan.gz.opt_param
else
# do parameter optimization by line search
%.param : %.ls.fa $(LSPAR)
	$(LINESEARCH) -fa $< -param $(LSPAR) -mf Makefile -of $@ -bindir $(PWD) 2> >(tee $@.log >&2)
endif
endif

# subset fastas prior to line search
%.ls.fa : %.train.fa
	cat $< | \
	$(FASTAPL) -e 'print ">", $$head, "\t", $$seq, "\n"' | \
	$(SHUF) -n $(LINESEARCH_INPUT_SIZE) | \
	$(PERL) -ane \
	'$$seq = pop @F; $$head = join(" ", @F); print $$head, "\n", $$seq, "\n";' > \
	$@

%.ls_sgdopt.fa : %.ls.fa
	ln -s $< $@

%.ls_sgdopt.affy : %.ls.affy
	ln -s $< $@

# link parameter files
%.train.param : %.param
	ln -sf $< $@

%.test.param : %.param
	ln -sf $< $@

# get test data
testclip.train.positives.fa : $(DATADIR)/testclip.train.positives.fa
	cp -f $< $@

testclip.train.negatives.fa : $(DATADIR)/testclip.train.negatives.fa
	cp -f $< $@

testclip.test.positives.fa : $(DATADIR)/testclip.test.positives.fa
	cp -f $< $@

testclip.test.negatives.fa : $(DATADIR)/testclip.test.negatives.fa
	cp -f $< $@

test_data_full_A.test.fa : $(DATADIR)/test_data_full_A.test.fa
	cp -f $< $@

test_data_full_A.train.fa : $(DATADIR)/test_data_full_A.train.fa
	cp -f $< $@

# create empty input
%.unknowns.fa :
	@echo ""
	@echo "using empty set of unknowns!"
	touch $@

# %.negatives.fa :
# 	@echo ""
# 	@echo "using empty set of negatives!"
# 	touch $@

%.unknowns.bed :
	@echo ""
	@echo "using empty set of unknowns!"
	touch $@

# %.negatives.bed :
# 	@echo ""
# 	@echo "using empty set of negatives!"
# 	touch $@

# compute performance measures
# remove unknowns, set negative class to 0 for perf
%.perf : %.predictions_class
	cat $< | awk '$$1!=0' | sed 's/^-1/0/g' | $(PERF) -confusion > $@

# plot precision-recall
%.prplot : %.predictions_class
	cat $< | sed 's/^-1/0/g' | $(PERF) -plot pr | awk 'BEGIN{p=1}/ACC/{p=0}{if (p) {print}}' > $@

%.prplot.svg : %.prplot
	cat $< | gnuplot -e "set ylabel 'precision'; set xlabel 'recall'; set terminal svg; set style line 1 linecolor rgb 'black'; plot [0:1] [0:1] '-' using 1:2 with lines;" > $@

# compute correlation: correlation \t pvalue
%.correlation : %.predictions_affy
	cat $< | $(RBIN) --slave -e 'require(stats); data=read.table("$<", col.names=c("prediction","measurement")); t <- cor.test(data$$measurement, data$$prediction, method="spearman", alternative="greater"); write.table(cbind(t$$estimate, t$$p.value), file="$@", col.names=F, row.names=F, quote=F, sep="\t")'

results_aucpr.csv : $(PERF_FILES)
	grep -H -e APR -e ROC $^ | tr ':' "\t" | $(RBIN) --slave -e 'require(reshape); d<-read.table("stdin", col.names=c("id","variable","value")); write.table( cast(d), file="", row.names=F, quote=F, sep="\t")' > $@

results_correlation.csv : $(CORRELATION_FILES)
	$(CAT_TABLES) $(CORRELATION_FILES) > $@

# convert bedGraph to bigWig
%.bw : %.bedGraph $(GENOME_BOUNDS)
	$(BEDGRAPH2BIGWIG) $*.bedGraph $(GENOME_BOUNDS) $@

# # do need genome bounds
# $(GENOME_BOUNDS) :
# 	@echo ""
# 	@echo "error: require genome boundaries $@" && exit 255

## phony target section
################################################################################
.PHONY: all ls cv classstats test clean distclean learningcurve

# do predictions and tests for all PROTEINS, summarize results
all: $(PERF_FILES) $(CORRELATION_FILES) results_aucpr.csv results_correlation.csv

# do parameter line search for all PROTEINS
ls : $(PARAM_FILES)

# do crossvalidation
cv : $(CV_FILES)

# train target
train : $(MODEL_FILES)

# test target
test : $(PERF_FILES) $(CORRELATION_FILES)

# compute nucleotide-wise margins
testpart : $(TESTPART_FILES)

# compute nucleotide-wise margins for genome-browser
testpart_coords : $(TESTPART_BIGWIG)

# see if additional data will help improve classification
learningcurve: $(LC_FILES)

# keep fasta, predictions and results
clean:
	-rm -rf log *.gspan *.gspan.gz *.threshold* *.feature *.affy *.feature_filtered \
	*.filter *.class

# delete all files
distclean: clean
	-rm -rf *.param *.perf *.predictions_class *.predictions_affy \
	*.predictions_svr *.predictions_sgd *.ls.fa *.log *.csv *model \
	*.sgeout *.class *.correlation *.cv *.cv.predictions \
	*.cv_svr *.model_* *.prplot *.prplot.svg $(LC_FILES) *.nt_margins* \
	*.vertex_margins *.vertex_dict

ifeq ($(EVAL_TYPE),CLIP)
# test various stuff
runtests: testclip.param testclip.train.cv \
	testclip.test.perf testclip.test.correlation \
	testclip.test.prplot.svg
endif
ifeq ($(EVAL_TYPE),RNACOMPETE)
# test various stuff
runtests: test_data_full_A.param test_data_full_A.train.cv \
	test_data_full_A.test.perf test_data_full_A.test.correlation \
	test_data_full_A.test.prplot.svg
endif

## miscellaneous rules
################################################################################

# load genome sizes from ucsc
%.tab :
	mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e \
	"select chrom, size from $*.chromInfo" | grep -v size > $@

# get sequence from bed using twoBitToFa
TWOBITTOFA:=/usr/local/ucsctools/2012-02/bin/twoBitToFa
%.positives.fa : %.positives.bed
	 $(TWOBITTOFA) -bed=$< $(GENOME) $@

# get sequence from bed using twoBitToFa
TWOBITTOFA:=/usr/local/ucsctools/2012-02/bin/twoBitToFa
%.negatives.fa : %.negatives.bed
	 $(TWOBITTOFA) -bed=$< $(GENOME) $@

# get sequence from bed using twoBitToFa
TWOBITTOFA:=/usr/local/ucsctools/2012-02/bin/twoBitToFa
%.unknowns.fa : %.unknowns.bed
	 $(TWOBITTOFA) -bed=$< $(GENOME) $@

## insert additional rules into this file
################################################################################
include EXPERIMENT_SPECIFIC_RULES


## old code for MEME search
################################################################################
# # binaries
# MEME_GETMARKOV:=/home/maticzkd/src/meme_4.7.0/local/bin/fasta-get-markov
# MEME:=/home/maticzkd/src/meme_4.7.0/local/bin/meme
# FASTAUID:=/usr/local/user/RNAtools/fastaUID.pl
# # perform meme oops (only one per sequence) search
# meme_oops: positives_unique.fa negatives_markov0.txt
# 	$(MEME) $< -mod oops -maxsites 3 -minw 5 -maxw 15 -bfile negatives_markov0.txt -dna -nmotifs 5 -maxsize 300000 -oc $@
#
# # perform meme zoops (zero or one per sequence) search
# meme_zoops: positives_unique.fa negatives_markov0.txt
# 	$(MEME) $< -mod zoops -maxsites 3 -minw 5 -maxw 15 -bfile negatives_markov0.txt -dna -nmotifs 5 -maxsize 300000 -oc $@
#
# # assign unique ids to fasta headers
# %_unique.fa: %.fa
# 	$(FASTAUID) -id pos_ < $< > $@
#
# # create background model for meme motif search
# negatives_markov0.txt: negatives.fa
# 	$(MEME_GETMARKOV) -norc < $< > $@
