############################ LOAD LIBRARIES ############################

# load libraries
library("plyr")
library("dplyr")
library("reshape2")
library("ggplot2")
library("phyloseq")
library("ANCOMBC")
library("breakaway")
library("DivNet")

############################ LOAD DATA ############################

print("LOADING DATA")

# set path according to the experiment
path_to_exp <- "/mnt/cinqueg/gabriele/work/microbiology/other_analyses/taxonomy_example_final_plots/"
# set paths
save_img <- paste0(path_to_exp, "figures/")
save_obj <- paste0(path_to_exp, "rdsObjs/")
save_taxa_figs <- paste0(save_img, "taxaFigs/")

# check if dirs exist
ifelse(!dir.exists(save_img), dir.create(save_img, recursive=TRUE), "dir exists!")
ifelse(!dir.exists(save_obj), dir.create(save_obj, recursive=TRUE), "dir exists!")
ifelse(dir.exists(save_taxa_figs), "dir exists!", dir.create(save_taxa_figs, recursive=TRUE))

# load asv table
asvs <- read.csv(paste0(path_to_exp, "asv_table.csv"), header=T, sep="\t")
metadata <- read.csv(paste0(path_to_exp, "metadata.csv"), header=T, sep="\t")
taxonomy <- read.csv(paste0(path_to_exp, "taxonomy.csv"), header=T, sep="\t")

############################ CREATE PHYLOSEQ OBJECTS ############################

# set rownames for asvs table
rownames(asvs) <- asvs$asv
# remove first column containing asvs names. it is now useless
asvs <- asvs[, -1]

# set rownamames for taxonomy
rownames(taxonomy) <- taxonomy$asv
taxonomy <- taxonomy[, -1]

# set rownames for metadata
rownames(metadata) <- metadata$sampleID

# load non-transformed object
phylo_data <- phyloseq(otu_table(asvs, taxa_are_rows=T), sample_data(metadata), tax_table(as(taxonomy, "matrix")))

# load metadata table
phylo_metadata <- sample_data(phylo_data)

############################ SELECT PARAMETERS FOR PLOTTING ############################

# loop through taxa levels, to aggregate appropriately
taxa_level <- "Order"
# get best n taxa for the selected taxa_level
n_best <- 3

############################ AGGREGATE ############################

# aggregate taxa at a specific taxonomy level
phylo_aggregated <- tax_glom(phylo_data, taxrank=taxa_level)

# get metadata for phylo_aggregated
phylo_aggregated_meta <- sample_data(phylo_aggregated)

############################ GET INFO NEEDED FOR THE PLOTS ############################

# define empty variable that will store all results for the final plot
control_and_warming_tmp <- list()
unidenti_taxa <- vector()
best_taxa <- vector()
other_taxa <- vector()

# storage for plots
all_plots <- list()

# loop through depths
for (a_depth in unique(phylo_aggregated_meta$Depth)) {

	print(paste0("Analysing ", a_depth))

	# subset by depth. remove samples that do not belong to a specific depth
	phylo_depth <- prune_samples(as.character(phylo_aggregated_meta$sampleID[which(phylo_aggregated_meta$Depth == a_depth)]), phylo_aggregated)

	# remove ASVs which are always zero for the subset under consideration
	phylo_depth <- prune_taxa(taxa_sums(phylo_depth) > 0, phylo_depth)
	
	# get metadata
	phylo_depth_meta <- sample_data(phylo_depth)

	# compute sums to obtain the contribution of all sampling sites
	# to a specific taxa, then order from most abundant to less abundant
	depth_sums <- sort(rowSums(otu_table(phylo_depth)), decreasing=T)

	# get ordered taxa names, using factors
	taxa_names <- tax_table(phylo_depth)[names(depth_sums), taxa_level]
	
	# find unidenti, best scoring taxa, and other taxa
	unidenti <- grep("_Unidentified", taxa_names)
	if (length(unidenti)>0){
		# get taxa_names
		unidenti_taxa <- taxa_names[unidenti]
		# remove unidentified taxa
		taxa_names <- taxa_names[-unidenti]
	}	
	# get best taxa from those that are identified
	best_taxa <- taxa_names[1:n_best]
	# get all other taxa, from those that are identified
	other_taxa <- taxa_names[(n_best+1):length(taxa_names)]
	
	# define some empty lists to store all results
	merged_tabs <- vector()
	a_treat_tab <- vector()
	
	# loop through treatments to compute abundances at treatment level
	for (a_treat in unique(phylo_depth_meta$Treatment)) {

		print(paste0("Analysing ", a_treat))
	
		# subset by depth
		phylo_treat <- prune_samples(as.character(phylo_depth_meta$sampleID[which(phylo_depth_meta$Treatment == a_treat)]), phylo_depth)

		# remove ASVs which are always zero for the subset under consideration
		phylo_treat <- prune_taxa(taxa_sums(phylo_treat) > 0, phylo_treat)

		for (taxa_type in c("best_taxa", "other_taxa", "unidenti_taxa")) {
		
			print(paste0("Analysing ", taxa_type))
			
			taxa_names_sub <- ""
				
			if (taxa_type=="best_taxa") {
				# get best taxa
				treat_sums <- rowSums(otu_table(phylo_treat)[which(tax_table(phylo_treat)[, taxa_level] %in% get(taxa_type)), ])
				# get taxa names
				taxa_names_sub <- tax_table(phylo_treat)[names(treat_sums), taxa_level]
			} else {
				
				# compute sums to obtain the contribution of all sampling sites
				# to a specific taxa, then order from most abundant to less abundant
				treat_sums <- sum(rowSums(otu_table(phylo_treat)[which(tax_table(phylo_treat)[, taxa_level] %in% get(taxa_type)), ]))
				# get taxa names
				if (taxa_type=="other_taxa") {
					taxa_names_sub <- "Other"
				} else {
					taxa_names_sub <- "Unidentified"
				}
			}

			# put columns together for phylo_depth samples and compute
			# the percentage of contribution each taxa is providing to
			# the total, to finally compare the taxa. by doing so, we
			# make sure the total will always sum up to 1, i.e. 100%
			a_treat_tab <- cbind.data.frame(T_level=taxa_names_sub, Treatment=a_treat, Depth=a_depth, Abundance=treat_sums, Relative=0)
			rownames(a_treat_tab) <- NULL
			colnames(a_treat_tab) <- c("T_level", "Treatment", "Depth", "Abundance", "Relative")
			merged_tabs <- rbind.data.frame(merged_tabs, a_treat_tab)
		}
		# add relative abundance once one treatment is all computed
		merged_tabs[which(merged_tabs$Treatment==a_treat), ]$Relative <- merged_tabs[which(merged_tabs$Treatment==a_treat), ]$Abundance/sum(merged_tabs[which(merged_tabs$Treatment==a_treat), ]$Abundance)
	}
	
	############################ SET FACTOR LEVELS ############################
	
	print("Organising for plots")

	# set levels in order
	last_two <- as.character(c(merged_tabs$T_level[!grepl("Other|Unidenti", merged_tabs$T_level)], merged_tabs$T_level[grepl("Other|Unidenti", merged_tabs$T_level)]))
	merged_tabs$T_level <- factor(merged_tabs$T_level, levels=unique(last_two))

	# re-set labels
	merged_tabs$Treatment <- ifelse(merged_tabs$Treatment=="W", "Warming", "Control")

	# merge depths together
	control_and_warming_tmp[[a_depth]] <- merged_tabs
	
	############################ PLOT WITH FACETS ############################
	
	print("plotting with facets")

	# get the merged_tabs plot ready
	all_plots[[a_depth]] <- merged_tabs %>%
		count(T_level, Treatment, Depth, wt=Relative, name="Relative") %>%
		ggplot() +
		geom_bar(aes(fill=T_level, y=Relative, x=Treatment), stat="identity", color="black") +
		scale_fill_viridis(discrete=T, option="turbo") +
		theme_bw() +
		theme(text=element_text(size=25, face="bold"), legend.position="none", strip.text=element_text(size=15), legend.key.height=unit(0.7,'cm'), axis.text.x=element_text(angle=50, hjust=1)) +
		guides(fill=guide_legend(ncol=1)) +
		xlab("") +
		ylab("") +
		facet_wrap(~T_level) +
		ggtitle(paste0("Taxa composition for ", a_depth, " at ", taxa_level, " level"))
}

############################ SET FACTORS LEVELS ############################

control_and_warming <- do.call(rbind.data.frame, control_and_warming_tmp)

# set levels for T_level
last_two <- as.character(c(control_and_warming$T_level[!grepl("^Other$|^Unidenti$", control_and_warming$T_level)], control_and_warming$T_level[grepl("^Other$|^Unidenti$", control_and_warming$T_level)]))
control_and_warming$T_level <- factor(control_and_warming$T_level, levels=unique(last_two))

# set levels for Depth
control_and_warming$Depth <- factor(control_and_warming$Depth, levels=unique(control_and_warming$Depth))

############################ PLOT STACKED ############################

print("plot stacked")

# get the plot ready
all_plots[["ribon"]] <- control_and_warming %>%
	count(T_level, Treatment, Depth, wt=Relative, name="Relative") %>%
	ggplot() +
	geom_bar(aes(fill=T_level, y=Relative, x=Treatment), position="fill", stat="identity", color="black") +
	scale_fill_viridis(discrete=T, option="turbo") +
	theme_bw() +
	theme(text=element_text(size=25, face="bold"), legend.position="right", legend.text=element_text(size=15), legend.key.height=unit(0.7,'cm'), axis.text.x=element_text(angle=50, hjust=1)) +
	guides(fill=guide_legend(ncol=1)) +
	xlab("") +
	ylab("") +
	ylim(0, 1) +
	facet_wrap(~Depth) +
	labs(fill=taxa_level)
