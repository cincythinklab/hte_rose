require(tidyverse)
require(ggplot2)
require(ggsci)
require(beanz)
require(gridExtra)
require(ggpmisc)
require(grid)
require(ggpubr)

####### Functions #######
empty_to_NA <- function(x) {
  ifelse(x == "", NA, x)
}

####### Load and prepare data #######
data_hte <- read.csv("data/rose_subset.csv")

data_hte <- data_hte %>%
	mutate_if(is.character, empty_to_NA) %>%
	mutate(
	  d90_mortality = ifelse(d90_mortality=="Dead", 1, 0),
    cluster1 = ifelse(subphenotype=="Subphenotype A", 1, 0),
		cluster2 = ifelse(subphenotype=="Subphenotype B", 1, 0),
		group = ifelse(intervention=="NMB", 1, 0)
	) %>%
	filter(!is.na(subphenotype)) %>%
	as.data.frame()


###################################################### WEAKLY INFORMATIVE PRIORS ######################################################

####### POOLED COHORT #######
hte_cluster <- bzGetSubgrpRaw(data_hte, var.resp = "d90_mortality", var.cov = c("cluster2","cluster1"), var.trt = "group", var.censor = NULL, resptype = "binary")
print(hte_cluster)

var_estvar <- c("Estimate","Variance")


#HTE (always variance)#
rst_ns <- bzCallStan("nse", dat.sub = hte_cluster, var.estvar = var_estvar, var.cov = c("cluster2","cluster1"), 
						par.pri = c(B = 0.177, MU = 0), chains = 4, iter = 11000, warmup = 1000, seed = 1000, cores = 4)

rst_sr <- bzCallStan("srs", dat.sub = hte_cluster, var.estvar = var_estvar, var.cov = c("cluster2","cluster1"), 
					   par.pri = c(B = 2.25, C = 0.177, D = 1, MU = 0), control = list(adapt_delta = 0.999), chains = 4, iter = 11000, warmup = 1000, seed = 1000, cores = 4)

tbl_sub_wi <- bzSummary(rst_sr, ref.stan.rst = rst_ns, ref.sel.grps = 1)
tbl_sub_wi$Mean <- as.numeric(as.character(tbl_sub_wi$Mean))
tbl_sub_wi$SD <- as.numeric(as.character(tbl_sub_wi$SD))
tbl_sub_wi$Q025 <- as.numeric(as.character(tbl_sub_wi$Q025))
tbl_sub_wi$Q975 <- as.numeric(as.character(tbl_sub_wi$Q975))

table_pooled_wi <- cbind(Group = c("All Patients", "Subphenotype A", "Subphenotype B"), 
				  OR = c(formatC(exp(tbl_sub_wi$Mean[3]), digits = 2, format = "f"), formatC(exp(tbl_sub_wi$Mean[1]), digits = 2, format = "f"), formatC(exp(tbl_sub_wi$Mean[2]), digits = 2, format = "f")),
				  LCI = c(formatC(exp(tbl_sub_wi$Q025[3]), digits = 2, format = "f"), formatC(exp(tbl_sub_wi$Q025[1]), digits = 2, format = "f"), formatC(exp(tbl_sub_wi$Q025[2]), digits = 2, format = "f")),
				  HCI = c(formatC(exp(tbl_sub_wi$Q975[3]), digits = 2, format = "f"), formatC(exp(tbl_sub_wi$Q975[1]), digits = 2, format = "f"), formatC(exp(tbl_sub_wi$Q975[2]), digits = 2, format = "f")))
table_pooled_wi
bzSummaryComp(rst_sr, sel.grps = 1:2, seed = 123) #33.3%#
#Figure#
tableplot_wi <- data.frame(Group = c("\u2264 0.90", "\u2264 0.97", "< 1.00", "> 1.00"),
					    All = c(paste0(formatC(pnorm(log(0.90), tbl_sub_wi$Mean[3], tbl_sub_wi$SD[3], lower.tail = T)*100, digits = 1, format = "f"), "%"), 
									paste0(formatC(pnorm(log(0.97), tbl_sub_wi$Mean[3], tbl_sub_wi$SD[3], lower.tail = T)*100, digits = 1, format = "f"), "%"),
									paste0(formatC(pnorm(log(0.9999), tbl_sub_wi$Mean[3], tbl_sub_wi$SD[3], lower.tail = T)*100, digits = 1, format = "f"), "%"),
									paste0(formatC(pnorm(log(1.00), tbl_sub_wi$Mean[3], tbl_sub_wi$SD[3], lower.tail = F)*100, digits = 1, format = "f"), "%")),
						Elderly = c(paste0(formatC(pnorm(log(0.90), tbl_sub_wi$Mean[1], tbl_sub_wi$SD[1], lower.tail = T)*100, digits = 1, format = "f"), "%"), 
									paste0(formatC(pnorm(log(0.97), tbl_sub_wi$Mean[1], tbl_sub_wi$SD[1], lower.tail = T)*100, digits = 1, format = "f"), "%"),
									paste0(formatC(pnorm(log(0.9999), tbl_sub_wi$Mean[1], tbl_sub_wi$SD[1], lower.tail = T)*100, digits = 1, format = "f"), "%"),
									paste0(formatC(pnorm(log(1.00), tbl_sub_wi$Mean[1], tbl_sub_wi$SD[1], lower.tail = F)*100, digits = 1, format = "f"), "%")),
						Young = c(paste0(formatC(pnorm(log(0.90), tbl_sub_wi$Mean[2], tbl_sub_wi$SD[2], lower.tail = T)*100, digits = 1, format = "f"), "%"), 
									paste0(formatC(pnorm(log(0.97), tbl_sub_wi$Mean[2], tbl_sub_wi$SD[2], lower.tail = T)*100, digits = 1, format = "f"), "%"),
									paste0(formatC(pnorm(log(0.9999), tbl_sub_wi$Mean[2], tbl_sub_wi$SD[2], lower.tail = T)*100, digits = 1, format = "f"), "%"),
									paste0(formatC(pnorm(log(1.00), tbl_sub_wi$Mean[2], tbl_sub_wi$SD[2], lower.tail = F)*100, digits = 1, format = "f"), "%")))
colnames(tableplot_wi) <- c("OR", "All\r\nPatients", "Subphenotype\r\nA", "Subphenotype\r\nB")
data_table_wi <- tibble(x = 1.7, y = 2.5, tableplot_wi = list(tableplot_wi))

p1 <- ggplot(data = data.frame(x = c(-1.5,1.7)), aes(x)) +
	stat_function(fun = dnorm, args = list(mean = as.numeric(as.character(tbl_sub_wi$Mean[1])), sd = as.numeric(as.character(tbl_sub_wi$SD[1]))), size = 1, aes(color = "Subphenotype A")) +
	stat_function(fun = dnorm, args = list(mean = as.numeric(as.character(tbl_sub_wi$Mean[2])), sd = as.numeric(as.character(tbl_sub_wi$SD[2]))), size = 1, aes(color = "Subphenotype B")) +
	stat_function(fun = dnorm, args = list(mean = as.numeric(as.character(tbl_sub_wi$Mean[3])), sd = as.numeric(as.character(tbl_sub_wi$SD[3]))), size = 1, aes(color = "All patients")) +
	geom_vline(xintercept = 0,linetype = 2) +
	scale_colour_manual("", values=c("black","blue","red"), labels=c("All Patients","Subphenotype A","Subphenotype B")) +
	theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(colour = "black", size = 18), axis.text = element_text(colour = "black", size = 13), legend.key=element_blank(), legend.background=element_blank(), legend.text=element_text(size=13), legend.position = "bottom", plot.title = element_text(hjust = 0.5, face="bold", size = 20)) +
	xlab("Odds Ratio for 90-Day Mortality") +
	ylab("Density") +
	scale_x_continuous(limits=c(-1.5,1.7), labels = round(c(exp(seq(-1.5, 1.7, by = 0.4))), 2), breaks = seq(-1.5, 1.7, by = 0.4)) +
	scale_y_continuous(expand = c(0,0)) +
	annotate("text", x = -1.3, y = 0.5, label = "Cisatracurium\r\nBetter", size = 5) +
	annotate("text", x = 1.0, y = 0.5, label = "Cisatracurium\r\nWorse", size = 5) +
	geom_table(data = data_table_wi, aes(x, y, label = tableplot_wi), size = 4) +
	geom_area(stat = "function", fun = dnorm, args = list(mean = as.numeric(as.character(tbl_sub_wi$Mean[1])), sd = as.numeric(as.character(tbl_sub_wi$SD[1]))), fill = "darkblue", xlim = c(0,1.7), alpha = 0.5) +
	geom_area(stat = "function", fun = dnorm, args = list(mean = as.numeric(as.character(tbl_sub_wi$Mean[1])), sd = as.numeric(as.character(tbl_sub_wi$SD[1]))), fill = "lightblue", xlim = c(-1.5,0), alpha = 0.5) +
	geom_area(stat = "function", fun = dnorm, args = list(mean = as.numeric(as.character(tbl_sub_wi$Mean[2])), sd = as.numeric(as.character(tbl_sub_wi$SD[2]))), fill = "red", xlim = c(0,1.7), alpha = 0.5) +
	geom_area(stat = "function", fun = dnorm, args = list(mean = as.numeric(as.character(tbl_sub_wi$Mean[2])), sd = as.numeric(as.character(tbl_sub_wi$SD[2]))), fill = "pink", xlim = c(-1.5,0), alpha = 0.5) +
	geom_area(stat = "function", fun = dnorm, args = list(mean = as.numeric(as.character(tbl_sub_wi$Mean[3])), sd = as.numeric(as.character(tbl_sub_wi$SD[3]))), fill = "gray30", xlim = c(0,1.7), alpha = 0.5) +
	geom_area(stat = "function", fun = dnorm, args = list(mean = as.numeric(as.character(tbl_sub_wi$Mean[3])), sd = as.numeric(as.character(tbl_sub_wi$SD[3]))), fill = "gray90", xlim = c(-1.5,0), alpha = 0.5)	

p1


###################################################### OPTIMISTIC PRIORS ######################################################

####### POOLED COHORT #######
hte_cluster <- bzGetSubgrpRaw(data_hte, var.resp = "d90_mortality", var.cov = c("cluster2","cluster1"), var.trt = "group", var.censor = NULL, resptype = "binary")
print(hte_cluster)

var_estvar <- c("Estimate","Variance")


#HTE (always variance)#
rst_ns <- bzCallStan("nse", dat.sub = hte_cluster, var.estvar = var_estvar, var.cov = c("cluster2","cluster1"), 
						par.pri = c(B = 0.030, MU = -0.287), chains = 4, iter = 11000, warmup = 1000, seed = 1000, cores = 4)

rst_sr <- bzCallStan("srs", dat.sub = hte_cluster, var.estvar = var_estvar, var.cov = c("cluster2","cluster1"), 
					   par.pri = c(B = 2.25, C = 0.030, D = 1, MU = -0.287), control = list(adapt_delta = 0.999), chains = 4, iter = 11000, warmup = 1000, seed = 1000, cores = 4)

tbl_sub_wi <- bzSummary(rst_sr, ref.stan.rst = rst_ns, ref.sel.grps = 1)
tbl_sub_wi$Mean <- as.numeric(as.character(tbl_sub_wi$Mean))
tbl_sub_wi$SD <- as.numeric(as.character(tbl_sub_wi$SD))
tbl_sub_wi$Q025 <- as.numeric(as.character(tbl_sub_wi$Q025))
tbl_sub_wi$Q975 <- as.numeric(as.character(tbl_sub_wi$Q975))

table_pooled_op <- cbind(Group = c("All Patients", "Subphenotype A", "Subphenotype B"), 
				  OR = c(formatC(exp(tbl_sub_wi$Mean[3]), digits = 2, format = "f"), formatC(exp(tbl_sub_wi$Mean[1]), digits = 2, format = "f"), formatC(exp(tbl_sub_wi$Mean[2]), digits = 2, format = "f")),
				  LCI = c(formatC(exp(tbl_sub_wi$Q025[3]), digits = 2, format = "f"), formatC(exp(tbl_sub_wi$Q025[1]), digits = 2, format = "f"), formatC(exp(tbl_sub_wi$Q025[2]), digits = 2, format = "f")),
				  HCI = c(formatC(exp(tbl_sub_wi$Q975[3]), digits = 2, format = "f"), formatC(exp(tbl_sub_wi$Q975[1]), digits = 2, format = "f"), formatC(exp(tbl_sub_wi$Q975[2]), digits = 2, format = "f")))
table_pooled_op
bzSummaryComp(rst_sr, sel.grps = 1:2, seed = 123) #36.4%#
#Figure#
tableplot_wi <- data.frame(Group = c("\u2264 0.90", "\u2264 0.97", "< 1.00", "> 1.00"),
					    All = c(paste0(formatC(pnorm(log(0.90), tbl_sub_wi$Mean[3], tbl_sub_wi$SD[3], lower.tail = T)*100, digits = 1, format = "f"), "%"), 
									paste0(formatC(pnorm(log(0.97), tbl_sub_wi$Mean[3], tbl_sub_wi$SD[3], lower.tail = T)*100, digits = 1, format = "f"), "%"),
									paste0(formatC(pnorm(log(0.9999), tbl_sub_wi$Mean[3], tbl_sub_wi$SD[3], lower.tail = T)*100, digits = 1, format = "f"), "%"),
									paste0(formatC(pnorm(log(1.00), tbl_sub_wi$Mean[3], tbl_sub_wi$SD[3], lower.tail = F)*100, digits = 1, format = "f"), "%")),
						Elderly = c(paste0(formatC(pnorm(log(0.90), tbl_sub_wi$Mean[1], tbl_sub_wi$SD[1], lower.tail = T)*100, digits = 1, format = "f"), "%"), 
									paste0(formatC(pnorm(log(0.97), tbl_sub_wi$Mean[1], tbl_sub_wi$SD[1], lower.tail = T)*100, digits = 1, format = "f"), "%"),
									paste0(formatC(pnorm(log(0.9999), tbl_sub_wi$Mean[1], tbl_sub_wi$SD[1], lower.tail = T)*100, digits = 1, format = "f"), "%"),
									paste0(formatC(pnorm(log(1.00), tbl_sub_wi$Mean[1], tbl_sub_wi$SD[1], lower.tail = F)*100, digits = 1, format = "f"), "%")),
						Young = c(paste0(formatC(pnorm(log(0.90), tbl_sub_wi$Mean[2], tbl_sub_wi$SD[2], lower.tail = T)*100, digits = 1, format = "f"), "%"), 
									paste0(formatC(pnorm(log(0.97), tbl_sub_wi$Mean[2], tbl_sub_wi$SD[2], lower.tail = T)*100, digits = 1, format = "f"), "%"),
									paste0(formatC(pnorm(log(0.9999), tbl_sub_wi$Mean[2], tbl_sub_wi$SD[2], lower.tail = T)*100, digits = 1, format = "f"), "%"),
									paste0(formatC(pnorm(log(1.00), tbl_sub_wi$Mean[2], tbl_sub_wi$SD[2], lower.tail = F)*100, digits = 1, format = "f"), "%")))
colnames(tableplot_wi) <- c("OR", "All\r\nPatients", "Subphenotype A", "Subphenotype B")
data_table_wi <- tibble(x = 1.5, y = 3, tableplot_wi = list(tableplot_wi))

p2 <- ggplot(data = data.frame(x = c(-1.5,1.5)), aes(x)) +
	stat_function(fun = dnorm, args = list(mean = as.numeric(as.character(tbl_sub_wi$Mean[1])), sd = as.numeric(as.character(tbl_sub_wi$SD[1]))), size = 1, aes(color = "Subphenotype A")) +
	stat_function(fun = dnorm, args = list(mean = as.numeric(as.character(tbl_sub_wi$Mean[2])), sd = as.numeric(as.character(tbl_sub_wi$SD[2]))), size = 1, aes(color = "Subphenotype B")) +
	stat_function(fun = dnorm, args = list(mean = as.numeric(as.character(tbl_sub_wi$Mean[3])), sd = as.numeric(as.character(tbl_sub_wi$SD[3]))), size = 1, aes(color = "All patients")) +
	geom_vline(xintercept = 0,linetype = 2) +
	scale_colour_manual("", values=c("black","blue","red"), labels=c("All Patients","Subphenotype A","Subphenotype B")) +
	theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(colour = "black", size = 18), axis.text = element_text(colour = "black", size = 13), legend.key=element_blank(), legend.background=element_blank(), legend.text=element_text(size=13), legend.position = "bottom", plot.title = element_text(hjust = 0.5, face="bold", size = 20)) +
	xlab("Odds Ratio for 90-Day Mortality") +
	ylab("Density") +
	scale_x_continuous(limits=c(-1.5,1.5), labels = round(c(exp(seq(-1.5, 1.5, by = 0.2))), 2), breaks = seq(-1.5, 1.5, by = 0.2)) +
	scale_y_continuous(expand = c(0,0)) +
	annotate("text", x = -1.3, y = 0.5, label = "Cisatracurium\r\nBetter", size = 5) +
	annotate("text", x = 1.0, y = 0.5, label = "Cisatracurium\r\nWorse", size = 5) +
	geom_table(data = data_table_wi, aes(x, y, label = tableplot_wi), size = 4) +
	geom_area(stat = "function", fun = dnorm, args = list(mean = as.numeric(as.character(tbl_sub_wi$Mean[1])), sd = as.numeric(as.character(tbl_sub_wi$SD[1]))), fill = "darkblue", xlim = c(0,1.5), alpha = 0.5) +
	geom_area(stat = "function", fun = dnorm, args = list(mean = as.numeric(as.character(tbl_sub_wi$Mean[1])), sd = as.numeric(as.character(tbl_sub_wi$SD[1]))), fill = "lightblue", xlim = c(-1.5,0), alpha = 0.5) +
	geom_area(stat = "function", fun = dnorm, args = list(mean = as.numeric(as.character(tbl_sub_wi$Mean[2])), sd = as.numeric(as.character(tbl_sub_wi$SD[2]))), fill = "red", xlim = c(0,1.5), alpha = 0.5) +
	geom_area(stat = "function", fun = dnorm, args = list(mean = as.numeric(as.character(tbl_sub_wi$Mean[2])), sd = as.numeric(as.character(tbl_sub_wi$SD[2]))), fill = "pink", xlim = c(-1.5,0), alpha = 0.5) +
	geom_area(stat = "function", fun = dnorm, args = list(mean = as.numeric(as.character(tbl_sub_wi$Mean[3])), sd = as.numeric(as.character(tbl_sub_wi$SD[3]))), fill = "gray30", xlim = c(0,1.5), alpha = 0.5) +
	geom_area(stat = "function", fun = dnorm, args = list(mean = as.numeric(as.character(tbl_sub_wi$Mean[3])), sd = as.numeric(as.character(tbl_sub_wi$SD[3]))), fill = "gray90", xlim = c(-1.5,0), alpha = 0.5)	

p2
###################################################### PESSIMISTIC PRIORS ######################################################

####### POOLED COHORT #######
hte_cluster <- bzGetSubgrpRaw(data_hte, var.resp = "d90_mortality", var.cov = c("cluster2","cluster1"), var.trt = "group", var.censor = NULL, resptype = "binary")
print(hte_cluster)

var_estvar <- c("Estimate","Variance")


#HTE (always variance)#
rst_ns <- bzCallStan("nse", dat.sub = hte_cluster, var.estvar = var_estvar, var.cov = c("cluster2","cluster1"), 
						par.pri = c(B = 0.012, MU = 0.183), chains = 4, iter = 11000, warmup = 1000, seed = 1000, cores = 4)

rst_sr <- bzCallStan("srs", dat.sub = hte_cluster, var.estvar = var_estvar, var.cov = c("cluster2","cluster1"), 
					   par.pri = c(B = 2.25, C = 0.012, D = 1, MU = 0.183), control = list(adapt_delta = 0.999), chains = 4, iter = 11000, warmup = 1000, seed = 1000, cores = 4)

tbl_sub_wi <- bzSummary(rst_sr, ref.stan.rst = rst_ns, ref.sel.grps = 1)
tbl_sub_wi$Mean <- as.numeric(as.character(tbl_sub_wi$Mean))
tbl_sub_wi$SD <- as.numeric(as.character(tbl_sub_wi$SD))
tbl_sub_wi$Q025 <- as.numeric(as.character(tbl_sub_wi$Q025))
tbl_sub_wi$Q975 <- as.numeric(as.character(tbl_sub_wi$Q975))

table_pooled_pe <- cbind(Group = c("All Patients", "Subphenotype A", "Subphenotype B"), 
				  OR = c(formatC(exp(tbl_sub_wi$Mean[3]), digits = 2, format = "f"), formatC(exp(tbl_sub_wi$Mean[1]), digits = 2, format = "f"), formatC(exp(tbl_sub_wi$Mean[2]), digits = 2, format = "f")),
				  LCI = c(formatC(exp(tbl_sub_wi$Q025[3]), digits = 2, format = "f"), formatC(exp(tbl_sub_wi$Q025[1]), digits = 2, format = "f"), formatC(exp(tbl_sub_wi$Q025[2]), digits = 2, format = "f")),
				  HCI = c(formatC(exp(tbl_sub_wi$Q975[3]), digits = 2, format = "f"), formatC(exp(tbl_sub_wi$Q975[1]), digits = 2, format = "f"), formatC(exp(tbl_sub_wi$Q975[2]), digits = 2, format = "f")))
table_pooled_pe
bzSummaryComp(rst_sr, sel.grps = 1:2, seed = 123) #37.4%#
#Figure#
tableplot_wi <- data.frame(Group = c("\u2264 0.90", "\u2264 0.97", "< 1.00", "> 1.00"),
					    All = c(paste0(formatC(pnorm(log(0.90), tbl_sub_wi$Mean[3], tbl_sub_wi$SD[3], lower.tail = T)*100, digits = 1, format = "f"), "%"), 
									paste0(formatC(pnorm(log(0.97), tbl_sub_wi$Mean[3], tbl_sub_wi$SD[3], lower.tail = T)*100, digits = 1, format = "f"), "%"),
									paste0(formatC(pnorm(log(0.9999), tbl_sub_wi$Mean[3], tbl_sub_wi$SD[3], lower.tail = T)*100, digits = 1, format = "f"), "%"),
									paste0(formatC(pnorm(log(1.00), tbl_sub_wi$Mean[3], tbl_sub_wi$SD[3], lower.tail = F)*100, digits = 1, format = "f"), "%")),
						Elderly = c(paste0(formatC(pnorm(log(0.90), tbl_sub_wi$Mean[1], tbl_sub_wi$SD[1], lower.tail = T)*100, digits = 1, format = "f"), "%"), 
									paste0(formatC(pnorm(log(0.97), tbl_sub_wi$Mean[1], tbl_sub_wi$SD[1], lower.tail = T)*100, digits = 1, format = "f"), "%"),
									paste0(formatC(pnorm(log(0.9999), tbl_sub_wi$Mean[1], tbl_sub_wi$SD[1], lower.tail = T)*100, digits = 1, format = "f"), "%"),
									paste0(formatC(pnorm(log(1.00), tbl_sub_wi$Mean[1], tbl_sub_wi$SD[1], lower.tail = F)*100, digits = 1, format = "f"), "%")),
						Young = c(paste0(formatC(pnorm(log(0.90), tbl_sub_wi$Mean[2], tbl_sub_wi$SD[2], lower.tail = T)*100, digits = 1, format = "f"), "%"), 
									paste0(formatC(pnorm(log(0.97), tbl_sub_wi$Mean[2], tbl_sub_wi$SD[2], lower.tail = T)*100, digits = 1, format = "f"), "%"),
									paste0(formatC(pnorm(log(0.9999), tbl_sub_wi$Mean[2], tbl_sub_wi$SD[2], lower.tail = T)*100, digits = 1, format = "f"), "%"),
									paste0(formatC(pnorm(log(1.00), tbl_sub_wi$Mean[2], tbl_sub_wi$SD[2], lower.tail = F)*100, digits = 1, format = "f"), "%")))
colnames(tableplot_wi) <- c("OR", "All\r\nPatients", "Subphenotype A", "Subphenotype B")
data_table_wi <- tibble(x = 1.7, y = 4, tableplot_wi = list(tableplot_wi))

p3 <- ggplot(data = data.frame(x = c(-1.5,1.7)), aes(x)) +
	stat_function(fun = dnorm, args = list(mean = as.numeric(as.character(tbl_sub_wi$Mean[1])), sd = as.numeric(as.character(tbl_sub_wi$SD[1]))), size = 1, aes(color = "Subphenotype A")) +
	stat_function(fun = dnorm, args = list(mean = as.numeric(as.character(tbl_sub_wi$Mean[2])), sd = as.numeric(as.character(tbl_sub_wi$SD[2]))), size = 1, aes(color = "Subphenotype B")) +
	stat_function(fun = dnorm, args = list(mean = as.numeric(as.character(tbl_sub_wi$Mean[3])), sd = as.numeric(as.character(tbl_sub_wi$SD[3]))), size = 1, aes(color = "All patients")) +
	geom_vline(xintercept = 0,linetype = 2) +
	scale_colour_manual("", values=c("black","blue","red"), labels=c("All Patients","Subphenotype A","Subphenotype B")) +
	theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(colour = "black", size = 18), axis.text = element_text(colour = "black", size = 13), legend.key=element_blank(), legend.background=element_blank(), legend.text=element_text(size=13), legend.position = "bottom", plot.title = element_text(hjust = 0.5, face="bold", size = 20)) +
	xlab("Odds Ratio for 90-Day Mortality") +
	ylab("Density") +
	scale_x_continuous(limits=c(-1.5,1.7), labels = round(c(exp(seq(-1.5, 1.7, by = 0.2))), 2), breaks = seq(-1.5, 1.7, by = 0.2)) +
	scale_y_continuous(expand = c(0,0)) +
	annotate("text", x = -1.3, y = 0.5, label = "Cisatracurium\r\nBetter", size = 5) +
	annotate("text", x = 1.0, y = 0.5, label = "Cisatracurium\r\nWorse", size = 5) +
	geom_table(data = data_table_wi, aes(x, y, label = tableplot_wi), size = 4) +
	geom_area(stat = "function", fun = dnorm, args = list(mean = as.numeric(as.character(tbl_sub_wi$Mean[1])), sd = as.numeric(as.character(tbl_sub_wi$SD[1]))), fill = "darkblue", xlim = c(0,1.7), alpha = 0.5) +
	geom_area(stat = "function", fun = dnorm, args = list(mean = as.numeric(as.character(tbl_sub_wi$Mean[1])), sd = as.numeric(as.character(tbl_sub_wi$SD[1]))), fill = "lightblue", xlim = c(-1.5,0), alpha = 0.5) +
	geom_area(stat = "function", fun = dnorm, args = list(mean = as.numeric(as.character(tbl_sub_wi$Mean[2])), sd = as.numeric(as.character(tbl_sub_wi$SD[2]))), fill = "red", xlim = c(0,1.7), alpha = 0.5) +
	geom_area(stat = "function", fun = dnorm, args = list(mean = as.numeric(as.character(tbl_sub_wi$Mean[2])), sd = as.numeric(as.character(tbl_sub_wi$SD[2]))), fill = "pink", xlim = c(-1.5,0), alpha = 0.5) +
	geom_area(stat = "function", fun = dnorm, args = list(mean = as.numeric(as.character(tbl_sub_wi$Mean[3])), sd = as.numeric(as.character(tbl_sub_wi$SD[3]))), fill = "gray30", xlim = c(0,1.7), alpha = 0.5) +
	geom_area(stat = "function", fun = dnorm, args = list(mean = as.numeric(as.character(tbl_sub_wi$Mean[3])), sd = as.numeric(as.character(tbl_sub_wi$SD[3]))), fill = "gray90", xlim = c(-1.5,0), alpha = 0.5)

p3
