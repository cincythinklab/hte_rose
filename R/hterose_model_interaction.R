require(tidyverse)
require(ggplot2)
require(brms)
require(ggpmisc)
require(gridExtra)

####### Functions #######
empty_to_NA <- function(x) {
  ifelse(x == "", NA, x)
}

multiplot <- function(..., cols = 1) {
  plots <- list(...)
  numPlots <- length(plots)
  grid.arrange(grobs = plots, ncol = cols)
}

####### Load and prepare data #######
data_hte <- read.csv("data/rose_subset.csv")

data_hte <- data_hte %>%
	mutate_if(is.character, empty_to_NA) %>%
	mutate(d90_mortality = ifelse(d90_mortality=="Dead", 1, 0),
		   cluster1 = ifelse(subphenotype=="Subphenotype A", 1, 0),
		   cluster2 = ifelse(subphenotype=="Subphenotype B", 1, 0),
		   subphenotype_no = ifelse(subphenotype=="Subphenotype A", 1, 0),
		   group = ifelse(intervention=="NMB", 1, 0)) %>%
	filter(!is.na(subphenotype)) %>%
	as.data.frame()

####################################### Interaction Age and Severity #######################################
prior_int <- c(set_prior("normal(0, 1.5)", class = "Intercept"),
			   set_prior("normal(0, 0.5)", class = "b", coef = "group1"),
			   set_prior("normal(0, 0.1)", class = "b", coef = "subphenotype_no1"),
			   set_prior("normal(0, 0.1)", class = "b", coef = "pao2_fio2_rt"))

data_hte$subphenotype_no <- as.factor(data_hte$subphenotype_no)
data_hte$group <- as.factor(data_hte$group)

model_pfr <- brm(d90_mortality ~ subphenotype_no * pao2_fio2_rt +
								 group * pao2_fio2_rt +
								 subphenotype_no * group,
								 prior = prior_int,
								 data = data_hte, 
								 family = bernoulli(link = "logit"), 
								 warmup = 1000, 
								 iter = 11000, 
								 chains = 4, 
								 seed = 123, 
								 cores = 4)	

model_int <- brm(d90_mortality ~ subphenotype_no * pao2_fio2_rt * group,
								 prior = prior_int,
								 data = data_hte, 
								 family = bernoulli(link = "logit"), 
								 warmup = 1000, 
								 iter = 11000, 
								 chains = 4, 
								 seed = 123, 
								 cores = 4)	

tableplot <- data.frame(Group = c("OR < 1.00", "OR > 1.00"),
					    All = c(paste0(formatC(pnorm(log(0.9999), summary(model_int)$fixed[8,1], summary(model_int)$fixed[8,2], lower.tail = T) * 100, digits = 1, format = "f"), "%"), 
								paste0(formatC(pnorm(log(1.00), summary(model_int)$fixed[8,1], summary(model_int)$fixed[8,2], lower.tail = F) * 100, digits = 1, format = "f"), "%")))
colnames(tableplot) <- c("Group", "Probability")
data_tableplot <- tibble(x = 1, y = 3.3, tableplot = list(tableplot))

#Figure 4#	
data_samples <- posterior_samples(model_pfr)[,c(4,6,7)] #Group 1; Group:PF; Sub:Group#

data_graph_int <- data.frame(expand.grid(pf_ratio = factor(c(rep(40, nrow(data_samples)), 
															 rep(50, nrow(data_samples)), 
															 rep(75, nrow(data_samples)), 
															 rep(100, nrow(data_samples)), 
															 rep(125, nrow(data_samples)), 
															 rep(150, nrow(data_samples)))), 
										 subphenotype = c("B", "A"))) 

OR = c(as.matrix(exp(data_samples[1] + 40*data_samples[2])),
	   as.matrix(exp(data_samples[1] + 50*data_samples[2])),
	   as.matrix(exp(data_samples[1] + 75*data_samples[2])),
	   as.matrix(exp(data_samples[1] + 100*data_samples[2])),
	   as.matrix(exp(data_samples[1] + 125*data_samples[2])),
	   as.matrix(exp(data_samples[1] + 150*data_samples[2])),
	   as.matrix(exp(data_samples[1] + 40*data_samples[2] + data_samples[3])),
	   as.matrix(exp(data_samples[1] + 50*data_samples[2] + data_samples[3])),
	   as.matrix(exp(data_samples[1] + 75*data_samples[2] + data_samples[3])),
	   as.matrix(exp(data_samples[1] + 100*data_samples[2] + data_samples[3])),
	   as.matrix(exp(data_samples[1] + 125*data_samples[2] + data_samples[3])),
	   as.matrix(exp(data_samples[1] + 150*data_samples[2] + data_samples[3])))

data_graph_int <- cbind(data_graph_int, OR)

p1 <- data_graph_int %>%
	  group_by(pf_ratio, subphenotype) %>%
	  summarise(mOR = median(OR), p025 = quantile(OR, p=0.025), p975 = quantile(OR, p=0.975)) %>%
	  ggplot(aes(x = pf_ratio, y = mOR, color = subphenotype)) +
	  geom_point(size = 3.5, position = position_dodge(width=0.5)) +
	  geom_errorbar(aes(ymin = p025, ymax = p975), position = position_dodge(width = 0.5), width = 0.5) +
	  xlab("PaO2 / FiO2") +
	  ylab("Odds Ratio (95% CrI)") +
	  scale_y_continuous(limits = c(0,4), breaks = seq(0, 4, 0.5)) +
	  scale_colour_manual("", values=c("red","blue"), labels=c("Subphenotype B","Subphenotype A")) +
	  annotate("text", x = -Inf, y = Inf, label = "A", size=7, fontface = "bold", hjust = -0.5, vjust = 1) +
	  geom_hline(yintercept = 1.0, lty = "dashed") +
	  geom_table(data = data_tableplot, aes(x, y, label = tableplot)) +
	  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(colour = "black", face = "bold"), axis.text = element_text(colour = "black"), axis.text.x = element_text(size=10), legend.key=element_blank(), legend.background=element_blank(), legend.text=element_text(size=10), legend.position = c(0.20,0.90), plot.title = element_text(hjust = 0.5, face="bold")) +
	  theme(aspect.ratio = 1)

p2 <- data_graph_int %>%
	  group_by(pf_ratio, subphenotype) %>%
	  summarise(mOR = median(pnorm(log(0.9999), mean(log(OR)), sd(log(OR)), lower.tail= T)*100)) %>%
	  ggplot(aes(x = pf_ratio, y = mOR, fill = subphenotype)) +
	  geom_bar(stat='identity', position = 'dodge') +
	  xlab("PaO2 / FiO2") +
	  ylab("Probability of Odds Ratio < 1.00") +
	  scale_fill_manual("", values=c("red","blue"), labels=c("Subphenotype B","Subphenotype A")) +
	  annotate("text", x = -Inf, y = Inf, label = "B", size=7, fontface = "bold", hjust = -0.5, vjust = 1) +
	  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(colour = "black", face = "bold"), axis.text = element_text(colour = "black"), axis.text.x = element_text(size=10), legend.key=element_blank(), legend.background=element_blank(), legend.text=element_text(size=10), legend.position = "none", plot.title = element_text(hjust = 0.5, face="bold")) +
	  theme(aspect.ratio = 1) +
	  scale_y_continuous(limits = c(0,100), breaks = seq(0,100, 10), expand = c(0,0))
	  
p3 <- data_graph_int %>%
	  group_by(pf_ratio, subphenotype) %>%
	  summarise(mOR = median(pnorm(log(0.90), mean(log(OR)), sd(log(OR)), lower.tail= T)*100)) %>%
	  ggplot(aes(x = pf_ratio, y = mOR, fill = subphenotype)) +
	  geom_bar(stat='identity', position = 'dodge') +
	  xlab("PaO2 / FiO2") +
	  ylab("Probability of Odds Ratio \u2264 0.90") +
	  scale_fill_manual("", values=c("red","blue"), labels=c("Subphenotype B","Subphenotype A")) +
	  annotate("text", x = -Inf, y = Inf, label = "C", size=7, fontface = "bold", hjust = -0.5, vjust = 1) +
	  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(colour = "black", face = "bold"), axis.text = element_text(colour = "black"), axis.text.x = element_text(size=10), legend.key=element_blank(), legend.background=element_blank(), legend.text=element_text(size=10), legend.position = "none", plot.title = element_text(hjust = 0.5, face="bold")) +
	  theme(aspect.ratio = 1) +
	  scale_y_continuous(limits = c(0,100), breaks = seq(0,100, 10), expand = c(0,0))

p4 <- data_graph_int %>%
	  group_by(pf_ratio, subphenotype) %>%
	  summarise(mOR = median(pnorm(log(1.10), mean(log(OR)), sd(log(OR)), lower.tail= F)*100)) %>%
	  ggplot(aes(x = pf_ratio, y = mOR, fill = subphenotype)) +
	  geom_bar(stat='identity', position = 'dodge') +
	  xlab("PaO2 / FiO2") +
	  ylab("Probability of Odds Ratio > 1.10") +
	  scale_fill_manual("", values=c("red","blue"), labels=c("Subphenotype B","Subphenotype A")) +
	  annotate("text", x = -Inf, y = Inf, label = "D", size=7, fontface = "bold", hjust = -0.5, vjust = 1) +
	  theme(panel.grid.major.y = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title = element_text(colour = "black", face = "bold"), axis.text = element_text(colour = "black"), axis.text.x = element_text(size=10), legend.key=element_blank(), legend.background=element_blank(), legend.text=element_text(size=10), legend.position = "none", plot.title = element_text(hjust = 0.5, face="bold")) +
	  theme(aspect.ratio = 1) +
	  scale_y_continuous(limits = c(0,100), breaks = seq(0,100, 10), expand = c(0,0))

multiplot(p1,p3,p2,p4, cols = 2)
