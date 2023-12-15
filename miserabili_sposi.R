library(tidyverse)
library(MASS) 
library(reshape2) 
library(reshape) 

data = read.csv("promessisposi/promessisposi.csv")
data2 = read.csv("imiserabili/imiserabili.csv")
L = c(473,466,254,271)
N = c(79,79,77,77)
k = 2 * L / N 
p = k / (N-1)
average_path_length = c(2.052,1.99,2.641,2.434)
average_clustering_coefficient = c(0.768,0.159,0.736,0.115)
diameter = c(4,3,5,5)

connected_components = c(NA,1,NA,1)
NG = c(NA,79,NA,77)


stats = data.frame(L, N, 
                  average_path_length, average_clustering_coefficient, diameter, 
                  connected_components, NG, k, p)
rownames(stats) <- c("ipromessisposi","random_sposi","imiserabili","random_miserabili")
colnames(stats) <- c("L","N","avg_PL","avg_cc","d","ncc","NG","k","p")


#I PROMESSI SPOSI
ggplot(data = data, aes(x = Degree, y = clustering)) +
  geom_point() + ggtitle("Figura 8: Regression Line Clustering Coefficient vs Degree 'I Promessi Sposi") +
  xlab("Degree") + ylab("Clustering coefficient") + 
  theme(plot.title = element_text(hjust = 0.5),axis.text.x=element_text(size=14),axis.text.y=element_text(size=15),
        text = element_text(size = 14)) +
  geom_smooth(method="lm" ) 
ggsave("promessisposi/regression_clustering_vs_degree.png")
model <- lm(clustering ~ Degree, data = data)
model
summary(model)

ggplot(data = data, aes(x = Degree, y = triangles)) +
  geom_point() + ggtitle("promessi sposi") +
  xlab("degree") + ylab("numero triangoli") + theme(plot.title = element_text(hjust = 0.5))

ggplot(data = data, aes(x = Degree, y = triangles)) +
  geom_point() + ggtitle("Figura 10: Regression Line Num. Triangles vs Degree 'I Promessi Sposi'") +
  xlab("Degree") + ylab("Numero di triangoli") + 
  theme(plot.title = element_text(hjust = 0.5),axis.text.x=element_text(size=14),axis.text.y=element_text(size=15),
        text = element_text(size = 14)) +
  geom_smooth(method="lm" ) 
ggsave("promessisposi/regression_triangles_vs_degree.png")

ggplot(data, aes(x=clustering)) + geom_histogram(binwidth=0.05) + 
  xlab("Clustering Coefficient") + ylab("Frequency") + 
  ggtitle("Figure 12: Clustering Coefficient distribution 'I Promessi Sposi'") + 
  theme(plot.title = element_text(hjust = 0.5),axis.text.x=element_text(size=14),axis.text.y=element_text(size=15),
        text = element_text(size = 14)) 
ggsave("promessisposi/clusteringhist.png")

ggplot(data, aes(Degree)) + geom_bar() + 
  xlab("Degree") + ylab("Frequency") + 
  ggtitle("Figure 14: Degree distribution 'I Promessi Sposi'") + 
  theme(plot.title = element_text(hjust = 0.5),axis.text.x=element_text(size=14),axis.text.y=element_text(size=15),
        text = element_text(size = 14)) 
ggsave("promessisposi/degreebar.png")

summary(data["Degree"])

frequency_degree = data.frame(table(cut(data$Degree,breaks=c(0,1,3,7,15,31,64))))

colnames(frequency_degree) = c("Bin","Frequency")

write.csv(frequency_degree,"promessiprosi/frequency_degree.csv")

lambda = 1 / mean(data$Degree)


#degree distribution goodness of fit with exponential
expected_probabilities = data.frame(dexp(x=1:64,rate=lambda) * length(data$Degree))
b = c(1,2,4,8,16,32,64)
freq <- c()
for (i in 1:6){
  a=b[i]
  c=b[i+1]-1
  print(a)
  print(c)
  print(expected_probabilities[a:c,])
  freq <- append(freq,sum(expected_probabilities[a:c,]))
}
frequency_degree$Expected_frequency = freq
df_combined <- cbind(frequency_degree$Frequency, frequency_degree$Expected_frequency)
colnames(df_combined) <- c("Frequency", "Expected_frequency")

df_long <- melt(cbind(frequency_degree["Bin"], df_combined),
                id.vars = "Bin",
                variable.name = "freq_type",
                value.name = "count")

# Create the barplot using ggplot2
ggplot(df_long, aes(x = Bin, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Figura 16: Degree expected vs observed", x = "Degree", y = "Count") +
  theme(plot.title = element_text(hjust = 0.5),axis.text.x=element_text(size=14),axis.text.y=element_text(size=15),
        text = element_text(size = 14)) 
ggsave("promessisposi/degree_expected_vs_obs.png")

frequency_degree$x_squared=(frequency_degree$Frequency-frequency_degree$Expected_frequency)**2/frequency_degree$Expected_frequency
sum(frequency_degree$x_squared)

x_squared = chisq.test(frequency_degree$Frequency, p = frequency_degree$Expected_frequency, rescale.p = TRUE)$statistic
x_squared
pchisq(q = x_squared, df = 4, lower.tail = FALSE)

#degree distribution goodness of fit with poisson
expected_probabilities = data.frame(dpois(x=1:64,lambda=mean(data$Degree)) * length(data$Degree))
b = c(1,2,4,8,16,32,64)
freq <- c()
for (i in 1:6){
  a=b[i]
  c=b[i+1]-1
  print(a)
  print(c)
  print(expected_probabilities[a:c,])
  freq <- append(freq,sum(expected_probabilities[a:c,]))
}
frequency_degree$Expected_frequency_pois = freq

df_combined_pois <- cbind(frequency_degree$Frequency, frequency_degree$Expected_frequency_pois)
colnames(df_combined_pois) <- c("Frequency", "Expected_frequency")

df_long_pois <- melt(cbind(frequency_degree["Bin"], df_combined_pois),
                id.vars = "Bin",
                variable.name = "freq_type",
                value.name = "count")

# Create the barplot using ggplot2
ggplot(df_long_pois, aes(x = Bin, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Degree expected vs observed", x = "Degree", y = "Count") +
  theme_bw()

frequency_degree$x_squared_pois=(frequency_degree$Frequency-frequency_degree$Expected_frequency_pois)**2/frequency_degree$Expected_frequency_pois
sum(frequency_degree$x_squared[0:3])

x_squared = chisq.test(frequency_degree$Frequency, p = frequency_degree$Expected_frequency_pois, rescale.p = TRUE)$statistic

pchisq(q = x_squared, df = 1, lower.tail = FALSE)


#degree vs pk log bins
frequency_degree$Pk = frequency_degree$Frequency / length(data$Degree)

ggplot(data = frequency_degree, aes(x = Bin, y = Pk),) +
  geom_point() + ggtitle("promessi sposi") +
  xlab("Degree") + ylab("Pk") + theme(plot.title = element_text(hjust = 0.5)) 

#degree vs pk log
probability_degree = data.frame(table(data$Degree) / length(data$Degree))
colnames(probability_degree) = c("Degree","Pk")
indx <- sapply(probability_degree, is.factor)
probability_degree[indx] <- lapply(probability_degree[indx], function(x) as.numeric(as.character(x)))
probability_degree$log_degree = log10(probability_degree$Degree)
probability_degree$log_Pk= log10(probability_degree$Pk)

ggplot(data = probability_degree, aes(x = log10(Degree), y = log10(Pk))) +
  geom_point() + ggtitle("promessi sposi") +
  xlab("Degree") + ylab("Pk") + theme(plot.title = element_text(hjust = 0.5))


#pk vs dergee log binning and log axes
mids = c(1,2.5,5.5,11.5,23.5,47.5)

hist_deg <- hist(data$Degree, breaks=c(0,1,3,7,15,31,64))
bins <- 2^seq(log2(min(data$Degree)), log2(max(data$Degree)), length.out = 6)
avg_deg <- sapply(1:length(hist_deg$counts), function(i) {
  if (i==length(hist_deg$counts)){
    mean(data$Degree[hist_deg$mids[i] <= data$Degree & data$Degree])
  }
  else{
    mean(data$Degree[hist_deg$mids[i] <= data$Degree & data$Degree < hist_deg$mids[i+1]])
  }
})

ggplot(data.frame(k = mids, 
                  Pk = hist_deg$density),
       aes(x = k, y = Pk)) +
  geom_point() +
  scale_x_continuous(trans="log10") +
  scale_y_continuous(trans="log10") +
  labs(x = "Degree (k)", y = "Fraction of nodes (Pk)")


#I MISERABILI
ggplot(data = data2, aes(x = Degree, y = clustering),) +
  geom_point() + ggtitle("miserabili") +
  xlab("degree") + ylab("clustering coefficient") + theme(plot.title = element_text(hjust = 0.5))

ggplot(data = data2, aes(x = Degree, y = clustering)) +
  geom_point() + ggtitle("Figura 9: Regression Line Clustering Coefficient vs Degree 'I Miserabili'") +
  xlab("Degree") + ylab("Clustering coefficient") + 
  theme(plot.title = element_text(hjust = 0.5),axis.text.x=element_text(size=14),axis.text.y=element_text(size=15),
        text = element_text(size = 14)) +
  geom_smooth(method="lm" ) 
ggsave("imiserabili/regression_clustering_vs_degree.png")

ggplot(data = data2, aes(x = Degree, y = triangles)) +
  geom_point() + ggtitle("Figura 11: Regression Line Num. Triangles vs Degree 'I Miserabili'") +
  xlab("Degree") + ylab("Numero di triangoli") + 
  theme(plot.title = element_text(hjust = 0.5),axis.text.x=element_text(size=14),axis.text.y=element_text(size=15),
        text = element_text(size = 14)) +
  geom_smooth(method="lm" ) 
ggsave("imiserabili/regression_triangles_vs_degree.png")

ggplot(data2, aes(x=clustering)) + geom_histogram(binwidth=0.05) + 
  xlab("Clustering Coefficient") + ylab("Frequency") + 
  ggtitle("Figure 13: Clustering Coefficient distribution 'I Miserabili'") + 
  theme(plot.title = element_text(hjust = 0.5),axis.text.x=element_text(size=14),axis.text.y=element_text(size=15),
      text = element_text(size = 14)) 
ggsave("imiserabili/clusteringhist.png")

ggplot(data2, aes(Degree)) + geom_bar() + 
  xlab("Degree") + ylab("Frequency") + 
  ggtitle("Figure 15: Degree distribution 'I Miserabili'") + 
  theme(plot.title = element_text(hjust = 0.5),axis.text.x=element_text(size=14),axis.text.y=element_text(size=15),
        text = element_text(size = 14)) 
ggsave("imiserabili/degreebar.png")

frequency_degree = data.frame(table(cut(data2$Degree,breaks=c(0,1,3,7,15,31,64))))

colnames(frequency_degree) = c("Bin","Frequency")

write.csv(frequency_degree,"imiserabili/frequency_degree.csv")

lambda = 1 / mean(data2$Degree)

#degree distribution goodness of fit with exponential
expected_probabilities = data.frame(dexp(x=1:64,rate=lambda) * length(data2$Degree))
b = c(1,2,4,8,16,32,64)
freq <- c()
for (i in 1:6){
  a=b[i]
  c=b[i+1]-1
  print(a)
  print(c)
  print(expected_probabilities[a:c,])
  freq <- append(freq,sum(expected_probabilities[a:c,]))
}
frequency_degree$Expected_frequency = freq
df_combined <- cbind(frequency_degree$Frequency, frequency_degree$Expected_frequency)
colnames(df_combined) <- c("Frequency", "Expected_frequency")

df_long <- melt(cbind(frequency_degree["Bin"], df_combined),
                id.vars = "Bin",
                variable.name = "freq_type",
                value.name = "count")

# Create the barplot using ggplot2
ggplot(df_long, aes(x = Bin, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Figura 17: Degree expected vs observed", x = "Degree", y = "Count") +
  theme(plot.title = element_text(hjust = 0.5),axis.text.x=element_text(size=14),axis.text.y=element_text(size=15),
        text = element_text(size = 14)) 
ggsave("imiserabili/degree_expected_vs_obs.png")

frequency_degree$x_squared=(frequency_degree$Frequency-frequency_degree$Expected_frequency)**2/frequency_degree$Expected_frequency
sum(frequency_degree$x_squared)

x_squared = chisq.test(frequency_degree$Frequency, p = frequency_degree$Expected_frequency, rescale.p = TRUE)$statistic
x_squared
pchisq(q = x_squared, df = 3, lower.tail = FALSE)