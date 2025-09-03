#########################
##########蜘蛛16S########
#########################
##########by WJ #########


library(ggplot2)
library(plotly)
library(Rtsne)
library(ggpubr)
library(uwot)
library("vegan")
library(ggsci)
library(stringr)
library(compositions)
library(ggrepel)


# 清空当前环境
rm(list = ls())
# 设置工作路径
setwd("E:/Project/13-蜘蛛微生物/0-data/")

###统计宿主样本情况###
sample_info <- read.csv("sample_info.csv", header = T)

grouped_data <- aggregate(sample_info[, 1], by = list(species = sample_info$Species, province = sample_info$province), FUN = length)
colnames(grouped_data)[3] <- "count"
filtered_data <- grouped_data[grouped_data$count > 5, ]
sort(table(filtered_data$species))

###整理微生物名字###
library(stringr)
library()
# 清空当前环境
rm(list = ls())
# 设置工作路径
setwd("E:/Project/13-蜘蛛微生物/0-data/")

data <- read.csv("Species_abundance_all.csv", header = T, row.names = 1)
data$new <- str_split_fixed(rownames(data), ";", Inf)[,ncol(str_split_fixed(rownames(data),";",Inf))] 
filtered_data <- data %>% 
  filter(str_detect(new, "^s__") & new != "s__")

rownames(filtered_data) <- filtered_data$new
filtered_data <- filtered_data[, !names(filtered_data) %in% "new", drop = FALSE]
#去除都为0的行
filtered_data <- filtered_data[rowSums(filtered_data) != 0,]

write.csv(filtered_data,"Species_abundance_all.csv")


#提取SC/GZ样品做分析
rm(list = ls())
level <- "Genus"
level <- "Species"
level <- "Phylum"
data <- read.csv(file = paste0("./", level, "_abundance_all.csv"), row.name = 1, header = T, fileEncoding = "UTF-8")
sample_info <- read.csv("./sample_info.csv", header = T)


group <- "Sichuan"
group <- "Guizhou"
group_sample <- sample_info$sample_name[sample_info$province == group]
group_abd <- data[,colnames(data) %in% group_sample]
#去除都为0的行
group_abd <- group_abd[rowSums(group_abd) != 0,]
#取出样品信息
# sample <- sample_info[sample_info$province == group,] 
# sample_info <- sample

#相对丰度
data_rel <- sweep(data, 2, colSums(data), FUN="/")
group_rel <- sweep(group_abd, 2, colSums(group_abd), FUN="/")

write.csv(data_rel,paste0(level,"_abundance_all_rel.csv"))

write.csv(group_abd,paste0(level,"_abundance_", group, ".csv"))
write.csv(group_rel,paste0(level,"_abundance_", group, "_rel.csv"))

# write.csv(sample,paste0("sample_info_", group, ".csv"))


#提取每个样本丰度最高的菌的名字
data <- read.csv("0-data/Genus_abundance_all.csv", header = T, row.names = 1)
max_genus <- apply(data, 2, function(x) rownames(data)[which.max(x)])

result_df <- data.frame(
  Sample = names(max_genus),
  Top_Genus = as.character(max_genus))


##############################
#######PCoA###################
##############################
library(vegan)
library(ggplot2)
library(ggrepel)

rm(list = ls())
setwd("E:/Project/13-蜘蛛微生物/0-data/")
data <- read.csv(file = "./Genus_abundance_all.csv", row.names = 1, header = T, fileEncoding = "UTF-8")
sample_info <- read.csv("./sample_info.csv", header = T)

data <- t(data)   #转置
distance <- vegdist(data, method = 'bray') #计算距离
pcoa <- cmdscale(distance, eig = T) #计算主坐标
##设置要显示的标签
points <- as.data.frame(pcoa$points)
points$name <- rownames(points)

points$group <- sample_info$province
points$group <- sample_info$foraging_strategy
##设置不显示标签的分类
#points$name<-gsub(points$name,pattern = 'Y-.+',replacement = '')
colnames(points) <- c("x", "y", 'name', 'group')
eig <- pcoa$eig
sum_eig=sum(pcoa$eig)
eig_percent=round(pcoa$eig/sum_eig*100,1)

#adonis
mm = adonis2(distance~group,points,permutations = 999,method='bray')
ado = paste0("adonis R2: ",round(mm$R2,2), "; P-value: ", mm$`Pr(>F)`)

# 画图PCoA
p <- ggplot(points, aes(x, y, fill = group)) +
  stat_ellipse(aes(fill = group), geom = 'polygon', level = 0.95, alpha = 0.3, show.legend = FALSE) +  
  geom_point(alpha = 1, size = 4, shape = 21, color = "black") +
  scale_fill_viridis_d() +
  #geom_text_repel(data = points, aes(x = x, y = y), label = points$name, size = 2) + ##显示标签
  theme_bw(base_rect_size = 1) +
  theme(legend.title = element_blank(),
        #panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black", size = 0),
        axis.text = element_text(color = "black", size = 12),
        axis.title = element_text(size = 12), legend.position = 'right',
        legend.text = element_text(size = 12), legend.box.background = element_rect(color = "black"),
        axis.ticks = element_line(color = "black"), legend.background = element_blank())+
  scale_fill_manual(values = c("#D87A78FF","#61C2A2FF","#6699CCFF")) +
  #scale_fill_manual(values = c("#BEAED4FF","#7FC97FFF")) +
  labs(x=paste("PCoA 1 (", eig_percent[1], "%)", sep=""),
       y=paste("PCoA 2 (", eig_percent[2], "%)", sep=""),
       title=ado)   

# facet_wrap(~points$group)+ # 按分组拆分PCA
p
ggsave(p, file = paste("foraging_pcoa", ".pdf", sep = ""), width = 5.5, height = 4, dpi = 600, path = "E:/Project/13-蜘蛛微生物/1-result/")





###alpha多样性结果箱线图绘制###
library(ggplot2)
library(tidyr)
library(ggpubr)
library(dplyr)
library(tibble)

rm(list = ls())
setwd("E:/Project/13-蜘蛛微生物/")
alpha_data <- read.csv("./0-data/alpha_diversity.csv", header = T)

data_long <- pivot_longer(alpha_data, 
                          cols = c(simpson, shannon, faith, chao1), 
                          names_to = "diversity_index", 
                          values_to = "value")


###不同省份###
pp <- ggplot(alpha_data, aes(x = province, y = simpson, fill = province)) +
  geom_boxplot() +
  #geom_jitter()+
  stat_compare_means(method = "wilcox.test") +
  scale_fill_manual(values = c("#BEAED4FF","#7FC97FFF")) +
  labs(title = "Simpson index") +
  theme(panel.grid = element_blank())+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),  # 设置标题字体大小为 18
    axis.title.x = element_text(size = 12),  # 设置 x 轴标题字体大小为 14
    axis.title.y = element_text(size = 12),  # 设置 y 轴标题字体大小为 14
    axis.text.x = element_text(size = 12),   # 设置 x 轴刻度标签字体大小为 12
    axis.text.y = element_text(size = 12),   # 设置 y 轴刻度标签字体大小为 12
    legend.title = element_text(size = 12),  # 设置图例标题字体大小为 14
    legend.text = element_text(size = 12)    # 设置图例文本字体大小为 12
  )

pp
ggsave(pp, file = paste("province_simpson", ".pdf", sep = ""), width = 3.5, height = 3.5, dpi = 600, path = "E:/Project/13-蜘蛛微生物/1-result/FIG/")


###不同摄食方式###
alpha_data$foraging_strategy <- as.factor(alpha_data$foraging_strategy)

#不同省份的不同摄食方式
sichuan_alpha <- alpha_data[alpha_data$province=="Sichuan",]
guizhou_alpha <- alpha_data[alpha_data$province=="Guizhou",]
sichuan_alpha$foraging_strategy <- as.factor(sichuan_alpha$foraging_strategy)
guizhou_alpha$foraging_strategy <- as.factor(guizhou_alpha$foraging_strategy)

my_comparisons <- list(c("burrowing","wandering"), c("wandering", "web-building"), c("burrowing", "web-building"))
kk <- ggplot(alpha_data, aes(x = foraging_strategy, y = shannon, fill = foraging_strategy)) +
  geom_boxplot() +
  #geom_jitter()+
  stat_compare_means(label.y = 15) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") +
  scale_fill_manual(values = c("#D87A78FF","#61C2A2FF","#6699CCFF")) +
  labs(title = "Shannon index") +
  theme(panel.grid = element_blank()) +
  theme_bw() +
  theme(legend.position = "none") +
  theme(plot.title = element_text(hjust = 0.5))+
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),  # 设置标题字体大小为 18
    axis.title.x = element_text(size = 12),  # 设置 x 轴标题字体大小为 14
    axis.title.y = element_text(size = 12),  # 设置 y 轴标题字体大小为 14
    axis.text.x = element_text(size = 12, angle = 20, hjust = 1),   # 设置 x 轴刻度标签字体大小为 12
    axis.text.y = element_text(size = 12),   # 设置 y 轴刻度标签字体大小为 12
    legend.title = element_text(size = 12),  # 设置图例标题字体大小为 14
    legend.text = element_text(size = 12)    # 设置图例文本字体大小为 12
  )

kk
ggsave(kk, file = paste("foraging_shannon", ".pdf", sep = ""), width = 3.5, height = 3.5, dpi = 600, path = "E:/Project/13-蜘蛛微生物/1-result/FIG")

###单个site不同捕食策略##################
site_list <- c("JGS","JLG","LG","LM","LQS","FJS","QRG","XS")
for (site in site_list) {
  site_name <- site
  site_alpha <- alpha_data[alpha_data$sampling_site==site_name,]
  site_alpha$foraging_strategy <- as.factor(site_alpha$foraging_strategy)
  my_comparisons <- list(c("burrowing","wandering"), c("wandering", "web-building"), c("burrowing", "web-building"))
  aa <- ggplot(site_alpha, aes(x = foraging_strategy, y = shannon, fill = foraging_strategy)) +
    geom_boxplot() +
    #geom_jitter()+
    stat_compare_means(label.y = 15) +
    stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") +
    scale_fill_manual(values = c("#D87A78FF","#61C2A2FF","#6699CCFF")) +
    labs(title = paste0(site_name, " Shannon index")) +
    theme(panel.grid = element_blank()) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(aa, file = paste0("E:/Project/13-蜘蛛微生物/1-result/8-site/2-alpha/", site_name, "_shannon_index.pdf"), height = 5, width = 5)
  
  bb <- ggplot(site_alpha, aes(x = foraging_strategy, y = simpson, fill = foraging_strategy)) +
    geom_boxplot() +
    #geom_jitter()+
    stat_compare_means(label.y = 3) +
    stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") +
    scale_fill_manual(values = c("#D87A78FF","#61C2A2FF","#6699CCFF")) +
    labs(title = paste0(site_name, " Simpson index")) +
    theme(panel.grid = element_blank()) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(bb, file = paste0("E:/Project/13-蜘蛛微生物/1-result/8-site/2-alpha/", site_name, "_simpson_index.pdf"), height = 5, width = 5)
}

##############不同site##################
alpha_data$province <- as.factor(alpha_data$province)
alpha_data$sampling_site <- as.factor(alpha_data$sampling_site)


#provin内多组检验
pro_results <- alpha_data %>%
  group_by(province) %>%
  do({
    # 对每个大组执行检验
    kruskal_result <- kruskal.test(shannon ~ sampling_site, data = .)
    # 提取p值并转换为显著性符号
    p_val <- kruskal_result$p.value
    # 返回结果数据框
    data.frame(p_value = p_val)
  }) 

#两两分组p值
pairwise_results <- pairwise.wilcox.test(
  x = alpha_data$shannon,
  g = alpha_data$sampling_site,
  p.adjust.method = "bonferroni",  # p值校正方法，控制多重检验误差
  exact = FALSE              # 样本量大时建议设为FALSE
)

pairwise_pvalues <- as.data.frame(pairwise_results$p.value) %>%
  rownames_to_column("group1") %>%
  pivot_longer(
    cols = -group1,
    names_to = "group2",
    values_to = "p_value"
  ) %>%
  filter(!is.na(p_value)) %>%  # 去除NA值（重复比较和自身比较）
  # 统一组对顺序（group1 < group2，避免重复）
  mutate(
    group1 = pmin(group1, group2),
    group2 = pmax(group1, group2)
  ) %>%
  distinct() %>%  # 去重
  arrange(p_value)  # 按p值从小到大排序

write.csv(pairwise_results$p.value, "sampling_site_wilcox.test.csv")

pp <- ggplot(alpha_data, aes(x = sampling_site, y = shannon, fill = sampling_site)) +
  geom_boxplot() +
  #geom_jitter()+
  facet_grid(~ province, scales = "free_x", space='free') + 
  scale_fill_brewer(palette = "Spectral") +
  labs(title = "Shnnon index") +
  theme(panel.grid = element_blank())+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),  # 设置标题字体大小为 18
    axis.title.x = element_text(size = 14),  # 设置 x 轴标题字体大小为 14
    axis.title.y = element_text(size = 14),  # 设置 y 轴标题字体大小为 14
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1),   # 设置 x 轴刻度标签字体大小为 12
    axis.text.y = element_text(size = 14),   # 设置 y 轴刻度标签字体大小为 12
    legend.title = element_text(size = 14),  # 设置图例标题字体大小为 14
    legend.text = element_text(size = 14) # 设置图例文本字体大小为 12
  )
  
pp

ggsave(pp, filename = "E:/Project/13-蜘蛛微生物/1-result/FIG/figure4/site_shannon.pdf", height = 5, width = 7)


wa <- alpha_data[alpha_data$foraging_strategy == "burrowing",]
alpha_data <- wa
#provin内多组检验
pro_results <- alpha_data %>%
  group_by(province) %>%
  do({
    # 对每个大组执行检验
    kruskal_result <- kruskal.test(shannon ~ sampling_site, data = .)
    # 提取p值并转换为显著性符号
    p_val <- kruskal_result$p.value
    # 返回结果数据框
    data.frame(p_value = p_val)
  }) 

#两两分组p值
pairwise_results <- pairwise.wilcox.test(
  x = alpha_data$shannon,
  g = alpha_data$sampling_site,
  p.adjust.method = "bonferroni",  # p值校正方法，控制多重检验误差
  exact = FALSE              # 样本量大时建议设为FALSE
)

pairwise_pvalues <- as.data.frame(pairwise_results$p.value) %>%
  rownames_to_column("group1") %>%
  pivot_longer(
    cols = -group1,
    names_to = "group2",
    values_to = "p_value"
  ) %>%
  filter(!is.na(p_value)) %>%  # 去除NA值（重复比较和自身比较）
  # 统一组对顺序（group1 < group2，避免重复）
  mutate(
    group1 = pmin(group1, group2),
    group2 = pmax(group1, group2)
  ) %>%
  distinct() %>%  # 去重
  arrange(p_value)  # 按p值从小到大排序

write.csv(pairwise_results$p.value, "sampling_site_wilcox.test_ambushing.csv")

pp <- ggplot(wa, aes(x = sampling_site, y = shannon, fill = sampling_site)) +
  geom_boxplot() +
  #geom_jitter()+
  facet_grid(~ province, scales = "free_x", space='free') + 
  scale_fill_brewer(palette = "Spectral") +
  labs(title = "Wandering Shnnon index") +
  theme(panel.grid = element_blank())+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))

pp
ggsave(pp, filename = "site_wan_shannon.pdf", height = 5, width = 7)


########################
###相对丰度柱状堆积图###
########################
library(psych)
library(reshape2)
library(ggplot2)
library(factoextra)
library(stringr)
# BiocManager::install("wesanderson")
library(wesanderson)
library(RColorBrewer)
library(ggpubr)
library(showtext)
library('ggnested')
library(ggalluvial)


# 清空当前环境
rm(list = ls())
# 设置工作路径
setwd("E:/Project/13-蜘蛛微生物/")

# 设置需要画最高的多少个门/属
n = 10

mycol21 <- c("#A6CEE3",
             "#FFFF99", "#8c510a", "#bf812d", "#D8D155","#6A3D9A","#b30000", 
             "#7A142C","#6a51a3", "#d73027", "#fc9272", "#9ecae1",
             "#E0367A","#B2DF8A","#FB9A99","#CAB2D6",   
             "#fecc11","#FF7F00","#E31A1C",
             "#1F78B4","#33A02C")
# 生成对应的颜色映射（top20及以下使用同一个色系）
col <- c(mycol21[1:1], mycol21[(22-n):21])

level <- "Genus"
data <- read.csv(paste0("./0-data/", level, "_abundance_all.csv") , header = T, row.names = 1)
sample_info <- read.csv("./0-data/sample_info.csv", header = T)

#SC GZ不同摄食
# location <- "Sichuan"
# data <- read.csv(paste0("./0-data/", level, "_abundance_", location,".csv") , header = T, row.names = 1)
# sample_info <- read.csv(paste0("./0-data/sample_info_",location, ".csv"), header = T)


# 如果输入的是count矩阵，转为丰度，丰度矩阵则不受影响
data_per <- as.data.frame(lapply(data, function(x) x / sum(x)))
row.names(data_per) <- row.names(data) #加一下行名
# data_per <- data
# 计算每个门水平的平均丰度 便于后续筛选                                 
data.ave <- apply(data_per, 1, FUN=mean)

# 选择丰度最高的n个门 剩下的放入others里
data.2 <- cbind(data_per, data.ave)[order(-data.ave),] #排个序
data.2 <- subset(data.2[1:n,], select=-data.ave)
# 统计others丰度
data.2 <- rbind(data.2, others=apply(data.2, 2, function(x){1-sum(x)}))
# 加一列行名 便于后续的长宽转换
data.2 <- cbind(dataID=row.names(data.2), data.2)

# 长宽转换
# 因子排个序
data.2$dataID <- factor(data.2$dataID, levels = rev(data.2$dataID))
data.gg <- melt(data.2, id.vars="dataID", variable.name="sample_name", value.name="Abundance")

# 添加分组信息
data.gg <- data.gg %>%
  merge(sample_info, by = "sample_name")
data.gg$foraging_strategy <- as.factor(data.gg$foraging_strategy)
data.gg$province <- as.factor(data.gg$province)
data.gg$enterotype <- as.factor(data.gg$enterotype)

# 分组柱图需要
data.foraging <- aggregate(Abundance ~ dataID + foraging_strategy, data = data.gg, FUN = mean)
data.province <- aggregate(Abundance ~ dataID + province, data = data.gg, FUN = mean)
data.enterotype <- aggregate(Abundance ~ dataID + enterotype, data = data.gg, FUN = mean)


# 绘图

# 按分组绘图data.province或data.foraging或data.enterotype 
ggplot(data.province, aes(x = province, y = Abundance, fill = dataID, stratum=dataID, alluvium=dataID)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7, size = 0.25) +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) + # 去掉坐标轴多余部分
  # scale_color_manual(values = col)+
  # facet_grid(~ enterotype) + 
  labs(x = "", y = "Relative Abundance") +
  # 设置图例显示的列数
  guides(fill = guide_legend(ncol = 1,keywidth = 0.9,keyheight = 0.9,title = level)) +
  theme_bw() +
  # 添加柱子之间的连线
  geom_alluvium(alpha = 0.4,
                width = 0.7,
                curve_type = "linear")+
  scale_fill_manual(values=col, labels = data.province$dataID%>%str_wrap(width = 36)) + # 颜色设置
  theme(axis.text.x = element_text(size = 12, vjust = 1.05, hjust = 1, angle = 45),
        # 设置固定宽高比例
        aspect.ratio = 1.5,
        # panel.grid =element_blank(),
        axis.title.y = element_text(size = 15),
        axis.ticks.x = element_blank(),
        axis.title = element_text(),
        axis.text = element_text(size = 15),
        plot.title = element_text(size = 15),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 15))

ggsave(filename = paste0("E:/Project/13-蜘蛛微生物/1-result/FIG/figure2/province_group_barplot_top",n,"_", level, ".pdf"), device="pdf", width=6, height=6)
ggsave(filename = paste0("E:/Project/13-蜘蛛微生物/1-result/3-SC&GZ不同摄食/3_bar/", location, "_group_barplot_top",n,"_", level, ".pdf"), device="pdf", width=10, height=10)

ggplot(data.foraging, aes(x = foraging_strategy, y = Abundance, fill = dataID, stratum=dataID, alluvium=dataID)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7, size = 0.25) +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) + # 去掉坐标轴多余部分
  # scale_color_manual(values = col)+
  # facet_grid(~ group, scales = "free_x", space='free') + 
  labs(x = "", y = "Relative Abundance") +
  # 设置图例显示的列数
  guides(fill = guide_legend(ncol = 1,keywidth = 0.9,keyheight = 0.9,title = level)) +
  theme_bw() +
  # 添加柱子之间的连线
  geom_alluvium(alpha = 0.4,
                width = 0.7,
                curve_type = "linear")+
  scale_fill_manual(values=col, labels = data.foraging$dataID%>%str_wrap(width = 36)) + # 颜色设置
  theme(axis.text.x = element_text(size = 12, vjust = 1.05, hjust = 1, angle = 45),
        # 设置固定宽高比例
        aspect.ratio = 1.7,
        # panel.grid =element_blank(),
        axis.title.y = element_text(size = 15),
        axis.ticks.x = element_blank(),
        axis.title = element_text(),
        axis.text = element_text(size = 15),
        plot.title = element_text(size = 15),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 15))

ggsave(filename = paste0("E:/Project/13-蜘蛛微生物/1-result/FIG/foraging_group_barplot_top",n,"_", level, "_", location, ".pdf"), device="pdf", width=6, height=6)

ggsave(filename = paste0("E:/Project/13-蜘蛛微生物/1-result/FIG/foraging_group_barplot_top",n,"_", level, ".pdf"), device="pdf", width=6, height=6)

ggplot(data.enterotype, aes(x = enterotype, y = Abundance, fill = dataID, stratum=dataID, alluvium=dataID)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7, size = 0.25) +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) + # 去掉坐标轴多余部分
  # scale_color_manual(values = col)+
  # facet_grid(~ group, scales = "free_x", space='free') + 
  labs(x = "", y = "Relative Abundance") +
  # 设置图例显示的列数
  guides(fill = guide_legend(ncol = 1,keywidth = 0.9,keyheight = 0.9,title = level)) +
  theme_bw() +
  # 添加柱子之间的连线
  geom_alluvium(alpha = 0.4,
                width = 0.7,
                curve_type = "linear")+
  scale_fill_manual(values=col, labels = data.foraging$dataID%>%str_wrap(width = 36)) + # 颜色设置
  theme(axis.text.x = element_text(size = 12, vjust = 1.05, hjust = 1, angle = 45),
        # 设置固定宽高比例
        aspect.ratio = 1.7,
        # panel.grid =element_blank(),
        axis.title.y = element_text(size = 15),
        axis.ticks.x = element_blank(),
        axis.title = element_text(),
        axis.text = element_text(size = 15),
        plot.title = element_text(size = 15),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 15))
ggsave(filename = paste0("E:/Project/13-蜘蛛微生物/1-result/FIG/figure5//enterotype_group_barplot_top",n,"_", level, ".pdf"), device="pdf", width=6, height=7)


# 按样品绘图
ggbarplot(data.gg, x = "sample_name", y="Abundance", fill="dataID",color = '#373535',
          legend="right", 
          legend.title="Top Genus" # main="Relative abundance per Phylum",
          # font.main = c(14,"bold", "black"), font.x = c(12, "bold"), font.y=c(12,"bold")
) + 
  scale_y_continuous(expand = expansion(mult = c(0,0)))+ # 去掉坐标轴多余部分
  # rotate_x_text() + 
  # scale_fill_manual(values=brewer.pal(11, 'Paired')) + # 颜色设置
  scale_fill_manual(values=col) +
  facet_grid(~ foraging_strategy, scales = "free_x", space='free') + 
  labs(x = "", y = "") + 
  # 设置图例显示的列数
  guides(fill = guide_legend(ncol = 1))+
  theme_bw() +
  # theme(axis.text.x=element_text(size = 6, family = "Arial", face = "bold", vjust = 1.2, hjust = 1, angle = 45)),
  theme(axis.text.x=element_text(size = 6, family = "Arial", vjust = 1.2, hjust = 1, angle = 45),
        # panel.grid =element_blank(),
        axis.ticks.x=element_blank(),
        axis.title = element_text(), 
        plot.title = element_text(), 
        legend.title = element_text()) 

ggsave(filename = paste0("samples_barplot_top",n,"_G.pdf"), device="pdf", width=49, height=8)


#############################
#####Lefse 差异分析##########
#############################
library(tidyverse)
library(microeco)
library(magrittr)

setwd("E:/Project/13-蜘蛛微生物/")
data <- read.csv("0-data/Genus_abundance_all.csv", header = T, row.names = 1)
sample_info <- read.csv("0-data/sample_info.csv", header = T)
otu_table <- data %>%
  mutate(otu_id = row_number())
rownames(otu_table) <- otu_table$otu_id
otu_table <- otu_table[, !names(otu_table) %in% "otu_id", drop = FALSE]
taxa_table <- as.data.frame(rownames(data))
colnames(taxa_table) <- "Species"
rownames(sample_info) <- sample_info$sample_name

dataset <- microtable$new(otu_table = otu_table,
                          sample_table = sample_info,
                          tax_table = taxa_table)
dataset$tidy_dataset()

lefse <- trans_diff$new(dataset = dataset, method = "lefse", group = "province", alpha = 0.01, lefse_subgroup = NULL)

# we show 30 taxa with the highest LDA (log10)
lefse$plot_diff_bar(use_number = 1:30, 
                    width = 0.8, 
                    group_order = c("burrowing", "wandering", "web-building"),
                    color_values = c("#D87A78FF","#61C2A2FF","#6699CCFF"),
                    xtext_size = 12, ytext_size = 12, axis_text_y = 12,
                    keep_prefix = FALSE
                    )
ggsave("E:/Project/13-蜘蛛微生物/1-result/FIG/lefse_species.pdf", height = 6, width = 6)
write.csv(lefse$res_diff %>% subset(LDA > 2), "E:/Project/13-蜘蛛微生物/1-result/2-不同摄食方式/4-lefse/lefse_phylum.csv")

lefse$plot_diff_bar(use_number = 1:30, 
                    width = 0.8, 
                    group_order = c("Guizhou","Sichuan"),
                    color_values = c("#BEAED4FF","#7FC97FFF"),
                    xtext_size = 12, ytext_size = 12, axis_text_y = 12,
                    keep_prefix = FALSE
                    )
ggsave("E:/Project/13-蜘蛛微生物/1-result/FIG/figure2/lefse_province_genus.pdf", height = 6, width = 6)

#######SC GZ不同摄食比较#########
setwd("E:/Project/13-蜘蛛微生物/")
location <- "Sichuan"
level <- "Genus"
data <- read.csv(paste0("0-data/", level, "_abundance_", location,".csv"), header = T, row.names = 1)
sample_info <- read.csv(paste0("0-data/sample_info_", location, ".csv"), header = T)
otu_table <- data %>%
  mutate(otu_id = row_number())
rownames(otu_table) <- otu_table$otu_id
otu_table <- otu_table[, !names(otu_table) %in% "otu_id", drop = FALSE]
taxa_table <- as.data.frame(rownames(data))
colnames(taxa_table) <- "Species"
rownames(sample_info) <- sample_info$sample_name

dataset <- microtable$new(otu_table = otu_table,
                          sample_table = sample_info,
                          tax_table = taxa_table)
dataset$tidy_dataset()

lefse <- trans_diff$new(dataset = dataset, method = "lefse", group = "foraging_strategy", alpha = 0.01, lefse_subgroup = NULL)

# we show 20 taxa with the highest LDA (log10)
lefse$plot_diff_bar(use_number = 1:30, 
                    width = 0.8, 
                    group_order = c("burrowing", "wandering", "web-building"),
                    color_values = c("#D87A78FF","#61C2A2FF","#6699CCFF"),
                    xtext_size = 12, ytext_size = 12, axis_text_y = 12,
                    keep_prefix = FALSE,
                    coord_flip = FALSE)

ggsave(paste0("E:/Project/13-蜘蛛微生物/1-result/FIG/figure2/", location, "_lefse_", level, ".pdf"), height = 4, width = 8)
write.csv(lefse$res_diff %>% subset(LDA > 2), paste0("E:/Project/13-蜘蛛微生物/1-result/3-SC&GZ不同摄食/4-lefse/", location, "_lefse_", level, ".csv"))


###############################
########菌群分型数目统计#######
###############################
#已添加样本肠道菌群分型到样本信息表

setwd("E:/Project/13-蜘蛛微生物/")
data <- read.csv("0-data/sample_info.csv", header = T)

province_count <- data %>%
  group_by(enterotype, province) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(category = "province")

foraging_strategy_count <- data %>%
  group_by(enterotype, foraging_strategy) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(category = "foraging_strategy")

family_count <- data %>%
  group_by(enterotype, Family) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(category = "family")


province_plot <- ggplot(province_count, aes(x = enterotype, y = count, fill = province)) +
  geom_col(position = "stack") +
  scale_fill_manual(values = c("#BEAED4FF","#7FC97FFF")) +
  geom_text(aes(label = count), position = position_stack(vjust = 0.5)) +
  labs(
    x = "Enterotype",
    y = "Count",
    fill = "Region"
  ) +
  theme_minimal()+
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),  # 设置标题字体大小为 18
    axis.title.x = element_text(size = 14),  # 设置 x 轴标题字体大小为 14
    axis.title.y = element_text(size = 14),  # 设置 y 轴标题字体大小为 14
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1),   # 设置 x 轴刻度标签字体大小为 12
    axis.text.y = element_text(size = 14),   # 设置 y 轴刻度标签字体大小为 12
    legend.title = element_text(size = 14),  # 设置图例标题字体大小为 14
    legend.text = element_text(size = 14)    # 设置图例文本字体大小为 12
  )

province_plot
ggsave(province_plot, filename = "E:/Project/13-蜘蛛微生物/1-result/FIG/figure5/province_enterotype.pdf", height = 5, width = 5)

foraging_plot <- ggplot(foraging_strategy_count, aes(x = enterotype, y = count, fill = foraging_strategy)) +
  geom_col(position = "stack") +
  scale_fill_manual(values = c("#D87A78FF","#61C2A2FF","#6699CCFF")) +
  geom_text(aes(label = count), position = position_stack(vjust = 0.5)) +
  labs(
    x = "Enterotype",
    y = "Count",
    fill = "Foraging strategy"
  ) +
  theme_minimal()+
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),  # 设置标题字体大小为 18
    axis.title.x = element_text(size = 14),  # 设置 x 轴标题字体大小为 14
    axis.title.y = element_text(size = 14),  # 设置 y 轴标题字体大小为 14
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1),   # 设置 x 轴刻度标签字体大小为 12
    axis.text.y = element_text(size = 14),   # 设置 y 轴刻度标签字体大小为 12
    legend.title = element_text(size = 14),  # 设置图例标题字体大小为 14
    legend.text = element_text(size = 14)    # 设置图例文本字体大小为 12
  )
  

foraging_plot
ggsave(foraging_plot, filename = "E:/Project/13-蜘蛛微生物/1-result/FIG/figure5/foraging_enterotype.pdf", height = 5, width = 5)

family_plot <- ggplot(family_count, aes(x = enterotype, y = count, fill = Family)) +
  geom_col(position = "stack") +
  scale_fill_viridis_d() +  #适合颜色多的，自动生成
  geom_text(aes(label = count), position = position_stack(vjust = 0.5)) +
  labs(
    x = "Enterotype",
    y = "Count",
    fill = "Host Family"
  ) +
  theme_minimal()

write.csv(family_count, "E:/Project/13-蜘蛛微生物/1-result/4-肠型分析/hostfamily_enterotype_count.csv")
ggsave(family_plot, filename = "E:/Project/13-蜘蛛微生物/1-result/4-肠型分析/hostfamily_enterotype.pdf", height = 20, width = 5)

###############################
########方差分解VPA#######
###############################
library(vegan)
library(dplyr)

setwd("E:/Project/13-蜘蛛微生物/")
otu <- read.csv("0-data/Genus_abundance_all.csv", header = T, row.names = 1)  #丰度矩阵
env_data <- read.csv("0-data/sample_info.csv", header = T)  #分组信息表

otu <- t(otu)

# 单变量RDA分析
rda_province <- rda(otu ~ province, data = env_data)
rda_foraging_strategy <- rda(otu ~ foraging_strategy, data = env_data)
rda_sampling_site <- rda(otu ~ sampling_site, data = env_data)
rda_family <- rda(otu ~ Family, data = env_data)
rda_genus <- rda(otu ~ Genus, data = env_data)
rda_species <- rda(otu ~ Species, data = env_data)
rda_sex <- rda(otu ~ sex, data = env_data)
rda_altitude <- rda(otu ~ altitude, data = env_data)
rda_season <- rda(otu ~ season, data = env_data)

# 输出RDA分析结果
summary(rda_province)
summary(rda_foraging_strategy)
summary(rda_sampling_site)
summary(rda_family)
summary(rda_genus)
summary(rda_species)
summary(rda_sex)
summary(rda_altitude)
summary(rda_season)




# 显著性检验（基于置换）
anova_rda_province <- anova(rda_province, permutations = 999)
anova_rda_foraging_strategy <- anova(rda_province, permutations = 999)
anova_rda_sampling_site <- anova(rda_sampling_site, permutations = 999)
anova_rda_family <- anova(rda_family, permutations = 999)
anova_rda_genus <- anova(rda_genus, permutations = 999)
anova_rda_species <- anova(rda_species, permutations = 999)
anova_rda_sex <- anova(rda_sex, permutations = 999)
anova_rda_altitude <- anova(rda_altitude, permutations = 999)
anova_rda_season <- anova(rda_season, permutations = 999)

# 输出显著性检验结果
print(anova_rda_province)
print(anova_rda_foraging_strategy)
print(anova_rda_family)
print(anova_rda_genus)
print(anova_rda_sex)

# 创建解释变量表
expl_host <- data.frame(
  family = env_data$Family,
  genus = env_data$Genus, 
  species = env_data$Species
)
expl_province <- data.frame(province = env_data$province)
expl_foraging_strategy <- data.frame(foraging_strategy = env_data$foraging_strategy)
expl_sex <- data.frame(sex = env_data$sex)
expl_altitude <- data.frame(altitude = env_data$altitude)
expl_sampling_site <- data.frame(sampling_site = env_data$sampling_site)
expl_season <- data.frame(season = env_data$season)

# 两两组合
p1 <- varpart(otu, expl_host, expl_province, transfo = "hel")
p2 <- varpart(otu, expl_host, expl_foraging_strategy, transfo = "hel")
p3 <- varpart(otu, expl_host, expl_sex, transfo = "hel")
p4 <- varpart(otu, expl_host, expl_altitude, transfo = "hel")
p5 <- varpart(otu, expl_host, expl_sampling_site, transfo = "hel")
p6 <- varpart(otu, expl_host, expl_season, transfo = "hel")
p7 <- varpart(otu, expl_province, expl_foraging_strategy, transfo = "hel")
p8 <- varpart(otu, expl_province, expl_sex, transfo = "hel")
p9 <- varpart(otu, expl_province, expl_altitude, transfo = "hel")
p10 <- varpart(otu, expl_province, expl_sampling_site, transfo = "hel")
p11 <- varpart(otu, expl_province, expl_season, transfo = "hel")
p12 <- varpart(otu, expl_foraging_strategy, expl_sex, transfo = "hel")
p13 <- varpart(otu, expl_foraging_strategy, expl_altitude, transfo = "hel")
p14 <- varpart(otu, expl_foraging_strategy, expl_sampling_site, transfo = "hel")
p15 <- varpart(otu, expl_foraging_strategy, expl_season, transfo = "hel")
p16 <- varpart(otu, expl6, expl7, transfo = "hel")
p17 <- varpart(otu, expl6, expl8, transfo = "hel")
p18 <- varpart(otu, expl6, expl9, transfo = "hel")
p19 <- varpart(otu, expl7, expl8, transfo = "hel")
p20 <- varpart(otu, expl7, expl9, transfo = "hel")
p21 <- varpart(otu, expl8, expl9, transfo = "hel")


data1111 <- data.frame(Combination = c("p1: lifestyle + site", "p2: lifestyle + family"), 
                   Contribution = c(summary(p1)$part[1, "Unique"], summary(p14$part[1, "Unique"]))) 



##################################
########功能绘图##################
##################################

#功能为注释的每个样本的相对丰度表
rm(list = ls())
# 设置工作路径
setwd("E:/Project/13-蜘蛛微生物/")

# 设置需要画最高的多少个功能
n = 10

mycol21 <- c("#A6CEE3",
             "#FFFF99", "#8c510a", "#bf812d", "#D8D155","#6A3D9A","#b30000", 
             "#7A142C","#6a51a3", "#d73027", "#fc9272", "#9ecae1",
             "#E0367A","#B2DF8A","#FB9A99","#CAB2D6",   
             "#fecc11","#FF7F00","#E31A1C",
             "#1F78B4","#33A02C")
# 生成对应的颜色映射（top20及以下使用同一个色系）
col <- c(mycol21[1:1], mycol21[(22-n):21])


data <- read.csv("./0-data/function_rel_abundance.csv" , header = T, row.names = 1)
sample_info <- read.csv("./0-data/sample_info.csv", header = T)

# 计算每个功能平均丰度 便于后续筛选                                 
data.ave <- apply(data, 1, FUN=mean)

# 选择丰度最高的n个门 剩下的放入others里
data.2 <- cbind(data, data.ave)[order(-data.ave),] #排个序
data.2 <- subset(data.2[1:n,], select=-data.ave)
# 统计others丰度
data.2 <- rbind(data.2, others=apply(data.2, 2, function(x){1-sum(x)}))
# 加一列行名 便于后续的长宽转换
data.2 <- cbind(dataID=row.names(data.2), data.2)

# 长宽转换
# 因子排个序
data.2$dataID <- factor(data.2$dataID, levels = rev(data.2$dataID))
data.gg <- melt(data.2, id.vars="dataID", variable.name="sample_name", value.name="Abundance")

# 添加分组信息
data.gg <- data.gg %>%
  merge(sample_info, by = "sample_name")
data.gg$foraging_strategy <- as.factor(data.gg$foraging_strategy)
data.gg$province <- as.factor(data.gg$province)

# 分组柱图需要
data.foraging <- aggregate(Abundance ~ dataID + foraging_strategy, data = data.gg, FUN = mean)
data.province <- aggregate(Abundance ~ dataID + province, data = data.gg, FUN = mean)


# 绘图

# 按分组绘图data.province或data.foraging或data.enterotype 
ggplot(data.province, aes(x = province, y = Abundance, fill = dataID, stratum=dataID, alluvium=dataID)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7, size = 0.25) +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) + # 去掉坐标轴多余部分
  labs(x = "", y = "Relative Abundance") +
  # 设置图例显示的列数
  guides(fill = guide_legend(ncol = 1,keywidth = 0.9,keyheight = 0.9)) +
  theme_bw() +
  # 添加柱子之间的连线
  geom_alluvium(alpha = 0.4,
                width = 0.7,
                curve_type = "linear")+
  scale_fill_manual(values=col, labels = data.province$dataID%>%str_wrap(width = 36)) + # 颜色设置
  theme(axis.text.x = element_text(size = 12, vjust = 1.05, hjust = 1, angle = 45),
        # 设置固定宽高比例
        aspect.ratio = 1.2,
        # panel.grid =element_blank(),
        axis.title.y = element_text(size = 15),
        axis.ticks.x = element_blank(),
        axis.title = element_text(),
        axis.text = element_text(size = 15),
        plot.title = element_text(size = 15),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 15))

ggsave(filename = paste0("E:/Project/13-蜘蛛微生物/1-result/7-function/province_function_barplot_top",n,".pdf"), device="pdf", width=10, height=10)


ggplot(data.foraging, aes(x = foraging_strategy, y = Abundance, fill = dataID, stratum=dataID, alluvium=dataID)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7, size = 0.25) +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) + # 去掉坐标轴多余部分
  labs(x = "", y = "Relative Abundance") +
  # 设置图例显示的列数
  guides(fill = guide_legend(ncol = 1,keywidth = 0.9,keyheight = 0.9)) +
  theme_bw() +
  # 添加柱子之间的连线
  geom_alluvium(alpha = 0.4,
                width = 0.7,
                curve_type = "linear")+
  scale_fill_manual(values=col, labels = data.foraging$dataID%>%str_wrap(width = 36)) + # 颜色设置
  theme(axis.text.x = element_text(size = 12, vjust = 1.05, hjust = 1, angle = 45),
        # 设置固定宽高比例
        aspect.ratio = 1.2,
        # panel.grid =element_blank(),
        axis.title.y = element_text(size = 15),
        axis.ticks.x = element_blank(),
        axis.title = element_text(),
        axis.text = element_text(size = 15),
        plot.title = element_text(size = 15),
        legend.text = element_text(size = 11),
        legend.title = element_text(size = 15))

ggsave(filename = paste0("E:/Project/13-蜘蛛微生物/1-result/7-function/foraging_function_barplot_top",n,".pdf"), device="pdf", width=10, height=10)


#########################
#######群落组装##########
#########################
# wd="/home/user001/data/zhizhu/raw/"
# com.file <- "table_raw.csv"
# tree.file <- "tree_new.nwk"
# treat.file <- "treat2col.csv"
# env.file <- "environment.csv"
# save.wd="/home/user001/data/zhizhu/raw/result/"
# prefix <- "ZZ"
# nworker <- 20
# 
# 
# library(iCAMP)
# library(ape)
# setwd(wd)
# comm=t(read.table(com.file, header = TRUE, sep = ",", row.names = 1,
#                   as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
#                   check.names = FALSE))
# tree=read.tree(file = tree.file)
# treat=read.table(treat.file, header = TRUE, sep = ",", row.names = 1,
#                  as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
#                  check.names = FALSE)
# env=read.table(env.file, header = TRUE, sep = ",", row.names = 1,
#                as.is = TRUE, stringsAsFactors = FALSE, comment.char = "",
#                check.names = FALSE) # skip this if you do not have env.file
# pd.wd <- paste0("/home/user001/data/zhizhu/raw")
# 
# setwd(save.wd)
# if (!file.exists("pd.desc"))
# {
#   pd.big = iCAMP::pdist.big(tree = tree, wd = save.wd, nworker = nworker)
# } else {
#   # 如果您已经在之前的运行中计算了系统发育距离矩阵
#   print(123)
#   pd.big = list()
#   pd.big$tip.label = read.csv(paste0(save.wd, "/pd.taxon.name.csv"), row.names = 1, stringsAsFactors = FALSE)[, 1]
#   pd.big$pd.wd = save.wd
#   pd.big$pd.file = "pd.desc"
#   pd.big$pd.name.file = "pd.taxon.name.csv"
# }
# 
# # 设置线程，不然服务器会炸，某个过程会调用全部线程
# # 设置 OpenBLAS 线程数
# Sys.setenv(OPENBLAS_NUM_THREADS = 1)
# # 设置 MKL 线程数
# Sys.setenv(MKL_NUM_THREADS = 1)
# 
# #icamp
# icres = iCAMP::icamp.big(comm = comm, tree = tree, rand = 1000, nworker = nworker, prefix = prefix,
#                          pd.desc = pd.big$pd.file, pd.spname = pd.big$tip.label,
#                          pd.wd = pd.big$pd.wd, bin.size.limit = 24, ds = 0.2,
#                          sig.index = "SES.RC")
# save(icres, file = paste0(prefix, ".icamp.SES.RC.rda"))
# 
# icbin=iCAMP::icamp.bins(icamp.detail = icres$detail,treat = treat,
#                         silent=FALSE, boot = TRUE,
#                         rand.time = 1000,between.group = TRUE)
# save(icbin,file = paste0(prefix,".iCAMP.Summary.rda")) 
# write.csv(icbin$Pt,file = paste0(prefix,".ProcessImportance_EachGroup.csv"),row.names = FALSE)
# write.csv(icbin$Ptk,file = paste0(prefix,".ProcessImportance_EachBin_EachGroup.csv"),row.names = FALSE)
# write.csv(icbin$Ptuv,file = paste0(prefix,".ProcessImportance_EachTurnover.csv"),row.names = FALSE)
# write.csv(icbin$BPtk,file = paste0(prefix,".BinContributeToProcess_EachGroup.csv"),row.names = FALSE)
# write.csv(data.frame(ID=rownames(icbin$Class.Bin),icbin$Class.Bin,stringsAsFactors = FALSE),
#           file = paste0(prefix,".Taxon_Bin.csv"),row.names = FALSE)
# write.csv(icbin$Bin.TopClass,file = paste0(prefix,".Bin_TopTaxon.csv"),row.names = FALSE)
# 
# 
# for (i in 1:3){
#   treat.use=treat[,i,drop=FALSE]
#   icamp.result=icres$CbMPDiCBraya
#   icboot=iCAMP::icamp.boot(icamp.result = icamp.result,treat = treat.use,rand.time = rand.time,
#                            compare = TRUE,silent = FALSE,between.group = TRUE,ST.estimation = TRUE)
#   save(icboot,file=paste0(prefix,".iCAMP.Boot.",colnames(treat)[i],".rda"))
#   write.csv(icboot$summary,file = paste0(prefix,".iCAMP.BootSummary.",colnames(treat)[i],".csv"),row.names = FALSE)
#   write.csv(icboot$compare,file = paste0(prefix,".iCAMP.Compare.",colnames(treat)[i],".csv"),row.names = FALSE)
# }
# 


##绘图
library(ggplot2)
library(tidyr)
library(ggrepel)
setwd("E:/Project/13-蜘蛛微生物/1-result/9-群落组装")
data <- read.csv("ZZ.ProcessImportance_EachGroup.csv", header = T)
foraging_strategy <- data[1:3,3:8]
long_foraging <- pivot_longer(foraging_strategy, cols = -Group, names_to = "Variable", values_to = "Value")
province <- data[10:11,3:8]
long_province <- pivot_longer(province, cols = -Group, names_to = "Variable", values_to = "Value")

custom_colors <- c("HeS" = "#C8A2C8", "HoS" = "#F7CAC9", "DL" = "#CCCCFF", "HD" = "#A020F0", "DR" = "#5A4FCF")

pp <- ggplot(long_foraging, aes(x="", y = Value, fill = Variable)) +
  geom_col() +
  geom_text(
    aes(label = paste0(round(Value * 100, 2), "%")),
    position = position_stack(vjust = 0.6),
    show.legend = FALSE,
    # segment.color = "black",
    # segment.size = 1,
    # min.segment.length = 0.1, # 设置线条最小长度
    # segment.curvature = -0.1, # 调整线条弯曲度
    # segment.ncp = 3, # 线条控制点数量
    # segment.angle = 20, # 线条角度
    size = 5
  ) +
  coord_polar("y", start = 0) + # 确保所有线条都显示
  theme_void() +
  scale_fill_manual(values = custom_colors) +
  facet_wrap(~Group, ncol = 3) +  # 按照Group变量分面绘制饼图,ncol设置3为3列
  theme(
    plot.title = element_text(size = 14),  # 设置标题文字大小
    legend.text = element_text(size = 14), # 设置图例文字大小
    strip.text = element_text(size = 14)   # 设置分面标签文字大小
  )
pp
ggsave(pp, filename = "E:/Project/13-蜘蛛微生物/1-result/FIG/figure5/Process_foraging_strategy.pdf", height = 4, width = 8)

kk <- ggplot(long_province, aes(x="", y = Value, fill = Variable)) +
  geom_col() +
  geom_text_repel(
    aes(label = paste0(round(Value * 100, 2), "%")),
    position = position_stack(vjust = 0.6),
    show.legend = FALSE,
    # segment.color = "black",
    # segment.size = 1,
    # min.segment.length = 0.1, # 设置线条最小长度
    # segment.curvature = -0.1, # 调整线条弯曲度
    # segment.ncp = 3, # 线条控制点数量
    # segment.angle = 20, # 线条角度
    size = 5
  ) +
  coord_polar("y", start = 0) + # 确保所有线条都显示
  theme_void() +
  scale_fill_manual(values = custom_colors) +
  facet_wrap(~Group, ncol = 2) +  # 按照Group变量分面绘制饼图,ncol设置3为3列
  theme(
    plot.title = element_text(size = 14),  # 设置标题文字大小
    legend.text = element_text(size = 14), # 设置图例文字大小
    strip.text = element_text(size = 14)   # 设置分面标签文字大小
  )

kk
ggsave(kk, filename = "E:/Project/13-蜘蛛微生物/1-result/FIG/figure5/Process_province.pdf", height = 3.5, width = 6)


##################################
########内共生菌统计绘图##########
##################################
library(dplyr)
library(ggplot2)
#province
data <- read.table("E:/Project/13-蜘蛛微生物/0-data/endosymbiont_province.txt", header = T)

# 将 genus 列转换为有序因子，顺序为数据中出现的顺序
data$genus <- factor(data$genus, levels = unique(data$genus))

kk <- ggplot(data, aes(x = genus, y = percent, fill = province)) +
  geom_col(position = "dodge") +
  labs(title = "Prevalence of endosymbiont",
       x = "Endosymbiont",
       y = "Percent (%)") +
  scale_fill_manual(values = c("SC" = "#7FC97FFF", "GZ" = "#BEAED4FF")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
 
ggsave(kk, filename = "province.pdf", height = 4, width = 6)

#foraging_strategy
data <- read.table("E:/Project/13-蜘蛛微生物/0-data/endosymbiont_for.txt", header = T)

# 将 genus 列转换为有序因子，顺序为数据中出现的顺序
data$genus <- factor(data$genus, levels = unique(data$genus))

kk <- ggplot(data, aes(x = genus, y = percent, fill = foraging_strategy)) +
  geom_col(position = "dodge") +
  labs(title = "Prevalence of endosymbiont",
       x = "Endosymbiont",
       y = "Percentage (%)") +
  scale_fill_manual(values = c("#D87A78FF","#61C2A2FF","#6699CCFF")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

kk

ggsave(kk, filename = "foraging_strategy.pdf", height = 4, width = 6)


##################################
############环境数据##############
##################################

#气温线
library(tidyr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
setwd("E:/Project/13-蜘蛛微生物/数据/环境/")
temp <- read.csv("Temp.csv", header = T)

data_long <- pivot_longer(temp, cols = starts_with("Temp"), names_to = "Month", values_to = "Temperature")
data_long$Month <- as.numeric(gsub("Temp", "", data_long$Month))

#原始数据单位为0.1°C，转为1°C
data_long$Temperature <- data_long$Temperature / 10

pp <- ggplot(data_long, aes(x = Month, y = Temperature, group = site, color = site, linetype = province)) +
  geom_line(linewidth = 1.2) +
  scale_linetype_manual(values = c("SC" = "solid", "GZ" = "twodash")) +
  scale_color_brewer(palette = "Spectral") +
  labs(x = "Month",
       y = paste("Temperature (", "\u00B0", "C)"),
       color = "site",
       linetype = "Region") +
  theme_minimal()+
  theme(
    axis.title = element_text(size = 14),  # 修改坐标轴标题字体大小
    axis.text = element_text(size = 12),   # 修改坐标轴刻度字体大小
    plot.title = element_text(size = 16),  # 修改图标题字体大小（如果有）
    legend.text = element_text(size = 14), # 修改图例文字字体大小
    legend.title = element_text(size = 14) # 修改图例标题字体大小
  )
pp
ggsave(pp, filename = "E:/Project/13-蜘蛛微生物/1-result/FIG/figure4/Temp.pdf", height = 4.5, width = 7)


#降水
#五年
rain <- read.csv("rain.csv", header = T)
data_long <- pivot_longer(rain, cols = starts_with("Month"), names_to = "Month", values_to = "Precipitation")
data_long$Month <- as.numeric(gsub("Month", "", data_long$Month))
#原始数据单位为0.1mm，转为1mm
data_long$Precipitation <- data_long$Precipitation / 10

pp <- ggplot(data_long, aes(x = Month, y = Precipitation, group = interaction(site, year), color = site, linetype = factor(province))) +
  geom_line(size = 1.2) +
  scale_color_brewer(palette = "Spectral") + 
  scale_linetype_manual(values = c("SC" = "solid", "GZ" = "twodash")) +
  labs(x = "Month",
       y = "Amount of precipitation (mm)",
       color = "site",
       linetype = "Province") +
  theme_minimal()+
  theme(
    axis.title = element_text(size = 14),  # 修改坐标轴标题字体大小
    axis.text = element_text(size = 12),   # 修改坐标轴刻度字体大小
    plot.title = element_text(size = 16),  # 修改图标题字体大小（如果有）
    legend.text = element_text(size = 14), # 修改图例文字字体大小
    legend.title = element_text(size = 14) # 修改图例标题字体大小
  )

pp
ggsave(pp, filename = "E:/Project/13-蜘蛛微生物/1-result/FIG/figure4/Rain_year.pdf", height = 4.5, width = 7)

#平均降雨
average_data <- data_long %>%
  group_by(province, site, Month) %>%
  summarize(Average_Precipitation = mean(Precipitation, na.rm = TRUE)) %>%
  ungroup()

kk <- ggplot(average_data, aes(x = Month, y = Average_Precipitation, group = site, color = site, linetype = factor(province))) +
  geom_line(size = 1.2) +
  scale_color_brewer(palette = "Spectral") + 
  scale_linetype_manual(values = c("SC" = "solid", "GZ" = "twodash")) +
  labs(x = "Month",
       y = "Amount of precipitation (mm)",
       color = "site",
       linetype = "Province") +
  theme_minimal()+
  theme(
    axis.title = element_text(size = 14),  # 修改坐标轴标题字体大小
    axis.text = element_text(size = 12),   # 修改坐标轴刻度字体大小
    plot.title = element_text(size = 16),  # 修改图标题字体大小（如果有）
    legend.text = element_text(size = 14), # 修改图例文字字体大小
    legend.title = element_text(size = 14) # 修改图例标题字体大小
  )

kk
ggsave(kk, filename = "E:/Project/13-蜘蛛微生物/1-result/FIG/figure4/Rain_avg.pdf", height = 4.5, width = 7)


#19bio因子
library(tidyverse)
library(ggpubr)
library(gridExtra)
library(grid)

bio <- read.csv("bio19.csv", header = T)
# 将数据转换为长格式
data_long <- bio %>%
  pivot_longer(cols = starts_with("bio"), names_to = "bio_factor", values_to = "value")

# 定义一个函数来绘制箱线图并进行差异检验
plot_and_test <- function(factor_name) {
  subset_data <- data_long %>% filter(bio_factor == factor_name)
  
  p <- ggplot(subset_data, aes(x = province, y = value, fill = province)) +
    geom_boxplot() +
    scale_fill_manual(values = c("SC" = "#7FC97FFF", "GZ" = "#BEAED4FF")) +
    labs(title = paste("Boxplot of", factor_name),
         x = "Region",
         y = factor_name) +
    theme_minimal() +
    theme(legend.position = "none")
  
  # 进行Wilcoxon秩和检验
  test_result <- wilcox.test(value ~ province, data = subset_data)
  
  # 在图上添加检验结果
  p <- p + annotate("text", x = 1.5, y = max(subset_data$value), 
                    label = paste("p-value =", round(test_result$p.value, 4)))
  
  return(p)
}

###画到一张上
# 创建一个PDF文件来保存所有的图
pdf("boxplots_and_tests.pdf", width = 10, height = 14)

# 批量生成并绘制19个因子的图
plots <- lapply(unique(data_long$bio_factor), plot_and_test)

# 将ggplot对象转换为grob对象
grob_list <- lapply(plots, ggplotGrob)

n_cols <- 4
n_rows <- ceiling(length(grob_list) / n_cols)

# 将所有图排列在一页PDF文件中
grid.arrange(grobs = grob_list, ncol = n_cols, nrow = n_rows)
# 关闭PDF文件
dev.off()

###每个因子单独一张
factor_names <- unique(data_long$bio_factor)
# 循环为每个因子生成单独的PDF文件
for (factor_name in factor_names) {
  # 创建PDF文件，以因子名为文件名
  pdf(paste0(factor_name, ".pdf")) 
  
  # 绘制当前因子的箱线图
  p <- plot_and_test(factor_name)
  ggsave(p, filename = paste0(factor_name, ".pdf"), height = 5, width = 4)
  
  # 关闭当前的PDF文件
  dev.off() 
}


###环境因子和微生物多样性相关性
library(ggplot2)
library(corrplot)
setwd("E:/Project/13-蜘蛛微生物/0-data/")
data <- read.csv("alpha_diversity.csv", header = T)

# 提取需要计算相关性的列
cor_data <- data[, c("shannon", "simpson", "faith", "chao1", "Temperature", "Precipitation")]

# 计算相关性矩阵
cor_matrix <- cor(cor_data)

# 绘制相关性矩阵图
corrplot(cor_matrix, method = "circle", type = "upper", 
         addCoef.col = "black", tl.col = "black", tl.srt = 45)

ggsave(filename = "E:/Project/13-蜘蛛微生物/1-result/FIG/figure4/cor.pdf", height = 4.5, width = 7)
# 计算 Shannon 指数与气温的相关性
shannon_temp_cor <- cor.test(data$shannon, data$Temperature)
print(shannon_temp_cor)

# 计算 Shannon 指数与降雨量的相关性
shannon_precip_cor <- cor.test(data$shannon, data$Precipitation)
print(shannon_precip_cor)

# 计算 Simpson 指数与气温的相关性
simpson_temp_cor <- cor.test(data$simpson, data$Temperature)
print(simpson_temp_cor)

# 计算 Simpson 指数与降雨量的相关性
simpson_precip_cor <- cor.test(data$simpson, data$Precipitation)
print(simpson_precip_cor)


#散点拟合
model <- lm(shannon ~ Precipitation, data = data)
p_value <- summary(model)$coefficients["Precipitation", "Pr(>|t|)"]

r_squared <- summary(model)$r.squared
adj_r_squared <- summary(model)$adj.r.squared

p <- ggplot(data, aes(x = Precipitation, y = shannon)) +
  geom_point(color = "#b0d5df") +
  geom_smooth(method = "lm", se = FALSE, color = "#1781b5") +
  annotate("text", x = max(data$Precipitation), y = max(data$shannon-1), 
           label = paste("p-value =", format(p_value), 
                         "\nR-squared =", format(r_squared, digits = 3),
                         "\nAdj R-squared =", format(adj_r_squared, digits = 3)), 
           hjust = 1) +
  labs(x = "Precipitation", y = "Shannon Index", 
       title = "Relationship between Shannon Index and Precipitation")+
  theme(
    axis.title = element_text(size = 14),  # 修改坐标轴标题字体大小
    axis.text = element_text(size = 12),   # 修改坐标轴刻度字体大小
    plot.title = element_text(size = 16),  # 修改图标题字体大小（如果有）
    legend.text = element_text(size = 14), # 修改图例文字字体大小
    legend.title = element_text(size = 14) # 修改图例标题字体大小
  )
  

p
ggsave(p, filename = "E:/Project/13-蜘蛛微生物/1-result/FIG/figure4/Precipitation_shannon.pdf", height = 4, width = 6)

model <- lm(simpson ~ Precipitation, data = data)
p_value <- summary(model)$coefficients["Precipitation", "Pr(>|t|)"]

r_squared <- summary(model)$r.squared
adj_r_squared <- summary(model)$adj.r.squared

p <- ggplot(data, aes(x = Precipitation, y = simpson)) +
  geom_point(color = "#d1c2d3") +
  geom_smooth(method = "lm", se = FALSE, color = "#525288") +
  annotate("text", x = max(data$Precipitation), y = max(data$simpson-0.1), 
           label = paste("p-value =", format(p_value),
                         "\nR-squared =", format(r_squared, digits = 3),
                         "\nAdj R-squared =", format(adj_r_squared, digits = 3)), 
           hjust = 1) +
  labs(x = "Precipitation", y = "Simpson Index", 
       title = "Relationship between Simpson Index and Precipitation")+
  theme(
    axis.title = element_text(size = 14),  # 修改坐标轴标题字体大小
    axis.text = element_text(size = 12),   # 修改坐标轴刻度字体大小
    plot.title = element_text(size = 16),  # 修改图标题字体大小（如果有）
    legend.text = element_text(size = 14), # 修改图例文字字体大小
    legend.title = element_text(size = 14) # 修改图例标题字体大小
  )

p
ggsave(p, filename = "E:/Project/13-蜘蛛微生物/1-result/FIG/figure4/Precipitation_simpson.pdf", height = 4, width = 6)


model <- lm(faith ~ Precipitation, data = data)
p_value <- summary(model)$coefficients["Precipitation", "Pr(>|t|)"]

r_squared <- summary(model)$r.squared
adj_r_squared <- summary(model)$adj.r.squared

p <- ggplot(data, aes(x = Precipitation, y = faith)) +
  geom_point(color = "#f1c4cd") +
  geom_smooth(method = "lm", se = FALSE, color = "#cc163a") +
  annotate("text", x = max(data$Precipitation), y = max(data$faith-30), 
           label = paste("p-value =", format(p_value), 
           "\nR-squared =", format(r_squared, digits = 3),
           "\nAdj R-squared =", format(adj_r_squared, digits = 3)),
           hjust = 1) +
  labs(x = "Precipitation", y = "Faith-PD Index", 
       title = "Relationship between Faith-PD Index and Precipitation")+
  theme(
    axis.title = element_text(size = 14),  # 修改坐标轴标题字体大小
    axis.text = element_text(size = 12),   # 修改坐标轴刻度字体大小
    plot.title = element_text(size = 16),  # 修改图标题字体大小（如果有）
    legend.text = element_text(size = 14), # 修改图例文字字体大小
    legend.title = element_text(size = 14) # 修改图例标题字体大小
  )

p
ggsave(p, filename = "E:/Project/13-蜘蛛微生物/1-result/FIG/Precipitation_faith.pdf", height = 4, width = 6)

model <- lm(chao1 ~ Precipitation, data = data)
p_value <- summary(model)$coefficients["Precipitation", "Pr(>|t|)"]

r_squared <- summary(model)$r.squared
adj_r_squared <- summary(model)$adj.r.squared

p <- ggplot(data, aes(x = Precipitation, y = chao1)) +
  geom_point(color = "#92b3a5") +
  geom_smooth(method = "lm", se = FALSE, color = "#2bae85") +
  annotate("text", x = max(data$Precipitation), y = max(data$chao1-300), 
           label = paste("p-value =", format(p_value),
                         "\nR-squared =", format(r_squared, digits = 3),
                         "\nAdj R-squared =", format(adj_r_squared, digits = 3)), 
           hjust = 1) +
  labs(x = "Precipitation", y = "chao1 Index", 
       title = "Relationship between chao1 Index and Precipitation")

p
ggsave(p, filename = "E:/Project/13-蜘蛛微生物/1-result/FIG/Precipitation_chao1.pdf", height = 4, width = 6)

model <- lm(shannon ~ NDVI, data = data)
p_value <- summary(model)$coefficients["NDVI", "Pr(>|t|)"]
p <- ggplot(data, aes(x = NDVI, y = shannon)) +
  geom_point(color = "#92b3a5") +
  geom_smooth(method = "lm", se = FALSE, color = "#2bae85") +
  annotate("text", x = max(data$NDVI), y = max(data$shannon), 
           label = paste("p-value =", format(p_value),
                         "\nR-squared =", format(r_squared, digits = 3),
                         "\nAdj R-squared =", format(adj_r_squared, digits = 3)), 
           hjust = 1) +
  labs(x = "NDVI", y = "shannon Index", 
       title = "Relationship between shannon Index and NDVI")+
  theme(
    axis.title = element_text(size = 14),  # 修改坐标轴标题字体大小
    axis.text = element_text(size = 12),   # 修改坐标轴刻度字体大小
    plot.title = element_text(size = 16),  # 修改图标题字体大小（如果有）
    legend.text = element_text(size = 14), # 修改图例文字字体大小
    legend.title = element_text(size = 14) # 修改图例标题字体大小
  )

p
ggsave(p, filename = "NDVI_shannon.pdf", height = 4, width = 6)


###地理距离#########
library(vegan)
library(geosphere)
library(ggplot2)
setwd("E:/Project/13-蜘蛛微生物/0-data/")
otu_table <- read.csv("Genus_abundance_all.csv", header = T, row.names = 1)
otu_table <- t(otu_table)
df_geo <- read.csv("site_geo.csv", header = T)
metadata <- read.csv("sample_info.csv", header = T)

# 1. 计算地理距离矩阵（假设df_geo包含经纬度）
geo_dist <- distm(df_geo[, c("long", "lat")], fun = distHaversine) # 单位：米
geo_dist <- as.dist(geo_dist / 1000) # 转换为千米

# 2. 计算Bray-Curtis距离
bray_dist <- vegdist(otu_table, method = "bray")

# 3. Mantel检验
mantel_result <- mantel(geo_dist, bray_dist, permutations = 9999)
print(mantel_result)

r_value <- mantel_result$statistic
p_value <- mantel_result$signif

# 4. PERMANOVA（省份分组）
# adonis_result <- adonis2(bray_dist ~ province, data = metadata, permutations = 9999)
# print(adonis_result)

# 5. 距离衰减图
df_plot <- data.frame(
  geo = as.vector(geo_dist),
  bray = as.vector(bray_dist)
)

gg <- ggplot(df_plot, aes(x = geo, y = 1 - bray)) +
  geom_point(color = "#f7cdbc", size = 1, alpha = 0.2) +
  geom_smooth(method = "lm",color = "#f2481b", se = FALSE, size = 2) +
  annotate("text", x = 600, y = 1, 
           label = paste("Mantel r =", round(r_value, 3), "\nP-value =", format(p_value, scientific = TRUE, digits = 3)),
           hjust = 1, vjust = 1, size = 4) +
  labs(x = "Geographic Distance (km)", y = "Microbial Similarity (1 - Bray-Curtis)") +
  theme(
    axis.title = element_text(size = 14),  # 修改坐标轴标题字体大小
    axis.text = element_text(size = 12),   # 修改坐标轴刻度字体大小
    plot.title = element_text(size = 16),  # 修改图标题字体大小（如果有）
    legend.text = element_text(size = 14), # 修改图例文字字体大小
    legend.title = element_text(size = 14) # 修改图例标题字体大小
  )

gg

ggsave("E:/Project/13-蜘蛛微生物/1-result/FIG/figure4/geo.png", plot = gg, height = 3.5, width = 5, dpi = 300)
