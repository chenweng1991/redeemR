#' Function to plot the mito depth summary
#'
#' This function allows you to plot both position-wise and cell-wise mito depth summary
#' @param depth The .depth file by function DepthSummary
#' @param name The plot name shown on top
#' @param w the Width of the plot, default=10
#' @param h the height of the plot default=3
#' @return directly out put the plot
#' @examples
#' plot_depth(DN1CD34_1.depth$Total,"Total")
#' @export
#' @import ggplot2
plot_depth<-function(depth=DN1CD34_1.depth,name="",w=10,h=3){
d1<-depth[[1]]
d2<-depth[[2]]
names(d1)<-c("pos","meanCov")
names(d2)<-c("cell","meanCov")
options(repr.plot.width=w, repr.plot.height=h)
p1<-ggplot(d1)+aes(pos,meanCov)+geom_point()+theme_bw()
p2<-ggplot(d1)+aes("cell",meanCov)+geom_violin()+geom_boxplot()+theme_bw()
grid.arrange(p1,p2,layout_matrix=rbind(c(1,1,1,1,1,1,2)),top=name)
}


#' Function to plot bulk level mutation signatures
#'
#' This function allows you to plot the mito mutation signatures
#' @param cell_variants a vector of variants formated as c('93_A_G''103_G_A''146_T_C'
#' @return p from ggplot2
#' @examples
#' MutationProfile.bulk(DN1CD34_1.Variants.feature.lst[[name]]$Variants
#' @export
#' @import ggplot2
MutationProfile.bulk<-function(cell_variants){  ## cell_variants look like c('93_A_G''103_G_A''146_T_C''150_C_T''152_T_C''182_C_T')
# Annotate with called variants
data(ref_all_long)
called_variants <- strsplit(cell_variants,"_") %>% sapply(.,function(x){paste(x[1],x[2],">",x[3],sep="")})
ref_all_long$called <- ref_all_long$variant %in% called_variants
# Compute changes in expected/observed
total <- dim(ref_all_long)[1]
total_called <- sum(ref_all_long$called)
prop_df <- ref_all_long %>% group_by(three_plot, group_change, strand) %>%
  dplyr::summarize(observed_prop_called = sum(called)/total_called, expected_prop = n()/total, n = n()) %>%
  mutate(fc_called = observed_prop_called/expected_prop)
prop_df$change_plot <- paste0(prop_df$group_change, "_", prop_df$three_plot)
# Visualize
p1 <- ggplot(prop_df, aes(x = change_plot, fill = strand, y = fc_called)) +
  geom_bar(stat = "identity", position = "dodge") + pretty_plot(fontsize = 8) + L_border() +
  theme(axis.title.x=element_blank(),
        axis.text.x =element_blank())+
  scale_fill_manual(values= c("firebrick", "dodgerblue3")) +
  theme(legend.position = "bottom") +
  scale_y_continuous(expand = c(0,0)) +
  geom_hline(yintercept = 1, linetype =2, color = "black") +
  labs(x = "Change in nucleotide", y = "Substitution Rate")+
   facet_grid(.~group_change,scales = "free",space="free")
return(p1)
}
