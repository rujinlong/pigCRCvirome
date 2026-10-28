diet_colors <- c(CTRL = "#357EBD", RL = "#D43F3A", COL = "#B87910")
diet_shapes <- c(CTRL = 16, RL = 17, COL = 15)
p_label <- function(p) ifelse(p < 0.001, "P < 0.001", sprintf("P = %.3f", p))
pub_theme <- function() theme_classic(base_size=8,base_family="Arial") +
 theme(plot.title=element_text(size=9,face="bold"),plot.subtitle=element_text(size=7),
 plot.title.position="plot",legend.position="bottom",legend.title=element_blank(),
 strip.background=element_blank(),strip.text=element_text(size=7.5),axis.text=element_text(color="black"),
 plot.tag=element_text(face="bold",size=11),plot.margin=margin(5,6,5,6))
save_pub_r <- function(plot,filename,width_mm=183,height_mm=210,dpi=600) {
 w <- width_mm/25.4; h <- height_mm/25.4
 svglite::svglite(paste0(filename,".svg"),width=w,height=h,system_fonts=list(sans="Arial")); print(plot); grDevices::dev.off()
 grDevices::cairo_pdf(paste0(filename,".pdf"),width=w,height=h,family="Arial"); print(plot); grDevices::dev.off()
 ragg::agg_tiff(paste0(filename,".tiff"),width=w,height=h,units="in",res=dpi,compression="lzw",background="white"); print(plot); grDevices::dev.off()
 ggsave(paste0(filename,".png"),plot,width=w,height=h,dpi=300,device=ragg::agg_png,bg="white")
}
