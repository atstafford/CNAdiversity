# DIFFERENCE IN PGA AND CLONALITY/FREQ OF ALTERATIONS BETWEEN AD AND CAR

# t-test to compare the average PGA between adenomas and carcinomas:
t.test(car.diversity$pga$prop.aneu, ad.diversity$pga$prop.aneu,
       alternative = 'two.sided', var.equal = FALSE)

# t-test and violin to compare the average PIC frac between adenomas and carcinomas:
t.test(car.diversity$pic.frac$pic.frac, ad.diversity$pic.frac$pic.frac,
       alternative = 'two.sided', var.equal = FALSE)

x <- rbind(car.diversity$pic.frac, ad.diversity$pic.frac)
x$type <- c(rep('car',81),rep('ad',19))
ggplot(x, aes(x=type, y=pic.frac, colour=type)) +
  geom_violin(trim=FALSE, size=1) +
  geom_jitter(shape=16, size = 5, position = position_jitter(0.1)) +
  stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.05, colour='black', size = 1)+
  theme(plot.margin = unit(c(t=0,r=0,b=0,l=0), "cm"),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black", size = 0.5),
        axis.title = element_text(size=10),
        axis.text = element_text(size=10),
        legend.position = "none") 

# carcinoma vs lung
t.test(car.diversity$pga$prop.aneu, tracerx.diversity$pga$prop.aneu,
       alternative = 'two.sided', var.equal = FALSE)

x <- rbind(car.diversity$pga[4], tracerx.diversity$pga[4])
x$type <- c(rep('CRC',car.info$noSamples),rep('NSCLC',tracerx.info$noSamples))
plot <- ggplot(x, aes(x=type, y=prop.aneu, colour=type)) +
  geom_violin(trim=FALSE, size=1) +
  geom_jitter(shape=16, size = 5, position = position_jitter(0.1)) +
  stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.05, colour='black', size = 1)+
  theme(plot.margin = unit(c(t=0,r=0,b=0,l=0), "cm"),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black", size = 0.5),
        axis.title = element_text(size=10),
        axis.text = element_text(size=10),
        legend.position = "none") 
jpeg('tempfig.jpeg', width = 750, height = 750)
plot
dev.off()

t.test(car.diversity$pic.frac, tracerx.diversity$pic.frac,
       alternative = 'two.sided', var.equal = FALSE)

x <- rbind(car.diversity$pic.frac, tracerx.diversity$pic.frac)
x$type <- c(rep('CRC',car.info$noPatients),rep('NSCLC',tracerx.info$noPatients))
plot <- ggplot(x, aes(x=type, y=pic.frac, colour=type)) +
  geom_violin(trim=FALSE, size=1) +
  geom_jitter(shape=16, size = 5, position = position_jitter(0.1)) +
  stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.05, colour='black', size = 1)+
  theme(plot.margin = unit(c(t=0,r=0,b=0,l=0), "cm"),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black", size = 0.5),
        axis.title = element_text(size=10),
        axis.text = element_text(size=10),
        legend.position = "none") 
jpeg('tempfig.jpeg', width = 750, height = 750)
plot
dev.off()

# train vs val
t.test(car.diversity$pga$prop.aneu, val.diversity$pga$prop.aneu,
       alternative = 'two.sided', var.equal = FALSE)

x <- rbind(car.diversity$pga[4], val.diversity$pga[4])
x$type <- c(rep('CRC',car.info$noSamples),rep('NSCLC',val.info$noSamples))
plot <- ggplot(x, aes(x=type, y=prop.aneu, colour=type)) +
  geom_violin(trim=FALSE, size=1) +
  geom_jitter(shape=16, size = 5, position = position_jitter(0.1)) +
  stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.05, colour='black', size = 1)+
  theme(plot.margin = unit(c(t=0,r=0,b=0,l=0), "cm"),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black", size = 0.5),
        axis.title = element_text(size=10),
        axis.text = element_text(size=10),
        legend.position = "none") 
jpeg('tempfig.jpeg', width = 750, height = 750)
plot
dev.off()

t.test(car.diversity$pic.frac, val.diversity$pic.frac,
       alternative = 'two.sided', var.equal = FALSE)

x <- rbind(car.diversity$pic.frac, val.diversity$pic.frac)
x$type <- c(rep('CRC',car.info$noPatients),rep('NSCLC',val.info$noPatients))
plot <- ggplot(x, aes(x=type, y=pic.frac, colour=type)) +
  geom_violin(trim=FALSE, size=1) +
  geom_jitter(shape=16, size = 5, position = position_jitter(0.1)) +
  stat_summary(fun.data="mean_sdl", geom="crossbar", width=0.05, colour='black', size = 1)+
  theme(plot.margin = unit(c(t=0,r=0,b=0,l=0), "cm"),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black", size = 0.5),
        axis.title = element_text(size=10),
        axis.text = element_text(size=10),
        legend.position = "none") 
jpeg('tempfig.jpeg', width = 750, height = 750)
plot
dev.off()


# t-test to compare the number of clonal CNAs between adenomas and carcinomas, as a proportion of PGA
t.test(car.clonality$patientClo$clonal/car.clonality$patientClo$CNA, ad.clonality$patientClo$clonal/ad.clonality$patientClo$CNA,
       alternative = 'two.sided',var.equal = FALSE)

# t-test to compare the number of subclonal CNAs between adenomas and carcinomas, as a proportion of PGA
t.test(car.clonality$patientClo$subclonal/car.clonality$patientClo$CNA, ad.clonality$patientClo$subclonal/ad.clonality$patientClo$CNA,
       alternative = 'two.sided',var.equal = FALSE)

# t-test to compare the number of clonal and subclonal CNAs in adenoma patients
t.test(ad.clonality$patientClo$subclonal/ad.clonality$patientClo$CNA,ad.clonality$patientClo$clonal/ad.clonality$patientClo$CNA,
       alternative = 'two.sided',var.equal = FALSE)

# t-test to compare the number of clonal and subclonal CNAs in carcinoma patients
t.test(car.clonality$patientClo$subclonal/car.clonality$patientClo$CNA,car.clonality$patientClo$clonal/car.clonality$patientClo$CNA,
       alternative = 'two.sided',var.equal = FALSE)



