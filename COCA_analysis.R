################################################################################
# This file has been adapted from the code by Przepiórkowski and Woźniak, 2023
# With modifications for the purpose of this study.
#############################################################################################################
### LOADING LIBRARIES
#############################################################################################################
required.packages = c("dplyr",
                      "xtable",
                      "effects",
                      "Cairo",
                      "gridExtra",
                      "janitor")
for(i in required.packages) {if(!require(i,character.only = TRUE)) install.packages(i)} 
for(i in required.packages) {library(i,character.only = TRUE)}
rm(required.packages, i)
#############################################################################################################
### LOADING DATA SOURCE
#############################################################################################################

setwd("C:\\Users\\DELL\\OneDrive\\Dokumenty\\BNP_PROJECT\\BNL\\split_COCA\\csv\\parse_raw\\ccp\\results\\analysis")
coords <- read.csv("coords.csv", header=TRUE, comment.char="", stringsAsFactors=TRUE)

# ensure there are no duplicates in data
coords <- distinct(coords)
print(nrow(coords))

# removing coulmns that are not analysed
coords <- coords %>% 
   select(-sentence.tree, -converted.from.file)

#############################################################################################################
### COMPUTING CONJUNCT LENGTH DIFFERENCE
#############################################################################################################

diff <- coords$R.words - coords$L.words
coords$words.diff <- abs(diff)
coords$words.shorter <- case_when(diff < 0 ~ "R",
                                  diff == 0 ~ "0",
                                  diff > 0 ~ "L")
coords$words.shorter <- factor(coords$words.shorter, levels=c("R","0","L"))

diff <- coords$R.chars - coords$L.chars
coords$chars.diff <- abs(diff)
coords$chars.shorter <- case_when(diff < 0 ~ "R",
                                  diff == 0 ~ "0",
                                  diff > 0 ~ "L")
coords$chars.shorter <- factor(coords$chars.shorter, levels=c("R","0","L"))


diff <- coords$R.syllables - coords$L.syllables
coords$syllables.diff <- abs(diff)
coords$syllables.shorter <- case_when(diff < 0 ~ "R",
                                      diff == 0 ~ "0",
                                      diff > 0 ~ "L")
coords$syllables.shorter <- factor(coords$syllables.shorter, levels=c("R","0","L"))


# diff <- coords$R.tokens - coords$L.tokens
# coords$tokens.diff <- abs(diff)
# coords$tokens.shorter <- case_when(diff < 0 ~ "R",
#                                       diff == 0 ~ "0",
#                                       diff > 0 ~ "L")
# coords$tokens.shorter <- factor(coords$tokens.shorter, levels=c("R","0","L"))

head(coords)

#############################################################################################################
### ANALYSIS
#############################################################################################################

################## number of coords ##################
total_length <- nrow(coords)
genre_count <- coords %>% 
   group_by(genre) %>% 
   count()
genre_count <- genre_count %>%
   arrange(desc(n)) %>% 
   mutate("%" = round(n/total_length*100,2)) %>% 
   adorn_totals("row")
names(genre_count) <- c("genre", "N", "%")
print(xtable(genre_count,
             label = "tab:coord:count"),
      align="lrr",
      include.rownames=FALSE,
      file = "tab_no_coord.tex"
)


################## different coordinations ##################
coordinators <- coords %>% 
   ungroup( ) %>% 
   select(conjunction.word) %>% 
   group_by(conjunction.word) %>% 
   count() %>% 
   arrange(desc(n)) %>% 
   mutate("%" = round(n/total_length*100,3)) 
coordinators <- coordinators %>% 
   ungroup() %>% 
   mutate(n_sum = cumsum(n)) 
coordinators <- coordinators %>% 
   mutate("sum%" = round(n_sum/total_length*100,3)) %>% 
   select(-n_sum)
write.csv(coordinators, "coordinators.csv")


################## analyzing number of coordination by governor side ##################
gov_side <- coords %>% 
   ungroup() %>% 
   group_by(governor.position) %>% 
   count() %>% 
   arrange(desc(n)) %>% 
   mutate("%" = round(n/total_length*100,2))
names(gov_side) <- c(colnames(gov_side)[1], "N", "%")
print(xtable(gov_side,
             label = "tab:gov_side"),
      align="lrr",
      include.rownames=FALSE,
      file = "tab_no_gov_side.tex"
)

#################################################
############# GENRE
#################################################

acad <- coords %>% 
   filter(genre == "acad")
blog <- coords %>% 
   filter(genre == "blog")
fic <- coords %>% 
   filter(genre == "fic")
mag <- coords %>% 
   filter(genre == "mag")
news <- coords %>% 
   filter(genre == "news")
tvm <- coords %>% 
   filter(genre == "tvm")
web <- coords %>% 
   filter(genre == "web")
genres = list("acad" = acad,
              "blog" = blog,
              "fic" = fic,
              "mag" = mag,
              "news" = news,
              "tvm" = tvm,
              "web" = web)

#####################################################
no.genres <- length(genres)

## Create data frames and LaTeX tables for each of the above four populations of coordinations:

tabs <- list()  #

addtorow <- list()
addtorow$pos <- list(0, 0)
addtorow$command <- c("& \\multicolumn{2}{c}{\\textls[50]{\\,median}} & \\multicolumn{2}{c}{\\textls[500]{mean}}  \\\\\n",
                      "& \\multicolumn{1}{c}{left} & \\multicolumn{1}{c}{\\hspace*{-1ex}right\\hspace*{-1ex}} & \\multicolumn{1}{c}{left} & \\multicolumn{1}{c}{right} & \\multicolumn{1}{c}{V} & \\multicolumn{1}{c}{$p$} \\\\\n")

for (i in 1:no.genres) {
   
   yy <- genres[[i]]
   
   tests <- list(
      wilcox.test(yy$L.chars, yy$R.chars, alternative="less", paired=TRUE, exact=FALSE, correct=FALSE),
      wilcox.test(yy$L.syllables, yy$R.syllables, alternative="less", paired=TRUE, exact=FALSE, correct=FALSE),
      wilcox.test(yy$L.words, yy$R.words, alternative="less", paired=TRUE, exact=FALSE, correct=FALSE))
   tabs[[i]] <- data.frame(measure = c('characters','syllables','words'),
                           med.1 = c(median(yy$L.chars), median(yy$L.syllables), median(yy$L.words)),
                           med.2 = c(median(yy$R.chars), median(yy$R.syllables), median(yy$R.words)),
                           mean.1 = c(mean(yy$L.chars), mean(yy$L.syllables), mean(yy$L.words)),
                           mean.2 = c(mean(yy$R.chars), mean(yy$R.syllables), mean(yy$R.words)),
                           V = c(tests[[1]]$statistic, tests[[2]]$statistic, tests[[3]]$statistic),
                           p = c(tests[[1]]$p.value, tests[[2]]$p.value, tests[[3]]$p.value))
   
   addtorow$pos <- append(addtorow$pos, rep(3*(i-1), 3))
   subtitle <- paste("\\multicolumn{7}{c}{\\emph{", names(genres)[i], " (N = ", format(nrow(yy),big.mark=",",scientific=FALSE), ")}} \\\\\n", sep="")
   addtorow$command <- append(addtorow$command, c("\\midrule\n", subtitle, ifelse(i>1, "\\midrule\n", "")))
}

## Combine these data frames and LaTeX tables into one data frame and one LaTeX table:

tab <- tabs[[1]]
for (i in 2:no.genres) tab <- rbind(tab,tabs[[i]])
print(
   xtable(tab,
          align="llrrrrrr",
          display=c("s","s","f","f","f","f","g","g"),
          digits=c(0,0,0,0,2,2,2,2),
          caption="Medians and means of lengths of left and right conjuncts in binary coordinations in COCA by genre.",
          label="tab:basic:genre"),
   size="\\setlength{\\tabcolsep}{4pt}",
   include.colnames=FALSE,
   include.rownames=FALSE,
   add.to.row=addtorow,
   file="tab_basic_genre.tex")

########################################################################
################## distribution of governors by genre ##################
########################################################################
# Created using code by Przepiórkowski and Woźniak, 2023

genre.names <- c("acad", "blog", "fic", "mag", "news", "tvm", "web")
coords <- coords %>% 
   filter(genre %in% genre.names)

oo <- coords %>% 
   filter(governor.position =="0")
ll <-  coords %>% 
   filter(governor.position =="L")
rr <-  coords %>% 
   filter(governor.position =="R")

ll.genres <- summary(droplevels(ll$genre))
rr.genres <- summary(droplevels(rr$genre))
oo.genres <- summary(droplevels(oo$genre))

ll.genres.perc <- 100*prop.table(table(ll$genre))
rr.genres.perc <- 100*prop.table(table(rr$genre))
oo.genres.perc <- 100*prop.table(table(oo$genre))

empty <- rep("~~~~",length(genre.names))

tab <- data.frame(empty, oo.genres, oo.genres.perc, empty, ll.genres, ll.genres.perc, empty, rr.genres, rr.genres.perc)[,c(2,4,6,8,10,12)]
rownames(tab) <- genre.names
addtorow <- list()
addtorow$pos <- list(0, 0, 0)
addtorow$command <- c("& \\multicolumn{6}{c}{\\textls[1000]{governor}}  \\\\\n",
                      "& \\multicolumn{2}{r}{\\textls[600]{none}} & \\multicolumn{2}{r}{\\textls[100]{on the left}} & \\multicolumn{2}{r}{\\textls[100]{on the right}}  \\\\\n",
                      "& \\# & \\% & \\# & \\% & \\# & \\% \\\\\n")
print(
   xtable(tab,
          align="lrrrrrr",
          display=c("s","d","f","d","f","d","f"),
          digits=c(0,0,2,0,2,0,2),
          caption="Numbers (\\#) and percentages (\\%) of coordinations of different COCA genres depending on the presence and position of the governor",
          label="tab:gov:genres"),
   include.colnames=FALSE,
   add.to.row=addtorow,
   file="tab_genres_gov.tex")


################## conj length ##################

######################################################################
## Code for generating (the LaTeX sources of) Table 4.5
## Created using code by Przepiórkowski and Woźniak, 2023
######################################################################

populations <- list("All coordinations"=coords,
                    "No governor"=oo,
                    "Governor on the left"=ll,
                    "Governor on the right"=rr)
no.populations <- length(populations)

## Create data frames and LaTeX tables for each of the above four populations of coordinations:

tabs <- list()  #

addtorow <- list()
addtorow$pos <- list(0, 0)
addtorow$command <- c("& \\multicolumn{2}{c}{\\textls[50]{\\,median}} & \\multicolumn{2}{c}{\\textls[500]{mean}}  \\\\\n",
                      "& \\multicolumn{1}{c}{left} & \\multicolumn{1}{c}{\\hspace*{-1ex}right\\hspace*{-1ex}} & \\multicolumn{1}{c}{left} & \\multicolumn{1}{c}{right} & \\multicolumn{1}{c}{V} & \\multicolumn{1}{c}{$p$} \\\\\n")

for (i in 1:no.populations) {
   
   yy <- populations[[i]]
   tests <- list(
      wilcox.test(yy$L.chars, yy$R.chars, alternative="less", paired=TRUE, exact=TRUE, correct=FALSE),
      wilcox.test(yy$L.syllables, yy$R.syllables, alternative="less", paired=TRUE, exact=TRUE, correct=FALSE),
      wilcox.test(yy$L.words, yy$R.words, alternative="less", paired=TRUE, exact=TRUE, correct=FALSE))
   tabs[[i]] <- data.frame(measure = c('characters','syllables','words'),
                           med.1 = c(median(yy$L.chars), median(yy$L.syllables), median(yy$L.words)),
                           med.2 = c(median(yy$R.chars), median(yy$R.syllables), median(yy$R.words)),
                           mean.1 = c(mean(yy$L.chars), mean(yy$L.syllables), mean(yy$L.words)),
                           mean.2 = c(mean(yy$R.chars), mean(yy$R.syllables), mean(yy$R.words)),
                           V = c(tests[[1]]$statistic, tests[[2]]$statistic, tests[[3]]$statistic),
                           p = c(tests[[1]]$p.value, tests[[2]]$p.value, tests[[3]]$p.value))
   
   addtorow$pos <- append(addtorow$pos, rep(3*(i-1), 3))
   subtitle <- paste("\\multicolumn{7}{c}{\\emph{", names(populations)[i], " (N = ", format(nrow(yy),big.mark=",",scientific=FALSE), ")}} \\\\\n", sep="")
   addtorow$command <- append(addtorow$command, c("\\midrule\n", subtitle, ifelse(i>1, "\\midrule\n", "")))
}

## Combine these data frames and LaTeX tables into one data frame and one LaTeX table:

tab <- tabs[[1]]
for (i in 2:no.populations) tab <- rbind(tab,tabs[[i]])
print(
   xtable(tab,
          align="llrrrrrr",
          display=c("s","s","f","f","f","f","g","g"),
          digits=c(0,0,0,0,2,2,2,2),
          caption="Medians and means of lengths of left and right conjuncts in binary coordinations in COCA",
          label="tab:basic"),
   size="\\setlength{\\tabcolsep}{4pt}",
   include.colnames=FALSE,
   include.rownames=FALSE,
   add.to.row=addtorow,
   file="tab_basic.tex")


######################################################################
## Code for generating (the LaTeX sources of) Table 4.7
## Created using code by Przepiórkowski and Woźniak, 2023
######################################################################
tests <- list(
   prop.test(x = c(nrow(ll[ll$chars.shorter=="L",]), nrow(rr[rr$chars.shorter=="L",])),
             n = c(nrow(ll[ll$chars.shorter!="0",]), nrow(rr[rr$chars.shorter!="0",])),
             correct=FALSE),
   prop.test(x = c(nrow(ll[ll$syllables.shorter=="L",]), nrow(rr[rr$syllables.shorter=="L",])),
             n = c(nrow(ll[ll$syllables.shorter!="0",]), nrow(rr[rr$syllables.shorter!="0",])),
             correct=FALSE),
   prop.test(x = c(nrow(ll[ll$words.shorter=="L",]), nrow(rr[rr$words.shorter=="L",])),
             n = c(nrow(ll[ll$words.shorter!="0",]), nrow(rr[rr$words.shorter!="0",])),
             correct=FALSE))

tab <- matrix(c(
   tests[[1]]$estimate[1], tests[[2]]$estimate[1], tests[[3]]$estimate[1], 
   nrow(ll[ll$chars.shorter!="0",]), nrow(ll[ll$syllables.shorter!="0",]), nrow(ll[ll$words.shorter!="0",]),
   tests[[1]]$estimate[2], tests[[2]]$estimate[2], tests[[3]]$estimate[2], 
   nrow(rr[rr$chars.shorter!="0",]), nrow(rr[rr$syllables.shorter!="0",]), nrow(rr[rr$words.shorter!="0",]),
   tests[[1]]$statistic,   tests[[2]]$statistic,   tests[[3]]$statistic,
   tests[[1]]$p.value,     tests[[2]]$p.value,     tests[[3]]$p.value),
   ncol=6)

addtorow <- list()
addtorow$pos <- list(0, 0, 0)
addtorow$command <- c("& \\multicolumn{4}{c}{\\textls[1000]{\\hspace{-2pt}governor}} & & \\\\\n",
                      "& \\multicolumn{2}{c}{\\textls[100]{on the left}} & \\multicolumn{2}{c}{on the right} \\\\\n",
                      "& \\multicolumn{1}{c}{prop} & \\multicolumn{1}{c}{$N$} & \\multicolumn{1}{r}{prop} & \\multicolumn{1}{c}{$N$} & \\multicolumn{1}{c}{$\\chi^2(1)$} & \\multicolumn{1}{c}{$p$} \\\\\n")
rownames(tab) <- c('characters','syllables','words')
print(
   xtable(tab,
          align="lrrrrrr",
          display=c("s","f","f","f","f","f","g"),
          digits=c(0,3,0,3,0,1,2),
          caption="Proportions of shorter conjuncts occurring on the left (vs.~right) depending on the position of the governor, in coordinations with conjuncts of different lengths",
          label="tab:basic:props:lr"),
   size="\\setlength{\\tabcolsep}{4pt}",
   include.colnames=FALSE,
   add.to.row=addtorow,
   file="tab_basic_props_lr.tex")


######################################################################
## Code for generating (the LaTeX sources of) Table 4.8
## Created using code by Przepiórkowski and Woźniak, 2023
######################################################################

tests <- list(
   prop.test(x = c(nrow(ll[ll$chars.shorter=="L",]), nrow(oo[oo$chars.shorter=="L",])),
             n = c(nrow(ll[ll$chars.shorter!="0",]), nrow(oo[oo$chars.shorter!="0",])),
             correct=FALSE),
   prop.test(x = c(nrow(ll[ll$syllables.shorter=="L",]), nrow(oo[oo$syllables.shorter=="L",])),
             n = c(nrow(ll[ll$syllables.shorter!="0",]), nrow(oo[oo$syllables.shorter!="0",])),
             correct=FALSE),
   prop.test(x = c(nrow(ll[ll$words.shorter=="L",]), nrow(oo[oo$words.shorter=="L",])),
             n = c(nrow(ll[ll$words.shorter!="0",]), nrow(oo[oo$words.shorter!="0",])),
             correct=FALSE),
   prop.test(x = c(nrow(oo[oo$chars.shorter=="L",]), nrow(rr[rr$chars.shorter=="L",])),
             n = c(nrow(oo[oo$chars.shorter!="0",]), nrow(rr[rr$chars.shorter!="0",])),
             correct=FALSE),
   prop.test(x = c(nrow(oo[oo$syllables.shorter=="L",]), nrow(rr[rr$syllables.shorter=="L",])),
             n = c(nrow(oo[oo$syllables.shorter!="0",]), nrow(rr[rr$syllables.shorter!="0",])),
             correct=FALSE),
   prop.test(x = c(nrow(oo[oo$words.shorter=="L",]), nrow(rr[rr$words.shorter=="L",])),
             n = c(nrow(oo[oo$words.shorter!="0",]), nrow(rr[rr$words.shorter!="0",])),
             correct=FALSE))

tab <- matrix(c(
   tests[[1]]$estimate[2], tests[[2]]$estimate[2], tests[[3]]$estimate[2], 
   nrow(oo[oo$chars.shorter!="0",]), nrow(oo[oo$syllables.shorter!="0",]), nrow(oo[oo$words.shorter!="0",]),
   tests[[1]]$statistic,   tests[[2]]$statistic,   tests[[3]]$statistic,
   tests[[1]]$p.value,     tests[[2]]$p.value,     tests[[3]]$p.value,
   tests[[4]]$statistic,   tests[[5]]$statistic,   tests[[6]]$statistic,
   tests[[4]]$p.value,     tests[[5]]$p.value,     tests[[6]]$p.value),
   ncol=6)

addtorow <- list()
addtorow$pos <- list(0, 0)
addtorow$command <- c("& \\multicolumn{2}{l}{no governor} & \\multicolumn{2}{c}{\\textls[50]{vs.~on the left}} & \\multicolumn{2}{c}{vs.~on the right} \\\\\n",
                      "& \\multicolumn{1}{c}{prop} & \\multicolumn{1}{c}{$N$} & \\multicolumn{1}{c}{$\\chi^2(1)$} & \\multicolumn{1}{c}{$p$} & \\multicolumn{1}{c}{$\\chi^2(1)$} & \\multicolumn{1}{c}{$p$} \\\\\n")
rownames(tab) <- c('characters','syllables','words')
print(
   xtable(tab,
          align="lrrrrrr",
          display=c("s","f","f","f","g","f","g"),
          digits=c(0,3,0,1,2,1,2),
          caption="Proportions of shorter conjuncts occurring on the left (vs.~right) in the absence of governor, in coordinations with conjuncts of different lengths",
          label="tab:basic:props:o"),
   size="\\setlength{\\tabcolsep}{4pt}",
   include.colnames=FALSE,
   add.to.row=addtorow,
   file="tab_basic_props_o.tex")




######################################################################
## Code for generating PNGs used in Figure 4.3
## and for the multifactorial binary logistic regression analysis
## Created using code by Przepiórkowski and Woźniak, 2023
######################################################################
xx <- coords


range.ch <- seq(1,120, 1)
range.sy <- seq(1,40, 1)
range.wo <- seq(1,20, 1)

xx.ch <- xx %>% 
   filter(chars.shorter != "0" & chars.diff <= range.ch[length(range.ch)]) %>% 
   select(governor.position, chars.diff, chars.shorter) %>% 
   transmute(governor = governor.position, diff = chars.diff, shorter = chars.shorter)

xx.sy <- xx %>% 
   filter(syllables.shorter != "0" & syllables.diff <= range.sy[length(range.sy)]) %>% 
   select(governor.position, syllables.diff, syllables.shorter) %>% 
   transmute(governor = governor.position, diff = syllables.diff, shorter = syllables.shorter)

xx.wo <- xx %>% 
   filter(words.shorter != "0" & words.diff <= range.wo[length(range.wo)]) %>% 
   select(governor.position, words.diff, words.shorter) %>% 
   transmute(governor = governor.position, diff = words.diff, shorter = words.shorter)

table(xx.ch) # checking that each bucket is sufficiently large
table(xx.sy) # checking that each bucket is sufficiently large
table(xx.wo) # checking that each bucket is sufficiently large

governors <- c("0", "L", "R")
gov.texts <- c("NO governor", "Governor on the LEFT", "Governor on the RIGHT")
differences <- list("characters"=xx.ch, "syllables"=xx.sy, "words"=xx.wo)

plots <- list()
gov.texts <- c("NO governor", "Governor on the LEFT", "Governor on the RIGHT")
differences <- list("characters"=xx.ch, "syllables"=xx.sy, "words"=xx.wo)


plots <- list()



###############################################################################
######  The graphs are drawn one by one and have to be manually saved as png
######  files. This decision was due to hardware limitations as automatic pdf
######  files could not be created.
###############################################################################


########################## CHARACTERS
d <- 1
d.name <- names(differences)[d]
cat(paste("\n---------- DIFFERENCES (BUCKETS) OF LENGTHS measured in ", d.name, " ----------\n", sep=""))
D <- differences[[d]]

# multifactorial analysis for a given measure (length counted in d.name):
# 
print(summary(mm <- glm(shorter ~ 1 + diff * governor, data=D, family=binomial, na.action=na.exclude)))
print(drop1(mm, test="Chisq"))
print(pairs(emmeans::emtrends(mm, ~ governor, var="diff"), adjust="none"))

########################## NO GOVERNOR
{
   g <- 1
   cat("\n", paste(rep("-",40), collapse=""), governors[g], paste(rep("-",40), collapse=""), "\n")
   g.name <- governors[g]
   
   # monofactorial analysis for a given measure (d.name) and governor presence/position (g.name):
   
   print(summary(m <- glm(shorter ~ 1 + diff, data=D[D$governor==g.name,], family=binomial, na.action=na.exclude)))
   print(drop1(m, test="Chisq"))
   ph <- data.frame(lend <- effect("diff", m))
   
   plot.name <- paste(g.name,d.name,sep="/")
   
   m.slope <- sprintf("%.2e", summary(m)$coefficients[2,1])
   m.prob <- sprintf("%.3g", summary(m)$coefficients[2,4])
   
   p <- plot(lend, type="response", ylim=c(0.5, 0.95), grid=TRUE,
             key=list(corner = c(0.01, 0.98), cex=.8, col=c("gray10","gray40"),
                      text=list(lab=c(paste("slope: ", m.slope, sep=""), paste("p: ", m.prob, sep="")))),
             main=paste(gov.texts[g], " (length in ", toupper(d.name), ")", sep=""),
             xlab=paste("absolute difference in",d.name),
             ylab="proportion of shorter left conjuncts"); p
   p
   
}

########################## LEFT GOVERNOR
{
   g <- 2
   cat("\n", paste(rep("-",40), collapse=""), governors[g], paste(rep("-",40), collapse=""), "\n")
   g.name <- governors[g]
   
   # monofactorial analysis for a given measure (d.name) and governor presence/position (g.name):
   
   print(summary(m <- glm(shorter ~ 1 + diff, data=D[D$governor==g.name,], family=binomial, na.action=na.exclude)))
   print(drop1(m, test="Chisq"))
   ph <- data.frame(lend <- effect("diff", m))
   
   plot.name <- paste(g.name,d.name,sep="/")
   
   m.slope <- sprintf("%.2e", summary(m)$coefficients[2,1])
   m.prob <- sprintf("%.3g", summary(m)$coefficients[2,4])
   
   p <- plot(lend, type="response", ylim=c(0.1, 0.95), grid=TRUE,
             key=list(corner = c(0.01, 0.98), cex=.8, col=c("gray10","gray40"),
                      text=list(lab=c(paste("slope: ", m.slope, sep=""), paste("p: ", m.prob, sep="")))),
             main=paste(gov.texts[g], " (length in ", toupper(d.name), ")", sep=""),
             xlab=paste("absolute difference in",d.name),
             ylab="proportion of shorter left conjuncts"); p
   p
}

########################## RIGHT GOVERNOR
{
   g <- 3
   cat("\n", paste(rep("-",40), collapse=""), governors[g], paste(rep("-",40), collapse=""), "\n")
   g.name <- governors[g]
   
   # monofactorial analysis for a given measure (d.name) and governor presence/position (g.name):
   
   print(summary(m <- glm(shorter ~ 1 + diff, data=D[D$governor==g.name,], family=binomial, na.action=na.exclude)))
   print(drop1(m, test="Chisq"))
   ph <- data.frame(lend <- effect("diff", m))
   
   plot.name <- paste(g.name,d.name,sep="/")
   
   m.slope <- sprintf("%.2e", summary(m)$coefficients[2,1])
   m.prob <- sprintf("%.3g", summary(m)$coefficients[2,4])
   
   p <- plot(lend, type="response", ylim=c(0.5, 0.95), grid=TRUE,
             key=list(corner = c(0.01, 0.98), cex=.8, col=c("gray10","gray40"),
                      text=list(lab=c(paste("slope: ", m.slope, sep=""), paste("p: ", m.prob, sep="")))),
             main=paste(gov.texts[g], " (length in ", toupper(d.name), ")", sep=""),
             xlab=paste("absolute difference in",d.name),
             ylab="proportion of shorter left conjuncts"); p
   
   p
}


##################################
####################### SYLLABLES
##################################
d <- 2
d.name <- names(differences)[d]
cat(paste("\n---------- DIFFERENCES (BUCKETS) OF LENGTHS measured in ", d.name, " ----------\n", sep=""))
D <- differences[[d]]

# multifactorial analysis for a given measure (length counted in d.name):
# 
print(summary(mm <- glm(shorter ~ 1 + diff * governor, data=D, family=binomial, na.action=na.exclude)))
print(drop1(mm, test="Chisq"))
print(pairs(emmeans::emtrends(mm, ~ governor, var="diff"), adjust="none"))

########################## NO GOVERNOR
{
   g <- 1
   cat("\n", paste(rep("-",40), collapse=""), governors[g], paste(rep("-",40), collapse=""), "\n")
   g.name <- governors[g]
   
   # monofactorial analysis for a given measure (d.name) and governor presence/position (g.name):
   
   print(summary(m <- glm(shorter ~ 1 + diff, data=D[D$governor==g.name,], family=binomial, na.action=na.exclude)))
   print(drop1(m, test="Chisq"))
   ph <- data.frame(lend <- effect("diff", m))
   
   plot.name <- paste(g.name,d.name,sep="/")
   
   m.slope <- sprintf("%.2e", summary(m)$coefficients[2,1])
   m.prob <- sprintf("%.3g", summary(m)$coefficients[2,4])
   
   p <- plot(lend, type="response", ylim=c(0.5, 0.95), grid=TRUE,
             key=list(corner = c(0.01, 0.98), cex=.8, col=c("gray10","gray40"),
                      text=list(lab=c(paste("slope: ", m.slope, sep=""), paste("p: ", m.prob, sep="")))),
             main=paste(gov.texts[g], " (length in ", toupper(d.name), ")", sep=""),
             xlab=paste("absolute difference in",d.name),
             ylab="proportion of shorter left conjuncts"); p
   p
}

########################## LEFT GOVERNOR
{
   g <- 2
   cat("\n", paste(rep("-",40), collapse=""), governors[g], paste(rep("-",40), collapse=""), "\n")
   g.name <- governors[g]
   
   # monofactorial analysis for a given measure (d.name) and governor presence/position (g.name):
   
   print(summary(m <- glm(shorter ~ 1 + diff, data=D[D$governor==g.name,], family=binomial, na.action=na.exclude)))
   print(drop1(m, test="Chisq"))
   ph <- data.frame(lend <- effect("diff", m))
   
   plot.name <- paste(g.name,d.name,sep="/")
   
   m.slope <- sprintf("%.2e", summary(m)$coefficients[2,1])
   m.prob <- sprintf("%.3g", summary(m)$coefficients[2,4])
   
   p <- plot(lend, type="response", ylim=c(0.5, 0.95), grid=TRUE,
             key=list(corner = c(0.01, 0.98), cex=.8, col=c("gray10","gray40"),
                      text=list(lab=c(paste("slope: ", m.slope, sep=""), paste("p: ", m.prob, sep="")))),
             main=paste(gov.texts[g], " (length in ", toupper(d.name), ")", sep=""),
             xlab=paste("absolute difference in",d.name),
             ylab="proportion of shorter left conjuncts"); p
   p
}

########################## RIGHT GOVERNOR
{
   g <- 3
   cat("\n", paste(rep("-",40), collapse=""), governors[g], paste(rep("-",40), collapse=""), "\n")
   g.name <- governors[g]
   
   # monofactorial analysis for a given measure (d.name) and governor presence/position (g.name):
   
   print(summary(m <- glm(shorter ~ 1 + diff, data=D[D$governor==g.name,], family=binomial, na.action=na.exclude)))
   print(drop1(m, test="Chisq"))
   ph <- data.frame(lend <- effect("diff", m))
   
   plot.name <- paste(g.name,d.name,sep="/")
   
   m.slope <- sprintf("%.2e", summary(m)$coefficients[2,1])
   m.prob <- sprintf("%.3g", summary(m)$coefficients[2,4])
   
   p <- plot(lend, type="response", ylim=c(0.5, 0.95), grid=TRUE,
             key=list(corner = c(0.01, 0.98), cex=.8, col=c("gray10","gray40"),
                      text=list(lab=c(paste("slope: ", m.slope, sep=""), paste("p: ", m.prob, sep="")))),
             main=paste(gov.texts[g], " (length in ", toupper(d.name), ")", sep=""),
             xlab=paste("absolute difference in",d.name),
             ylab="proportion of shorter left conjuncts"); p
   p
}


##################################
####################### WORDS
##################################
d <- 3
d.name <- names(differences)[d]
cat(paste("\n---------- DIFFERENCES (BUCKETS) OF LENGTHS measured in ", d.name, " ----------\n", sep=""))
D <- differences[[d]]

# multifactorial analysis for a given measure (length counted in d.name):
# 
print(summary(mm <- glm(shorter ~ 1 + diff * governor, data=D, family=binomial, na.action=na.exclude)))
print(drop1(mm, test="Chisq"))
print(pairs(emmeans::emtrends(mm, ~ governor, var="diff"), adjust="none"))

########################## NO GOVERNOR
{
   g <- 1
   cat("\n", paste(rep("-",40), collapse=""), governors[g], paste(rep("-",40), collapse=""), "\n")
   g.name <- governors[g]
   
   # monofactorial analysis for a given measure (d.name) and governor presence/position (g.name):
   
   print(summary(m <- glm(shorter ~ 1 + diff, data=D[D$governor==g.name,], family=binomial, na.action=na.exclude)))
   print(drop1(m, test="Chisq"))
   ph <- data.frame(lend <- effect("diff", m))
   
   plot.name <- paste(g.name,d.name,sep="/")
   
   m.slope <- sprintf("%.2e", summary(m)$coefficients[2,1])
   m.prob <- sprintf("%.3g", summary(m)$coefficients[2,4])
   
   p <- plot(lend, type="response", ylim=c(0.5, 0.95), grid=TRUE,
             key=list(corner = c(0.01, 0.98), cex=.8, col=c("gray10","gray40"),
                      text=list(lab=c(paste("slope: ", m.slope, sep=""), paste("p: ", m.prob, sep="")))),
             main=paste(gov.texts[g], " (length in ", toupper(d.name), ")", sep=""),
             xlab=paste("absolute difference in",d.name),
             ylab="proportion of shorter left conjuncts"); p
   p
   
}

########################## LEFT GOVERNOR
{
   g <- 2
   cat("\n", paste(rep("-",40), collapse=""), governorss[g], paste(rep("-",40), collapse=""), "\n")
   g.name <- governors[g]
   
   # monofactorial analysis for a given measure (d .name) and governor presence/position (g.name):
   
   print(summary(m <- glm(shorter ~ 1 + diff, data=D[D$governor==g.name,], family=binomial, na.action=na.exclude)))
   print(drop1(m, test="Chisq"))
   ph <- data.frame(lend <- effect("diff", m))
   
   plot.name <- paste(g.name,d.name,sep="/")
   
   m.slope <- sprintf("%.2e", summary(m)$coefficients[2,1])
   m.prob <- sprintf("%.3g", summary(m)$coefficients[2,4])
   
   p <- plot(lend, type="response", ylim=c(0.5, 0.95), grid=TRUE,
             key=list(corner = c(0.01, 0.98), cex=.8, col=c("gray10","gray40"),
                      text=list(lab=c(paste("slope: ", m.slope, sep=""), paste("p: ", m.prob, sep="")))),
             main=paste(gov.texts[g], " (length in ", toupper(d.name), ")", sep=""),
             xlab=paste("absolute difference in",d.name),
             ylab="proportion of shorter left conjuncts"); p
   p
}

########################## RIGHT GOVERNOR
{
   g <- 3
   cat("\n", paste(rep("-",40), collapse=""), governors[g], paste(rep("-",40), collapse=""), "\n")
   g.name <- governors[g]
   
   # monofactorial analysis for a given measure (d.name) and governor presence/position (g.name):
   
   print(summary(m <- glm(shorter ~ 1 + diff, data=D[D$governor==g.name,], family=binomial, na.action=na.exclude)))
   print(drop1(m, test="Chisq"))
   ph <- data.frame(lend <- effect("diff", m))
   
   plot.name <- paste(g.name,d.name,sep="/")
   
   m.slope <- sprintf("%.2e", summary(m)$coefficients[2,1])
   m.prob <- sprintf("%.3g", summary(m)$coefficients[2,4])
   
   p <- plot(lend, type="response", ylim=c(0.5, 0.95), grid=TRUE,
             key=list(corner = c(0.01, 0.98), cex=.8, col=c("gray10","gray40"),
                      text=list(lab=c(paste("slope: ", m.slope, sep=""), paste("p: ", m.prob, sep="")))),
             main=paste(gov.texts[g], " (length in ", toupper(d.name), ")", sep=""),
             xlab=paste("absolute difference in",d.name),
             ylab="proportion of shorter left conjuncts"); p
   p
}

