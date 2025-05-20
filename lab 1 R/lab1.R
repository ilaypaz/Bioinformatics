library(stringr)
#q1
b="How much wood would a woodchuck chuck if a woodchuck could chuck wood? He would chuck, he would, as much as he could, and chuck as much wood As a woodchuck would if a woodchuck could chuck wood"
b2 = str_split(b, " ")
length(b)
str_count(b, "(?i)wood")
#Q2
bases=c("A","T","G","C")
set.seed(12)
dna =paste(sample(bases,10000, replace = T), collapse = "")
dna3=dna2=dna1=dna
dna1==dna2
dna2==dna3
dna3==dna1
((str_count(dna3, "(?i)[GC]"))/10000)*100
((str_count(dna1, "(?i)[GC]"))/10000)*100
((str_count(dna2, "(?i)[GC]"))/10000)*100
#q3
list <- vector(mode = "list", length = 4)
list
#q4
names(list) <- c("E1", "E2", "E3", "E4")
names(list)
#q5
v=(1:10)
l=(30:39)
c=(2:11)
mat <-cbind(v,l,c)
list[[1]] = mat
list
#q6
list[[2]] = log10(mat)
df = as.data.frame(mat)
list
#q7
list[[3]] = df
df
list
df1 = as.data.frame(mat)
#q8
df1$New_column <- letters[1:10]
list[[4]] <- df1
list
#q9
list[[4]] <- as.data.frame(lapply(df1, function(x) {
  if (is.numeric(x)) {
    log10(x)  
  } else {
    x         
  }
}))


list
