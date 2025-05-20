install.packages("tidyverse")
library(tidyverse)
#1
data(mtcars)
#2 
typeof(mtcars)
mtcars
#3
car_data<- mtcars %>% mutate(Car_model=row.names(mtcars))%>% select(Car_model, everything())
car_data
#4. The column am represents the car’s transmission, 0=Automatic and 1 = Manual write a conditional statement that converts 0 and 1 to Automatic and Manual respectively. Do it within a mutate function. Move am to be the 2nd column of the dataframe. Do it in a single line of code
car_data<-car_data%>%mutate(am=if_else(am==0,'Automatic','Manual'))%>%select(Car_model,am,everything())
car_data
#5. Remove the columns vs, gear and carb from car_data using a Tidyverse function
car_data<-car_data%>%select(Car_model,am,mpg,cyl,disp,hp,drat,wt)%>%mutate(vs=NULL,gear=NULL,carb=NULL)
car_data
#6. Use ggplot to plot a box plot of mpg, with cyl as the grouping category (if you’re getting an error thing what types of values are categorical values). Add jitter to the plot, make sure the size of the dots is 2. Make sure that the x axis title is Number of cylinders and the y axis title is Miles per gallon. Brownie points for any other cosmetic modification of the plot.
car_data$cyl <- as.factor(car_data$cyl)
plot <- ggplot(data = car_data, aes(x = cyl, y = mpg)) +
  geom_boxplot()+
  geom_jitter(size = 2, aes(color = mpg))+
  xlab("Number of cylinders")+
  ylab("Miles per gallon")
plot
#7. Export the plot as pdf using a ggplot function, attach the file to your submission.
ggsave(plot = plot, filename = "plot.pdf", units = "in", width = 11.69, height = 8.27)
#8. Write a for loop that generates three dataframes, the dataframes needs to have 2 numerical columns and 1 character columns and 100 rows. The loop should store each dataframe as an element of a list and save each dataframe as a csv file. The names should be dataframe1.csv - dataframe3.csv. You can use the paste function with the loop counting variable to make sure the names are not hard coded. Hint, make the list before you loop.
listt<-vector(mode="list", length=3)

for(i in 1:3){
  df<- data.frame
  v=(1:100)
  l=(101:200)
  c <- rep(letters, length.out = 100)
  mat<-cbind(v,l,c)
  df=as.data.frame(mat)
  listt[[i]]<-df
  write.csv(df, file = paste("dataframe", i, ".csv", sep = ""), row.names = T)
}
listt
#9
imp.data <- function(Pattern= ".csv", Path = getwd(), Recursive = FALSE){
  files <- dir(path = Path, recursive = Recursive, pattern = Pattern)
  i = 1
  f_list <- vector(mode = "list", length = length(files))
  
  for (i in seq_along(files)) {
    f_list[[i]] <- read.csv(file = files[i], row.names = 1, stringsAsFactors = F) %>% mutate_if(is.numeric, log10)
    
  }
  names(f_list) <- files
  #Note return
  return(f_list)
}
f_list<- imp.data()
f_list
