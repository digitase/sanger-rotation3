
library(ggplot2)
library(reshape2)

# n hps in single
# Counter({'st22': 1285, 'st239': 237, 'st8': 1186})
# n in 2
# Counter({'st22': 226, 'st239': 129, 'st8': 287})
# n in 3
# Counter({'st22': 44, 'st239': 44, 'st8': 44})
df <- data.frame(
    x=c('ST22', 'ST8', 'ST239'),
    st3=c(44, 44, 44),
    st2=c(226, 287, 129),
    st1=c(1285, 1186, 237)
)

g <- ggplot(melt(df)) +
    geom_col(aes(x, value, fill=variable), position = position_stack(reverse = TRUE)) + 
    coord_flip()
print(g)

