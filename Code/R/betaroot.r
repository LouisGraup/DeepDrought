# analyze root beta parameter relationships

library(tidyverse)

Y = function(beta, d) {
  # root distribution function describing
  # cumulative root fraction (Y) by depth (d, in cm)
  # from Gale & Grigal (1987)
  
  Y = 1 - beta ^ d
}

# depth in cm
d = seq(0, 200, .1)

# beta and corresponding root fraction
beta = 0.97
root_frac = Y(beta, d)

root_df = data.frame(depth=d, root_frac)

ggplot(root_df, aes(root_frac, depth))+geom_line()+
  scale_x_reverse()+scale_y_reverse()+theme_bw()


betas = c(0.9, 0.95, 0.975, 0.99, 0.995)
roots = sapply(betas, Y, d)

roots_df = data.frame(d, roots)
colnames(roots_df) = c("depth", betas)

roots_long = roots_df %>% pivot_longer(-depth, names_to="beta", values_to="root_frac")
roots_long$beta = as.numeric(roots_long$beta)

ggplot(roots_long, aes(root_frac, depth, group=beta, color=beta))+geom_line()+
  scale_x_reverse()+scale_y_reverse()+theme_bw()


