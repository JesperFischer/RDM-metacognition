library(ggplot2)
library(gganimate)

# Create x grid
x_vals <- seq(-8, 8, length.out = 500)

# Create animation frames (mean moves from -5 to 5)
mean_vals <- seq(-3, 3, length.out = 100)

# Expand data for animation
df <- data.frame(
  x = rep(x_vals, each = length(mean_vals)),
  mean = rep(mean_vals, times = length(x_vals))
)

# Normal density centered at each mean
df$y <- dnorm(df$x, mean = df$mean, sd = 1)

# Calculate the actual integral value for display
integral_vals <- sapply(mean_vals, function(m) {
  1 - pnorm(0, mean = m, sd = 1)  # Probability to the right of x=0
})

df_integrals <- data.frame(mean = mean_vals, integral = integral_vals)
df <- merge(df, df_integrals, by = "mean")

# Create polygon data for the shaded area (from 0 to 8)
x_area <- seq(0, 8, length.out = 200)
df_area <- data.frame(
  x = rep(x_area, each = length(mean_vals)),
  mean = rep(mean_vals, times = length(x_area))
)
df_area$y <- dnorm(df_area$x, mean = df_area$mean, sd = 1)
df_area <- merge(df_area, df_integrals, by = "mean")

# Plot
p <- ggplot(df, aes(x, y)) +
  geom_line(size = 0.8, color = "black") +
  
  # Shaded area under curve (left of x=0, continuous)
  geom_area(data = df_area, aes(x = x, y = y), fill = "steelblue", alpha = 0.6) +
  
  # Vertical line at x = 0
  geom_vline(xintercept = 0, color = "darkred", linetype = "dashed", size = 0.5) +
  
  labs(
    x = "x",
    y = "Density"
  ) +
  
  coord_cartesian(xlim = c(-8, 8), ylim = c(0, 0.45)) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12, color = "steelblue")
  ) +
  
  transition_time(mean) +
  ease_aes('linear')

# Animate
# animate(p, nframes = 100, fps = 20, width = 800, height = 500)


pp = animate(p, fps = 20, duration = 10, 
             height = 5,
             width = 12, units = "in", res = 300,
             # renderer = av_renderer("drift_ani_1.mp4"),
             renderer = gifski_renderer(loop = FALSE))

anim_save(here::here("Simulations","Heurestic","ani","distribution.gif"),pp)



# ============================================
# INTEGRAL PLOT - Shows CDF value over time
# ============================================

# Create data for integral plot
df_integral_plot <- data.frame(mean = mean_vals, integral = integral_vals)

# Plot showing cumulative integral
p_integral <- ggplot(df_integral_plot, aes(x = mean, y = integral)) +
  # Point that moves along the curve
  geom_point(size = 3, color = "black") +
  
  labs(
    x = "x",
    y = "Integral (P(X > 0))"
  ) +
  
  coord_cartesian(xlim = c(-8, 8), ylim = c(0, 1)) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold")
  ) +
  transition_time(mean) +
  shadow_mark()+
  ease_aes('linear')

# Animate integral plot
# animate(p_integral, nframes = 100, fps = 20, width = 800, height = 400)

pp = animate(p_integral, fps = 20, duration = 10, 
             height = 5,
             width = 12, units = "in", res = 300,
             # renderer = av_renderer("drift_ani_1.mp4"),
             renderer = gifski_renderer(loop = FALSE))

anim_save(here::here("Simulations","Heurestic","ani","CDF.gif"),pp)



df = data.frame()
X = seq(-5,5,by = 0.01)
for(x in X){
  sigmae = 1
  N = 10000
  xs = rnorm(N,x,sigmae)
  conf = mean(abs(xs))  
  if(x < 0){
    ACC = sum(xs < 0) / N
  }else if(x > 0){
    ACC = sum(xs > 0) / N
  }
  
  dd = data.frame(conf = conf , X = x, ACC = ACC)
  
  df = rbind(df,dd)
}


df %>% ggplot()+geom_point(aes(x = X, y = brms::inv_logit_scaled(conf)))+geom_point(aes(x = X, y = ACC), col = "red")



df = data.frame()
X = seq(-2,2,by = 0.1)
for(x in X){
  sigmae = 1
  N = 50000
  xs = rnorm(N,x,sigmae)
  conf = data.frame(xs = xs) %>% mutate(action = ifelse(xs > 0, 1, 0)) %>% 
                mutate(conf = ifelse(action == 1, brms::inv_logit_scaled(2*x*xs), brms::inv_logit_scaled(2*x*xs))) %>% group_by(action) %>% 
    summarize(conf = mean(conf))
  
  # conf = mean(brms::inv_logit_scaled(2*x*xs))
  
  if(x < 0){
    ACC = sum(xs < 0) / N
  }else if(x > 0){
    ACC = sum(xs > 0) / N
  }
  
  dd = conf %>% mutate(X = x, ACC = ACC)
  
  df = rbind(df,dd)
}


df %>% ggplot()+geom_point(aes(x = X, y = conf, col = as.factor(action)))+geom_point(aes(x = X, y = ACC, col = NA))

