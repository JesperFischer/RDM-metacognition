pacman::p_load(ggplot2)


###### Simulate SDT ####################
X = 1
D = 1

sigma_e = 0.6

XD = X*D

e = seq(-4, 4, length.out = 2000)

y = dnorm(e, mean = XD, sd = sigma_e)

df = data.frame(e=e, y =y)

ggplot(df)+geom_vline(xintercept = 0, linewidth = 1.5, linetype = "dashed")+
  geom_line(aes(x= e, y = y), linewidth = 2, col = 'royalblue4') + theme_classic(base_size = 30)+
  labs(x= " ")+
  theme(axis.title.y = element_blank(),axis.text.y  = element_blank(),axis.ticks.y = element_blank())+
  geom_vline(xintercept = 0, linewidth = 1.5, linetype = "dashed")+
  scale_x_continuous(breaks = 0)+
  geom_ribbon(data = subset(df, e >= 0),aes(x = e, ymin = 0, ymax = y),
    fill = "lightskyblue",
    alpha = 0.25)+
  scale_y_continuous(expand = c(0, 0))


  
##### Simulate ranom Gaussian walk ########################

drift      <- 0.2      
noise_sd   <- 2  
n_steps    <- 500      
upper_b    <- 40      
lower_b    <- -40       

set.seed(1238)
steps = drift + rnorm(n_steps, mean = 0, sd = noise_sd)
evidence = cumsum(steps)

df = data.frame(
  t = 1:n_steps,
  evidence = evidence
)


ggplot(df, aes(t, evidence)) +
  geom_line(linewidth = 1.2, colour = "steelblue") +
  geom_hline(yintercept = upper_b, linetype = "dashed", colour = "black") +
  geom_hline(yintercept = lower_b, linetype = "dashed", colour = "black") +
  theme_classic(base_size = 30) +
  labs(x = "Time", y = " ")+
  theme(axis.title.y = element_blank(),axis.text.y  = element_blank(),axis.ticks.y = element_blank())

