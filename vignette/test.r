devtools::load_all("C:/repository/concentrationgames")
devtools::document("C:/repository/concentrationgames")

xs = define_xs(0, 10, 5, 50, 8, 50)
excess.mass = 20
H = 15
b = 0.75

res1 = distribute_mass(excess.mass, xs, H, b)
with(res1, print(c(LOB.mass, ROB.mass, C.mass)))

us = 0.2
k = 0.4
w = 0.03

res2 = distribute_mass(excess.mass, xs, H, b, rouse, us = 0.2, k = 0.4, w = 0.03)
with(res2, print(c(LOB.mass, ROB.mass, C.mass)))


library(ggplot2)
ggplot(mass_to_bedchange(deposit_mass(res2, 0.5, 0.01))) + 
  aes(x = y, y = z + dz*1e4) + geom_line(linetype = "dashed") +
  geom_line(aes(y = z))


ggplot(mass_to_bedchange(deposit_mass(res1, 0.5, 0.05, 30))) + 
  aes(x = y, y = z + dz*1e4) + geom_line(linetype = "dashed") +
  geom_line(aes(y = z))



