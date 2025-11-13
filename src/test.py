import pygame
import numpy as np

pygame.init()
WIDTH, HEIGHT = 800, 800
WIN = pygame.display.set_mode((WIDTH, HEIGHT))
pygame.display.set_caption("Planet Simulation")

white = (255,255,255)
class Planet:
    AU = 149.6e6*1000
    G = 6.67428e-11
    SCALE = 250/AU # 1AU = 100px
    TIMESTEP = 3600*24
    def __init__(self, x, y, radius, color, mass):
        self.x = x
        self.y = y
        self.radius = radius
        self.color = color
        self.mass = mass

        self.orbit = []
        self.sun = False
        self.distance_to_sun = 0

        self.x_vel = 0
        self.y_vel = 0

    def draw(self, win):
        x = self.x*self.SCALE+WIDTH/2
        y = self.y*self.SCALE+HEIGHT/2
        pygame.draw.circle(win, self.color, (x, y), self.radius)





def main():
    run = True
    clock = pygame.time.Clock()
    sun = Planet(0,0,30,(255,255,255),(255,255,255))


    while run:
        clock.tick(60)
        WIN.fill(white)
        pygame.display.update()
        for event in pygame.event.get(): # get all keypresses etc
            if event.type == pygame.QUIT:
                run = False
    pygame.quit()

main()

