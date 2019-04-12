# The aim is to write the code in python first then trasncribe to C++
# This model essentially counts the charge on the inner sphere due to protons leaving a positive charge and holes leaving to be filled by electrons on the outer sphere
# make the outer circle radius 4 and the inner circle radius 1

import numpy as np
import tkinter as tk
import random

root = tk.Tk()
myCanvas = tk.Canvas(root, width=1500, height=1500,  borderwidth=0, highlightthickness=0, bg="white")
myCanvas.pack()

def create_circle(x, y, r, canvasName): #center coordinates, radius
    x0 = x - r
    y0 = y - r
    x1 = x + r
    y1 = y + r
    return canvasName.create_oval(x0, y0, x1, y1)

create_circle(750, 750, 250, myCanvas)
create_circle(750, 750, 100, myCanvas)

# Due to how graphics draws, the origin is at (750,750)

class proton:
    def __init__(self):
        # intialise with some random radius that makes the proton outside the outer sphere and have some randome angle phi. We can then convert them to cartesian.
        self.r = round(random.uniform(250.01,750), 3)
        self.phi = round(random.uniform(0,2*np.pi),3)
        self.vx = round(np.random.normal(0,500),3)
        self.vy = round(np.random.normal(-200,100),3)
        self.x = round(self.r * np.cos(self.phi) + 750,3)
        self.y = round(self.r * np.sin(self.phi) + 750,3)
    def project(self):
        t = 15
        self.x2 = self.x + self.vx * t
        self.y2 = self.y + self.vy * t   

#### Apply all function to apply same method to a list of objects ####
def apply_on_all(seq):
    for obj in seq:
        obj.project()
######################################################################

instancelist = [ proton() for i in range(29)]
apply_on_all(instancelist)
for obj in instancelist:
    myCanvas.create_line(obj.x,obj.y,obj.x2,obj.y2, fill = "blue")
root.mainloop()

