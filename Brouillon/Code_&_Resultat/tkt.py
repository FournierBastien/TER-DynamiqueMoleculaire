# coding: utf-8
 
from tkinter import * 

fenetre = Tk()

value = DoubleVar()
scale = Scale(fenetre, variable=value)
scale.pack()
fenetre.mainloop()
