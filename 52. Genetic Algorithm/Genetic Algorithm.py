from tkinter import *
from tkinter import ttk
import string, random

letters = string.ascii_lowercase + ' ' + string.ascii_uppercase


def update_label():
    global sample01, sample
    sample.set(''.join(random.choice(letters) for i in range(20)))



# Main window 
root = Tk()
root.geometry('980x760')
root.title('Genetic Algorithm')

# Definition of the frame inside the window
# and the grid to localize the elements
frm = ttk.Frame(root, padding=10)
frm.grid()


title = ttk.Label(frm, text="Genetic Algorithm\n", font=('Arial', 40))
title.grid(column=0, row=0)

Best_result = ttk.Label(frm, text="The Best result is", font=('Arial', 20)).grid(column=0, row=2, sticky=W)
Best = ttk.Label(frm, text="XFOJWOSDFNAÑSKND", font=('Arial', 30), anchor=NW).grid(column=0, row=3)

# Horizontal space
Label(frm, text='', font=('Arial', 20)).grid(column=1, row=2, padx=150)


# Samples visualization
sample = StringVar()
sample.set(''.join(random.choice(letters) for i in range(20)))

sample01 = ttk.Label(frm, textvariable=sample, font=('Arial', 20)).grid(column=2, row=2)
sample02 = ttk.Label(frm, text=sample, font=('Arial', 20)).grid(column=2, row=3)
sample03 = ttk.Label(frm, text=sample, font=('Arial', 20)).grid(column=2, row=4)
sample04 = ttk.Label(frm, text=sample, font=('Arial', 20)).grid(column=2, row=5)
sample05 = ttk.Label(frm, text=sample, font=('Arial', 20)).grid(column=2, row=6)
sample06 = ttk.Label(frm, text=sample, font=('Arial', 20)).grid(column=2, row=7)
sample07 = ttk.Label(frm, text=sample, font=('Arial', 20)).grid(column=2, row=8)
sample08 = ttk.Label(frm, text=sample, font=('Arial', 20)).grid(column=2, row=9)
sample09 = ttk.Label(frm, text=sample, font=('Arial', 20)).grid(column=2, row=10)
sample10 = ttk.Label(frm, text=sample, font=('Arial', 20)).grid(column=2, row=11)
sample11 = ttk.Label(frm, text=sample, font=('Arial', 20)).grid(column=2, row=12)
sample12 = ttk.Label(frm, text=sample, font=('Arial', 20)).grid(column=2, row=13)
sample13 = ttk.Label(frm, text=sample, font=('Arial', 20)).grid(column=2, row=14)
sample14 = ttk.Label(frm, text=sample, font=('Arial', 20)).grid(column=2, row=15)
sample15 = ttk.Label(frm, text=sample, font=('Arial', 20)).grid(column=2, row=16)
sample16 = ttk.Label(frm, text=sample, font=('Arial', 20)).grid(column=2, row=17)
sample17 = ttk.Label(frm, text=sample, font=('Arial', 20)).grid(column=2, row=18)
sample18 = ttk.Label(frm, text=sample, font=('Arial', 20)).grid(column=2, row=19)
sample19 = ttk.Label(frm, text=sample, font=('Arial', 20)).grid(column=2, row=20)
sample20 = ttk.Label(frm, text=sample, font=('Arial', 20)).grid(column=2, row=21)
sample21 = ttk.Label(frm, text=sample, font=('Arial', 20)).grid(column=2, row=22)
# Vertical space
ttk.Label(frm, text='', font=('Arial', 20)).grid(column=2, row=23)


quit_button = ttk.Button(frm, text="Quit", command=root.destroy)
quit_button.grid(column=2, row=24)
start_button = ttk.Button(frm, text="Start", command=update_label)
start_button.grid(column=0, row=24)

# ttk Style
s = ttk.Style()
s.configure('TButton', foreground='black')
s.configure('TLabel', foreground='black', anchor='NW')

# Main Loop for the window
root.mainloop()