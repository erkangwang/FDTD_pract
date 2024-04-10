from FDTD_cell import fdtd_cell_1D
import numpy as np
from numpy.random import random
import matplotlib.pyplot as plt
import time

def print_hi(name):
    # Use a breakpoint in the code line below to debug your script.
    print(f'Hi, {name}')  # Press Ctrl+F8 to toggle the breakpoint.


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    print_hi('FDTD_practice')

x = np.array([[1, 2, 3], [4, 5, 6]])
print(x)


##creat class and object
#from turtle import Turtle, Screen
#timmy = Turtle()
#print(timmy)
#my_screen = Screen()
#print(my_screen.canvheight)
#my_screen.exitonclick()

# np array manipulation
my_array = np.array([[1, 2, 3, 4],
                     [2, 3, 4, 5]])
print(f'array has {my_array.shape[0]} rows and {my_array.shape[1]} columns')

a = np.arange(10, 30, 2)
print(a[-3:])
print(np.flip(a))

z = random((3, 3, 3))
print(z.shape)

num_time_point=601
t_space_ind = np.linspace(0, num_time_point-1, num=num_time_point, dtype=int)
src_ind = 50
#E_amp = np.linspace(0, 100, num=100)
aa = fdtd_cell_1D(t_space_ind, src_ind)
aa = aa.reshape((602, 101))
print(aa.shape)

x = aa[0]
y = aa[1]
# to run GUI event loop
plt.ion()
# here we are creating sub plots
figure, ax = plt.subplots(figsize=(10, 5))
plt.ylim(-6, 6)
#plt.xlim(0,10)
line1, = ax.plot(x, y)
# setting title
plt.title("1D FDTD", fontsize=20)
# setting x-axis label and y-axis label
plt.xlabel("X-axis")
plt.ylabel("Y-axis")

for t in t_space_ind:

    # creating new Y values
    new_y = aa[t]
    # updating data values
    line1.set_xdata(x)
    line1.set_ydata(new_y)
    # drawing updated values
    figure.canvas.draw()

    # This will run the GUI event
    # loop until all UI events
    # currently waiting have been processed
    figure.canvas.flush_events()

    time.sleep(0.01)
#plt.show()




































