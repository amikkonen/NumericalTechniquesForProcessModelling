# -*- coding: utf-8 -*-
"""
Simple real time plotted and logger. Read data from serial.

Serial data is provided from Arduino based infrared temperature sensor. See 
    verticalPlateDemo.ino

Easy to modify for other similar data plotting from serial.

Antti Mikkonen, a.mikkonen@iki.fi, winter 2019.

"""

import serial
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time
import os


def animate(i, ser, ax, T_reds, T_infs, times, first_time):

    # Read data from serial, perform some string formating
    parts = ser.readline().decode("utf-8").replace("\n","").split(",")
    # Infrared temperature
    T_red = float(parts[0])
    # Sensor temperature
    T_inf = float(parts[1])
    # Append to storage
    T_reds.append(T_red) 
    T_infs.append(T_inf)
    times.append((time.time()-first_time)/60) # minutes
    
    # Print
    print(round(times[-1],2), T_red, T_inf)    
    
    # Write to log
    with open("log", "a") as ofile:
        ofile.write(str(times[-1]) + ", " +  str(T_red) + ", " +  str(T_inf) + "\n")
    
    # Plotting
    ax.clear()
    ax.plot(times,T_reds,"r", label="$T$")
    ax.plot(times,T_infs,"b", label="$T_{\infty}$")
    ax.legend(frameon=False)
    ax.set_xlabel("$t$ $(\mathrm{min})$")
    ax.set_ylabel("$T$ $(\mathrm{^\circ C})$")
    ax.set_xlim(times[0], None)
    ax.grid(True)
    
def main():
    
    # Initialize serial communication
    # This setup for Ubuntu and Arduino Uno
    # Likely different path on other hardware/operating system
    ser = serial.Serial('/dev/ttyACM3', 9600)

    # Initialize  plots and storages
    fig, ax = plt.subplots()
    T_reds = []
    T_infs = []    
    times = []
    first_time = time.time()

    # Clean old data from file
    if os.path.exists("log"):
        os.remove("log")
    
    # Animate results and write to log
    ani = animation.FuncAnimation(fig, animate, 
                                  fargs=(ser, ax, T_reds, T_infs, 
                                         times, first_time), 
                                  interval=1)
    plt.show()
    
if __name__ == "__main__":
    main()
