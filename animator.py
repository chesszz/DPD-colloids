"""
Animation of Elastic collisions with Gravity

author: Jake Vanderplas
email: vanderplas@astro.washington.edu
website: http://jakevdp.github.com
license: BSD
Please feel free to use and modify this, but keep the above information. Thanks!
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

ANIMATE         = False
TILE            = False
PLOT_TEMP       = True
PLOT_MAXBOLTZ   = True

#------------------------------------------------------------
### READ FROM DYN FILE ###

# Read dyn file to get the simulation results
with open("dyn.out","r") as f:
    dyn_data = f.readlines()

# Chop off first line which is the column header
dyn_data = dyn_data[1:]

# Split each line by commas, and then convert each entry from str->float
for index, line in enumerate(dyn_data):
    dyn_data[index] = line.split(",")
    for inner_index, entry in enumerate(dyn_data[index]):
        dyn_data[index][inner_index] = float(entry)

# Convert to np array to use multiple indexing syntax
dyn_data = np.array(dyn_data)

#------------------------------------------------------------
### READ FROM INPUT FILE ###

NUM_INPUTS = 9

with open("inputs.in", "r") as g:
    inputs = g.readlines()
# Only have 9 inputs
input_data = inputs[:NUM_INPUTS]
input_data = list(map(float, input_data))

# Skip 1 line because it's a space
labels = inputs[NUM_INPUTS+1:]
# Chop out the part before the first space
for index, string in enumerate(labels):
    string = string.strip()
    bef,space,aft = string.partition(" ")
    labels[index] = aft

input_dict = dict(zip(labels,input_data))

N_WATER = int(input_dict["N_WATER"])
BOX_SIZE = input_dict["BOX_SIZE"]
N_STEPS = int(input_dict["N_STEPS"])
TIME_STEP = input_dict["TIME_STEP"]

#------------------------------------------------------------
### READ FROM TEMP FILE ###

# Read temp file to get the temp results
with open("temp.out","r") as t:
    temp_data = t.readlines()

# Split by space and grab the last part, which is the temp, and then convert to float.
temp = [float(temp_entry.split()[-1]) for temp_entry in temp_data]
times = np.arange(0, N_STEPS * TIME_STEP, TIME_STEP)

#------------------------------------------------------------

if ANIMATE:
    # set up figure and animation
    fig = plt.figure()
    #fig.subplots_adjust(left=0, right=1, bottom=0, top=1)

    if TILE:
        limits = (-1.2 * BOX_SIZE, 2.2 * BOX_SIZE)
    else:
        limits = (-0.2 * BOX_SIZE, 1.2 * BOX_SIZE)

    ax = fig.add_subplot(111, aspect='equal', autoscale_on=False, xlim=limits, ylim=limits)

    # particles holds the locations of the particles
    particles, = ax.plot([], [], 'bo', ms=6)

    # rect is the box edge
    rect = plt.Rectangle((0, 0), BOX_SIZE, BOX_SIZE, ec='none', lw=2, fc='none')
    ax.add_patch(rect)

    # Convert limits to array to do scalar multiplication, then unpack into the position args for text
    text_pos = 0.95 * np.array(limits)
    time_text = ax.text(*text_pos, '', fontsize=12)

    def init():
        """initialize animation"""
        global rect
        particles.set_data([], [])
        rect.set_edgecolor('none')
        time_text.set_text('')
        return particles, rect, time_text

    def animate(i):
        """perform animation step"""
        global rect, ax, fig, dyn_data

        # ms = int(fig.dpi * 2 * box.size * fig.get_figwidth()
                 # / np.diff(ax.get_xbound())[0])
        
        # update pieces of the animation
        rect.set_edgecolor('k')

        if TILE: 
            ## PLOTS TILINGS OF THE MAIN CELL
            particles.set_data(
                list(dyn_data[N_WATER*i:N_WATER*i+N_WATER, 3]) +
                list(dyn_data[N_WATER*i:N_WATER*i+N_WATER, 3]) +
                list(dyn_data[N_WATER*i:N_WATER*i+N_WATER, 3]) +
                [a-BOX_SIZE for a in dyn_data[N_WATER*i:N_WATER*i+N_WATER, 3]] +
                [a-BOX_SIZE for a in dyn_data[N_WATER*i:N_WATER*i+N_WATER, 3]] +
                [a-BOX_SIZE for a in dyn_data[N_WATER*i:N_WATER*i+N_WATER, 3]] +
                [a+BOX_SIZE for a in dyn_data[N_WATER*i:N_WATER*i+N_WATER, 3]] +
                [a+BOX_SIZE for a in dyn_data[N_WATER*i:N_WATER*i+N_WATER, 3]] +
                [a+BOX_SIZE for a in dyn_data[N_WATER*i:N_WATER*i+N_WATER, 3]] ,
                
                list(dyn_data[N_WATER*i:N_WATER*i+N_WATER, 4])                 +
                [a+BOX_SIZE for a in dyn_data[N_WATER*i:N_WATER*i+N_WATER, 4]] +
                [a-BOX_SIZE for a in dyn_data[N_WATER*i:N_WATER*i+N_WATER, 4]] +
                list(dyn_data[N_WATER*i:N_WATER*i+N_WATER, 4])                 +
                [a+BOX_SIZE for a in dyn_data[N_WATER*i:N_WATER*i+N_WATER, 4]] +
                [a-BOX_SIZE for a in dyn_data[N_WATER*i:N_WATER*i+N_WATER, 4]] +
                list(dyn_data[N_WATER*i:N_WATER*i+N_WATER, 4])                 +
                [a+BOX_SIZE for a in dyn_data[N_WATER*i:N_WATER*i+N_WATER, 4]] +
                [a-BOX_SIZE for a in dyn_data[N_WATER*i:N_WATER*i+N_WATER, 4]] 
            )
        else:
            ## ONLY PLOTS THE MAIN CELL

            # 3 and 4 represent the 4th and 5th column of dyn_data - x and y coords
            # First part is a window that jumps by N_WATER rows every timestep
            # and has width of N_WATER, to include all particles at a particular
            # timestep value.
            particles.set_data( dyn_data[N_WATER*i:N_WATER*i+N_WATER, 3],
                                dyn_data[N_WATER*i:N_WATER*i+N_WATER, 4])

        if TILE:
            markersize = 82/BOX_SIZE
        else:
            markersize = 200/BOX_SIZE

        particles.set_markersize(markersize)
        time_text.set_text('Time = {0} / {1}'.format(i, N_STEPS))
        
        return particles, rect, time_text

    ani = animation.FuncAnimation(fig, animate, frames=N_STEPS,
                                  interval=1, blit=True, init_func=init)


    # save the animation as an mp4.  This requires ffmpeg or mencoder to be
    # installed.  The extra_args ensure that the x264 codec is used, so that
    # the video can be embedded in html5.  You may need to adjust this for
    # your system: for more information, see
    # http://matplotlib.sourceforge.net/api/animation_api.html
    # ani.save('particle_box.mp4', fps=30, extra_args=['-vcodec', 'libx264'])

if PLOT_TEMP:
    plt.figure()
    ax = plt.subplot()
    ax.set_ylim([0, 3])

    ax.plot(times, temp)
    ax.plot(times, np.ones(len(times)))

    half_temp = temp[len(temp)//2:]
    avg_half_temp = sum(half_temp) / len(half_temp)
    ax.text(0, 0, 'Average later half temp: {0}'.format(avg_half_temp), fontsize=12)
    print("Average temperature over the later half of the data was {0}.".format(avg_half_temp))

if PLOT_MAXBOLTZ:
    plt.figure()
    # Take mod of velocity vectors
    v_list = np.sqrt(dyn_data[-N_WATER:, 6]**2+ dyn_data[-N_WATER:, 7]**2 + dyn_data[-N_WATER:, 8]**2)

    # Initial velocity list
    # v_list = np.sqrt(dyn_data[:N_WATER, 6]**2+ dyn_data[:N_WATER, 7]**2 + dyn_data[:N_WATER, 8]**2)

    # Have at most 50 bins, at least 5 bins, otherwise N_WATER / 8 bins
    num_bins = min(50, max(5, N_WATER//8))
    plt.hist(v_list, bins=num_bins, normed=True)

    x = np.linspace(0, 1.2 * max(v_list), 100)
    
    # Calculate the theoretical velocity
    y = 4 * np.pi * x**2 * (1 / (2 * np.pi))**1.5 * np.exp(-1 * x**2 / 2)

    plt.plot(x,y)
    
plt.show()
plt.close()
