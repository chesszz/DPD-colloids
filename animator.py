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
PLOT_CROSSFLOW  = True

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
V_SHEAR = input_dict["V_SHEAR"]

#------------------------------------------------------------
### READ FROM DYN FILE ###

# Read the entire dyn file to get the simulation results for all time
if ANIMATE:
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

# Only read the last time step of dyn results
if PLOT_MAXBOLTZ or PLOT_CROSSFLOW:

    # Add 1 to account for the column header row
    num_lines = N_WATER * N_STEPS + 1
    # Grab only the last N_WATER steps
    with open("dyn.out","r") as f:
        dyn_data_last = [line for index, line in enumerate(f) if index > (num_lines - N_WATER - 1)]

    # Split each line by commas, and then convert each entry from str->float
    for index, line in enumerate(dyn_data_last):
        dyn_data_last[index] = line.split(",")
        for inner_index, entry in enumerate(dyn_data_last[index]):
            dyn_data_last[index][inner_index] = float(entry)

    # Convert to np array to use multiple indexing syntax
    dyn_data_last = np.array(dyn_data_last)
    

#------------------------------------------------------------
### READ FROM TEMP FILE ###

if PLOT_TEMP:
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
        #limits = (-1.2 * BOX_SIZE, 2.2 * BOX_SIZE)
        limits = (-BOX_SIZE, 2 * BOX_SIZE)
    else:
        limits = (-0.2 * BOX_SIZE, 1.2 * BOX_SIZE)

    ax_animate = fig.add_subplot(111, aspect='equal', autoscale_on=False, xlim=limits, ylim=limits)

    # particles holds the locations of the particles
    particles, = ax_animate.plot([], [], 'bo', ms=6)
    centres, = ax_animate.plot([], [], 'ko', ms=2)

    # rect is the box edge
    rect = plt.Rectangle((0, 0), BOX_SIZE, BOX_SIZE, ec='none', lw=2, fc='none')
    ax_animate.add_patch(rect)

    if TILE:
        rect_top = plt.Rectangle((0, BOX_SIZE), BOX_SIZE, BOX_SIZE, ec='none', lw=2, fc='none')
        ax_animate.add_patch(rect_top)

        rect_tl = plt.Rectangle((-BOX_SIZE, BOX_SIZE), BOX_SIZE, BOX_SIZE, ec='none', lw=2, fc='none')
        ax_animate.add_patch(rect_tl)

        rect_tl2 = plt.Rectangle((-2*BOX_SIZE, BOX_SIZE), BOX_SIZE, BOX_SIZE, ec='none', lw=2, fc='none')
        ax_animate.add_patch(rect_tl2)

        rect_tr = plt.Rectangle((BOX_SIZE, BOX_SIZE), BOX_SIZE, BOX_SIZE, ec='none', lw=2, fc='none')
        ax_animate.add_patch(rect_tr)

        rect_bot = plt.Rectangle((0, -BOX_SIZE), BOX_SIZE, BOX_SIZE, ec='none', lw=2, fc='none')
        ax_animate.add_patch(rect_bot)

        rect_bl = plt.Rectangle((-BOX_SIZE, -BOX_SIZE), BOX_SIZE, BOX_SIZE, ec='none', lw=2, fc='none')
        ax_animate.add_patch(rect_bl)

        rect_br = plt.Rectangle((BOX_SIZE, -BOX_SIZE), BOX_SIZE, BOX_SIZE, ec='none', lw=2, fc='none')
        ax_animate.add_patch(rect_br)

        rect_bl2 = plt.Rectangle((-BOX_SIZE, -BOX_SIZE), BOX_SIZE, BOX_SIZE, ec='none', lw=2, fc='none')
        ax_animate.add_patch(rect_bl2)

    # Convert limits to array to do scalar multiplication, then unpack into the position args for text
    text_pos = 0.95 * np.array(limits)
    time_text = ax_animate.text(*text_pos, '', fontsize=12)

    def init():
        """initialize animation"""
        return particles, centres, rect, rect_top, rect_tl, rect_tl2, rect_tr, rect_bot, rect_bl, rect_br, rect_bl2, time_text

    def animate(i):
        """perform animation step"""
        global rect, rect_top, ax_animate, fig, dyn_data

        # ms = int(fig.dpi * 2 * box.size * fig.get_figwidth()
                 # / np.diff(ax_animate.get_xbound())[0])
        
        # update pieces of the animation
        rect.set_edgecolor('k')

        xmin, xmax = ax_animate.get_xlim()
        range_x = xmax - xmin

        if TILE: 

            t = i * TIME_STEP

            rect_top.set_edgecolor('k')
            rect_tl.set_edgecolor('k')
            rect_tl2.set_edgecolor('k')
            rect_tr.set_edgecolor('k')
            rect_bot.set_edgecolor('k')
            rect_bl.set_edgecolor('k')
            rect_bl2.set_edgecolor('k')
            rect_br.set_edgecolor('k')

            rect_top.set_x((            t * V_SHEAR + BOX_SIZE) % (3*BOX_SIZE) - BOX_SIZE)
            rect_bot.set_x((          - t * V_SHEAR + BOX_SIZE) % (3*BOX_SIZE) - BOX_SIZE)
            rect_tl.set_x( (-BOX_SIZE + t * V_SHEAR + BOX_SIZE) % (3*BOX_SIZE) - BOX_SIZE)
            rect_tl2.set_x((            t * V_SHEAR)           % BOX_SIZE      - 2*BOX_SIZE)
            rect_bl.set_x( (-BOX_SIZE - t * V_SHEAR + BOX_SIZE) % (3*BOX_SIZE) - BOX_SIZE)
            rect_bl2.set_x((          - t * V_SHEAR           ) % BOX_SIZE     - 2*BOX_SIZE)
            rect_tr.set_x( ( BOX_SIZE + t * V_SHEAR + BOX_SIZE) % (3*BOX_SIZE) - BOX_SIZE)
            rect_br.set_x( ( BOX_SIZE - t * V_SHEAR + BOX_SIZE) % (3*BOX_SIZE) - BOX_SIZE)
            
            ## PLOTS TILINGS OF THE MAIN CELL
            particles.set_data(
                list(dyn_data[N_WATER*i:N_WATER*i+N_WATER, 3]) +                                                # mid
                [(a+t*V_SHEAR-xmin)%range_x+xmin          for a in dyn_data[N_WATER*i:N_WATER*i+N_WATER, 3]] +  # top                                       
                [(a-t*V_SHEAR-xmin)%range_x+xmin          for a in dyn_data[N_WATER*i:N_WATER*i+N_WATER, 3]] +  # bot
                [a-BOX_SIZE                               for a in dyn_data[N_WATER*i:N_WATER*i+N_WATER, 3]] +  # left
                [(a-BOX_SIZE+t*V_SHEAR-xmin)%range_x+xmin for a in dyn_data[N_WATER*i:N_WATER*i+N_WATER, 3]] +  # top left
                [(a-BOX_SIZE-t*V_SHEAR-xmin)%range_x+xmin for a in dyn_data[N_WATER*i:N_WATER*i+N_WATER, 3]] +  # bot left
                [a+BOX_SIZE                               for a in dyn_data[N_WATER*i:N_WATER*i+N_WATER, 3]] +  # right
                [(a+BOX_SIZE+t*V_SHEAR-xmin)%range_x+xmin for a in dyn_data[N_WATER*i:N_WATER*i+N_WATER, 3]] +  # top right
                [(a+BOX_SIZE-t*V_SHEAR-xmin)%range_x+xmin for a in dyn_data[N_WATER*i:N_WATER*i+N_WATER, 3]] ,  # bot right
                
                list(dyn_data[N_WATER*i:N_WATER*i+N_WATER, 4])                 +                                # mid
                [a+BOX_SIZE for a in dyn_data[N_WATER*i:N_WATER*i+N_WATER, 4]] +                                # top
                [a-BOX_SIZE for a in dyn_data[N_WATER*i:N_WATER*i+N_WATER, 4]] +                                # bot
                list(dyn_data[N_WATER*i:N_WATER*i+N_WATER, 4])                 +                                # left
                [a+BOX_SIZE for a in dyn_data[N_WATER*i:N_WATER*i+N_WATER, 4]] +                                # top left
                [a-BOX_SIZE for a in dyn_data[N_WATER*i:N_WATER*i+N_WATER, 4]] +                                # bot left
                list(dyn_data[N_WATER*i:N_WATER*i+N_WATER, 4])                 +                                # right
                [a+BOX_SIZE for a in dyn_data[N_WATER*i:N_WATER*i+N_WATER, 4]] +                                # top right
                [a-BOX_SIZE for a in dyn_data[N_WATER*i:N_WATER*i+N_WATER, 4]]                                  # bot right
            )

            centres.set_data(
                list(dyn_data[N_WATER*i:N_WATER*i+N_WATER, 3]) +                                                # mid
                [(a+t*V_SHEAR-xmin)%range_x+xmin          for a in dyn_data[N_WATER*i:N_WATER*i+N_WATER, 3]] +  # top                                       
                [(a-t*V_SHEAR-xmin)%range_x+xmin          for a in dyn_data[N_WATER*i:N_WATER*i+N_WATER, 3]] +  # bot
                [a-BOX_SIZE                               for a in dyn_data[N_WATER*i:N_WATER*i+N_WATER, 3]] +  # left
                [(a-BOX_SIZE+t*V_SHEAR-xmin)%range_x+xmin for a in dyn_data[N_WATER*i:N_WATER*i+N_WATER, 3]] +  # top left
                [(a-BOX_SIZE-t*V_SHEAR-xmin)%range_x+xmin for a in dyn_data[N_WATER*i:N_WATER*i+N_WATER, 3]] +  # bot left
                [a+BOX_SIZE                               for a in dyn_data[N_WATER*i:N_WATER*i+N_WATER, 3]] +  # right
                [(a+BOX_SIZE+t*V_SHEAR-xmin)%range_x+xmin for a in dyn_data[N_WATER*i:N_WATER*i+N_WATER, 3]] +  # top right
                [(a+BOX_SIZE-t*V_SHEAR-xmin)%range_x+xmin for a in dyn_data[N_WATER*i:N_WATER*i+N_WATER, 3]] ,  # bot right
                
                list(dyn_data[N_WATER*i:N_WATER*i+N_WATER, 4])                 +                                # mid
                [a+BOX_SIZE for a in dyn_data[N_WATER*i:N_WATER*i+N_WATER, 4]] +                                # top
                [a-BOX_SIZE for a in dyn_data[N_WATER*i:N_WATER*i+N_WATER, 4]] +                                # bot
                list(dyn_data[N_WATER*i:N_WATER*i+N_WATER, 4])                 +                                # left
                [a+BOX_SIZE for a in dyn_data[N_WATER*i:N_WATER*i+N_WATER, 4]] +                                # top left
                [a-BOX_SIZE for a in dyn_data[N_WATER*i:N_WATER*i+N_WATER, 4]] +                                # bot left
                list(dyn_data[N_WATER*i:N_WATER*i+N_WATER, 4])                 +                                # right
                [a+BOX_SIZE for a in dyn_data[N_WATER*i:N_WATER*i+N_WATER, 4]] +                                # top right
                [a-BOX_SIZE for a in dyn_data[N_WATER*i:N_WATER*i+N_WATER, 4]]                                  # bot right
            )
        else:
            ## ONLY PLOTS THE MAIN CELL

            # 3 and 4 represent the 4th and 5th column of dyn_data - x and y coords
            # First part is a window that jumps by N_WATER rows every timestep
            # and has width of N_WATER, to include all particles at a particular
            # timestep value.
            particles.set_data( dyn_data[N_WATER*i:N_WATER*i+N_WATER, 3],
                                dyn_data[N_WATER*i:N_WATER*i+N_WATER, 4])
            centres.set_data( dyn_data[N_WATER*i:N_WATER*i+N_WATER, 3],
                                dyn_data[N_WATER*i:N_WATER*i+N_WATER, 4])

        if TILE:
            markersize = 82/BOX_SIZE
        else:
            markersize = 260/BOX_SIZE

        particles.set_markersize(markersize)
        time_text.set_text('Time = {0} / {1}'.format(i, N_STEPS))
        
        return particles, centres, rect, rect_top, rect_tl, rect_tl2, rect_tr, rect_bot, rect_bl, rect_br, rect_bl2, time_text

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

    quart_temp = temp[3*len(temp)//4:]
    avg_quart_temp = sum(quart_temp) / len(quart_temp)
    ax.text(0, 0, 'Average last quarter temp: {0:.5}'.format(avg_quart_temp), fontsize=12)

if PLOT_MAXBOLTZ:
    plt.figure()
    # Take mod of velocity vectors
    v_list = np.sqrt(dyn_data_last[:, 6]**2+ dyn_data_last[:, 7]**2 + dyn_data_last[:, 8]**2)

    # Have at most 50 bins, at least 5 bins, otherwise N_WATER / 8 bins
    num_bins = min(50, max(5, N_WATER//8))
    plt.hist(v_list, bins=num_bins, normed=True)

    x = np.linspace(0, 1.2 * max(v_list), 100)
    
    # Calculate the theoretical velocity
    y = 4 * np.pi * x**2 * (1 / (2 * np.pi))**1.5 * np.exp(-1 * x**2 / 2)

    plt.plot(x,y)

if PLOT_CROSSFLOW:
    num_boundaries = 21
    num_bins = num_boundaries - 1
    y_bins = np.linspace(0, BOX_SIZE, num_boundaries)
    y_points = np.linspace(0, BOX_SIZE, num_bins)
    
    # Create an array with the y-pos and the x-vel (4 and 6 positions)
    y_vx_list = dyn_data_last[:, [4,6]]
    # Takes the first index of each row, sort them, and then index the original
    # array using that order. In other words, sorts the rows by their y-pos.
    y_vx_list = y_vx_list[y_vx_list[:,0].argsort()]

    # Tells us which bin each row belongs in
    bin_indices = np.searchsorted(y_bins, y_vx_list[:,0]) - 1

    xvel_entries = np.zeros([num_bins, 2])
    xvel_avg = np.zeros(num_bins)

    # Sum up the y-vel for each bin by referring to the bin_indices list
    # y_bins_entries now has a list of total y-vel in that bin + number of particles in that bin
    for index, row in enumerate(y_vx_list):
        bin_belonging = bin_indices[index]
        xvel_entries[bin_belonging, 0] += row[1]
        xvel_entries[bin_belonging, 1] += 1

    count_list = xvel_entries[:, 1]

    # Compute the average y-vel for each bin
    for index, row in enumerate(xvel_entries):
        if row[1] == 0:
            xvel_avg[index] = 0
        else:
            xvel_avg[index] = row[0] / row[1]

    plt.figure()
    
    ax1 = plt.subplot(211)
    ax1.plot(y_points, xvel_avg, 'bo')
    plt.title("Average x velocity as a function of y position")

    ax2 = plt.subplot(212)
    ax2.plot(y_points, count_list, 'ro')
    ax2.plot(y_points, np.ones(len(y_points)) * sum(count_list) / len(count_list), 'b-')
    plt.title("Number of particles as a function of y position")
    ax2.set_ylim([0, sum(count_list)/4])
    
plt.show()
plt.close()
