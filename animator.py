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
from matplotlib.patches import Circle

ANIMATE         = True
TILE            = True
PLOT_TEMP       = True
PLOT_MAXBOLTZ   = False
PLOT_CROSSFLOW  = False

NUM_INPUTS      = 15
R_WATER         = 0.5

R_PARTICLE_SOFT = 2.5   # Particles / Water objects attract at this shell
R_PARTICLE      = 2.225 # Particles repel at this shell
R_PARTICLE_INNER= 1.95  # Water objects repel at this shell

#------------------------------------------------------------

def parse_whole_output(filename):
    with open(filename,"r") as f:
        data = f.readlines()

    # Chop off first line which is the column header
    data = data[1:]

    # Split each line by commas, and then convert each entry from str->float
    for index, line in enumerate(data):
        data[index] = line.split(",")
        for inner_index, entry in enumerate(data[index]):
            data[index][inner_index] = float(entry)

    # Convert to np array to use multiple indexing syntax and return
    return np.array(data)

def parse_output_final_timestep(filename, n_obj, n_step):

     # Add 1 to account for the column header row
    num_lines = n_obj * n_step + 1
    # Grab only the last N_WATER steps
    with open("water.out","r") as f:
        data_last = [line for index, line in enumerate(f) if index > (num_lines - N_WATER - 1)]

    # Split each line by commas, and then convert each entry from str->float
    for index, line in enumerate(data_last):
        data_last[index] = line.split(",")
        for inner_index, entry in enumerate(data_last[index]):
            data_last[index][inner_index] = float(entry)

    # Convert to np array to use multiple indexing syntax and return
    return np.array(data_last)

# Takes in a data list and outputs the x and y coordinates of the objects in the
# list. N represents the number of objects, i is the time at which we want these
# positions. The replication of the lsits with various displacements represent
# the tilings. 2 and 3 represent the fact that x and y positions are represented
# as in the 2 and 3 entries of each row. 
def data_tiling(i, data_list, N, range_x, xmin):

    if N == 0:
        return ([],[])

    t = (i + wat_data[0, 0]) * TIME_STEP

    x = (list(data_list[N*i:N*i+N, 2])                                              +  # mid
        [(a+t*V_SHEAR-xmin)%range_x+xmin          for a in data_list[N*i:N*i+N, 2]] +  # top                                       
        [(a-t*V_SHEAR-xmin)%range_x+xmin          for a in data_list[N*i:N*i+N, 2]] +  # bot
        [a-BOX_SIZE                               for a in data_list[N*i:N*i+N, 2]] +  # left
        [(a-BOX_SIZE+t*V_SHEAR-xmin)%range_x+xmin for a in data_list[N*i:N*i+N, 2]] +  # top left
        [(a-BOX_SIZE-t*V_SHEAR-xmin)%range_x+xmin for a in data_list[N*i:N*i+N, 2]] +  # bot left
        [a+BOX_SIZE                               for a in data_list[N*i:N*i+N, 2]] +  # right
        [(a+BOX_SIZE+t*V_SHEAR-xmin)%range_x+xmin for a in data_list[N*i:N*i+N, 2]] +  # top right
        [(a+BOX_SIZE-t*V_SHEAR-xmin)%range_x+xmin for a in data_list[N*i:N*i+N, 2]])   # bot right
                
    y = (list(data_list[N*i:N*i+N, 3])                +                                # mid
        [a+BOX_SIZE for a in data_list[N*i:N*i+N, 3]] +                                # top
        [a-BOX_SIZE for a in data_list[N*i:N*i+N, 3]] +                                # bot
        list(data_list[N*i:N*i+N, 3])                 +                                # left
        [a+BOX_SIZE for a in data_list[N*i:N*i+N, 3]] +                                # top left
        [a-BOX_SIZE for a in data_list[N*i:N*i+N, 3]] +                                # bot left
        list(data_list[N*i:N*i+N, 3])                 +                                # right
        [a+BOX_SIZE for a in data_list[N*i:N*i+N, 3]] +                                # top right
        [a-BOX_SIZE for a in data_list[N*i:N*i+N, 3]])                                 # bot right

    return (x,y)

def data_nontile(i, data_list, N):

    if N == 0:
        return ([],[])

    return (data_list[N*i : N*i+N, 2],
            data_list[N*i : N*i+N, 3])

#------------------------------------------------------------
### READ FROM INPUT FILE ###

with open("inputs.in", "r") as g:
    inputs = g.readlines()

# Only have NUM_INPUTS inputs, grab and convert to float
input_data = inputs[:NUM_INPUTS]
input_data = list(map(float, input_data))

# Skip 1 line because it's a space, the lines below are the labels
labels = inputs[NUM_INPUTS + 1:]

# Chop out the part before the first space
for index, string in enumerate(labels):
    string = string.strip()
    bef,space,aft = string.partition(" ")
    labels[index] = aft

# Associate each input float with its label
input_dict = dict(zip(labels,input_data))

# Grab the relevant ones and give them nicknames. Some have to be ints to be
# indexed.
N_WATER     = int(input_dict["N_WATER"])
N_PARTICLES = int(input_dict["N_PARTICLES"])
BOX_SIZE    = input_dict["BOX_SIZE"]
N_STEPS     = int(input_dict["N_STEPS"])
TIME_STEP   = input_dict["TIME_STEP"]
V_SHEAR     = input_dict["V_SHEAR"]

#------------------------------------------------------------
### READ FROM WATER FILE ###

# Read the entire water file to get the simulation results for all time
if ANIMATE:
    wat_data  = parse_whole_output("water.out")
    part_data = parse_whole_output("particles.out")

# Only read the last time step of water results
if PLOT_MAXBOLTZ or PLOT_CROSSFLOW:
   wat_data_last = parse_output_final_timestep("water.out", N_WATER, N_STEPS)
    
#------------------------------------------------------------
### READ FROM TEMP FILE ###

if PLOT_TEMP:
    # Read temp file to get the temp results
    with open("temp.out","r") as t:
        temp_data = t.readlines()

    # Split by space and grab the last part, which is the temp, and then convert
    # to float.
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

    ax_animate = fig.add_subplot(111, aspect='equal', autoscale_on=False, 
                                xlim=limits, ylim=limits)

    # water holds the locations of the water, centres_w is a small dot at its 
    # centre.
    if TILE:
        num_water_circles = 9 * N_WATER
        num_particles_circles = 9 * N_PARTICLES
    else:
        num_water_circles = N_WATER
        num_particles_circles = N_PARTICLES

    water = [Circle((0, 0), R_WATER, fc='blue') for _ in range(num_water_circles)]
    particles = [Circle((0, 0), R_PARTICLE, fc='brown') for _ in range(num_particles_circles)]
    particles_shell = [Circle((0, 0), R_PARTICLE_SOFT, fc='brown', alpha=0.5) for _ in range(num_particles_circles)]
    particles_inner = [Circle((0, 0), R_PARTICLE_INNER, fc='yellow', alpha=0.5) for _ in range(num_particles_circles)]

    for patch in water:
        ax_animate.add_patch(patch)
    for patch in particles:
        ax_animate.add_patch(patch)
    for patch in particles_shell:
        ax_animate.add_patch(patch)
    for patch in particles_inner:
        ax_animate.add_patch(patch)

    centres_w, = ax_animate.plot([], [], 'ko')
    centres_p, = ax_animate.plot([], [], 'yo')

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

        rect_bl2 = plt.Rectangle((-BOX_SIZE, -BOX_SIZE), BOX_SIZE, BOX_SIZE, ec='none', lw=2, fc='none')
        ax_animate.add_patch(rect_bl2)

        rect_br = plt.Rectangle((BOX_SIZE, -BOX_SIZE), BOX_SIZE, BOX_SIZE, ec='none', lw=2, fc='none')
        ax_animate.add_patch(rect_br)

    # Convert limits to array to do scalar multiplication, then unpack into the position args for text
    # Times 0.95 to be almost at the border but not quite
    text_pos = 0.95 * np.array(limits)
    time_text = ax_animate.text(*text_pos, '', fontsize=12)

    def init():
        """Initialize animation"""
        rect.set_edgecolor('k')

        if TILE:
            rect_top.set_edgecolor('k')
            rect_tl.set_edgecolor('k')
            rect_tl2.set_edgecolor('k')
            rect_tr.set_edgecolor('k')
            rect_bot.set_edgecolor('k')
            rect_bl.set_edgecolor('k')
            rect_bl2.set_edgecolor('k')
            rect_br.set_edgecolor('k')

        centres_w.set_markersize(2) 
        centres_p.set_markersize(2)

        if TILE:
            return tuple(water) + tuple(particles) +  tuple(particles_shell) + tuple(particles_inner) + (centres_w, centres_p, rect, rect_top, rect_tl, rect_tl2, rect_tr, rect_bot, rect_bl, rect_bl2, rect_br, time_text)
        else:
            return tuple(water) + tuple(particles) +  tuple(particles_shell) + tuple(particles_inner) + (centres_w, centres_p, rect, time_text)


    def animate(i):
        """perform animation step"""

        xmin, xmax = ax_animate.get_xlim()
        range_x = xmax - xmin

        t = (i + wat_data[0, 0]) * TIME_STEP

        if TILE: 
            
            # Set the shifting positions of the boxes
            rect_top.set_x((            t * V_SHEAR + BOX_SIZE) % (3*BOX_SIZE) - BOX_SIZE)
            rect_bot.set_x((          - t * V_SHEAR + BOX_SIZE) % (3*BOX_SIZE) - BOX_SIZE)
            rect_tl.set_x ((-BOX_SIZE + t * V_SHEAR + BOX_SIZE) % (3*BOX_SIZE) - BOX_SIZE)
            rect_tl2.set_x((            t * V_SHEAR)            % BOX_SIZE     - 2*BOX_SIZE)
            rect_bl.set_x ((-BOX_SIZE - t * V_SHEAR + BOX_SIZE) % (3*BOX_SIZE) - BOX_SIZE)
            rect_bl2.set_x((          - t * V_SHEAR           ) % BOX_SIZE     - 2*BOX_SIZE)
            rect_tr.set_x (( BOX_SIZE + t * V_SHEAR + BOX_SIZE) % (3*BOX_SIZE) - BOX_SIZE)
            rect_br.set_x (( BOX_SIZE - t * V_SHEAR + BOX_SIZE) % (3*BOX_SIZE) - BOX_SIZE)
            
            # Calculate data of the water & particles using tiling
            # Get data from the function, which returns in a tuple of 2 lists 
            # representing the x and y values. 
            water_xy_lists = data_tiling(i, wat_data, N_WATER, range_x, xmin)
            particles_xy_lists = data_tiling(i, part_data, N_PARTICLES, range_x, xmin)

        else: # Only plots the main cell   
            # Calculate data of the water & particles using nontile
            water_xy_lists = data_nontile(i, wat_data, N_WATER)
            particles_xy_lists = data_nontile(i, part_data, N_PARTICLES)

        # Use zip (with unpacking of the tuple) to make it into ordered (x,y) 
        # pairs instead, and then call tuple() to convert from generator to 
        # tuple.
        water_xy = tuple(zip(*water_xy_lists))
        particles_xy = tuple(zip(*particles_xy_lists))

        # Assign these (x,y) pairs to the Circle objects
        for index, patch in enumerate(water):
            patch.center = water_xy[index]
        for index, patch in enumerate(particles):
            patch.center = particles_xy[index]
        for index, patch in enumerate(particles_shell):
            patch.center = particles_xy[index]
        for index, patch in enumerate(particles_inner):
            patch.center = particles_xy[index]

        # These are simpler to assign, just assign the tuple of lists directly
        # to the points that are plotted.
        centres_w.set_data(water_xy_lists)
        centres_p.set_data(particles_xy_lists)

        time_text.set_text('Time = {0} / {1}'.format(i, N_STEPS))
        
        if TILE:
            return tuple(water) + tuple(particles) +  tuple(particles_shell) + tuple(particles_inner) + (centres_w, centres_p, rect, rect_top, rect_tl, rect_tl2, rect_tr, rect_bot, rect_bl, rect_bl2, rect_br, time_text)
        else:
            return tuple(water) + tuple(particles) +  tuple(particles_shell) + tuple(particles_inner) + (centres_w, centres_p, rect, time_text)


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
    v_list = np.sqrt(wat_data_last[:, 5]**2+ wat_data_last[:, 6]**2 + wat_data_last[:, 7]**2)

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
    
    # Create an array with the y-pos and the x-vel (3 and 5 positions)
    y_vx_list = wat_data_last[:, [3,5]]
    # Takes the first index of each row, sort them, and then index the original
    # array using that order. In other words, sorts the rows by their y-pos.
    y_vx_list = y_vx_list[y_vx_list[:,0].argsort()]

    # Tells us which bin each row belongs in
    bin_indices = np.searchsorted(y_bins, y_vx_list[:,0]) - 1

    xvel_entries = np.zeros([num_bins, 2])
    xvel_avg = np.zeros(num_bins)

    # Sum up the y-vel for each bin by referring to the bin_indices list
    # y_bins_entries now has a list of total y-vel in that bin + number of water blobs in that bin
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
    plt.title("Number of water as a function of y position")
    ax2.set_ylim([0, sum(count_list)/4])
    
plt.show()
plt.close()
