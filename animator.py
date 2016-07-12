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

ANIMATE         = False # Plot the particles moving around. 
TILE            = False # Plot the replicated units
ANIMATE_WATER   = False # Plot the water objects (since there are so many)
PLOT_TEMP       = True # Plot the temperature-time curve 
PLOT_MAXBOLTZ   = False # Plot the velocity distribution of water at the last time step
PLOT_CROSSFLOW  = False # Plot the x-velocity distribution of water against y-position at the last time step
PLOT_VISC       = True # Plot the viscosity-time curve
PLOT_SQ_DISP_W  = False # Plot the squared-displacement of water
PLOT_SQ_DISP_P  = False  # Plot the squared-displacement of particles
PLOT_PRESSURE   = False  # Plot the pressure-time curve

NUM_INPUTS      = 15    # How many inputs there are to be read from inputs.in

R_WATER         = 0.5   # Radius of water -- set to be 0.5 by definition of length scale
R_PARTICLE_SOFT = 2.5   # Particles / Water objects attract at this shell
R_PARTICLE      = 2.225 # Particles repel at this shell
R_PARTICLE_INNER= 1.95  # Water objects repel at this shell

#------------------------------------------------------------

# Reads the output file, strips out the first row, and returns all the remaining
# data in a numpy array with first index representing rows and 2nd index
# representing columns
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

# Reads the output file, and ONLY returns the n_obj rows of data in a numpy 
# array with first index representing rows and 2nd index representing columns.
# This corresponds with only reading the data from the last time step.
def parse_output_final_timestep(filename, n_obj):

    # Count how many lines in this file
    num_lines = 0
    with open(filename, "r") as f:
        while f.readline() != "":
            num_lines += 1
            
    # Grab only the last n_obj steps
    with open(filename,"r") as f:
        data_last = [line for index, line in enumerate(f) if index > (num_lines - n_obj - 1)]

    # Split each line by commas, and then convert each entry from str->float
    for index, line in enumerate(data_last):
        data_last[index] = line.split(",")
        for inner_index, entry in enumerate(data_last[index]):
            data_last[index][inner_index] = float(entry)

    # Convert to np array to use multiple indexing syntax and return
    return np.array(data_last)

# Takes in a data list and outputs the x and y coordinates of the objects in the
# list. N represents the number of objects, i is the time at which we want these
# positions. The replication of the lists with various displacements represent
# the tilings. 2 and 3 represent the fact that x and y positions are represented
# as in the 2 and 3 entries of each row. 
def data_tiling(i, data_list, N, range_x, xmin):

    if N == 0:
        return ([],[])

    try:
        t = (i + part_data[0, 0]) * TIME_STEP
    except (IndexError, NameError):
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

# Takes in a data_list which holds the details of each particle at each time
# and outputs the x and y coordinates of all objects at time i. 2 and 3 
# represent the fact that x and y positions are represented as in the 2nd and 
# 3rd entries of each row. 
def data_nontile(i, data_list, N):

    if N == 0:
        return ([],[])

    return (data_list[N*i : N*i+N, 2],
            data_list[N*i : N*i+N, 3])

# Plots the squared displacemnts of particles from a particular data list
# ax specifies which axes to plot on, colour is a string for the colour of the
# lines that we should draw
def plot_sq_disp(data_list, N, ax, colour):
    # Repeat this loop for every particle present
    for i in range(N):

        # Take initial data from row i --> Particle i
        # 2, 3, 4 refer to the x, y, z coordiates being in these respective
        # columns of the data
        x0 = data_list[i, 2]
        y0 = data_list[i, 3]
        z0 = data_list[i, 4]

        x_prev = x0
        y_prev = y0
        z_prev = z0

        pass_x = 0
        pass_y = 0
        pass_z = 0

        disp_sq = []
        
        # Skip with a stride of N to take data only for the i-th object
        for row in data_list[i::N]:

            # Confirm it's object i
            assert(row[1] == i)

            x = row[2] + pass_x * BOX_SIZE
            y = row[3] + pass_y * BOX_SIZE
            z = row[4] + pass_z * BOX_SIZE

            # If x moved off the right side, then x_prev would be large while 
            # x would be ~0, therefore (x - x_prev) < 0, and we need to add 1 to
            # our new x value to reflect that it moved off the right.
            if (abs(x - x_prev)) > 0.8 * BOX_SIZE:
                pass_x -= np.sign(x - x_prev)
                x = row[2] + pass_x * BOX_SIZE
            if (abs(y - y_prev)) > 0.8 * BOX_SIZE:
                pass_y -= np.sign(y - y_prev)
                y = row[3] + pass_y * BOX_SIZE
            if (abs(z - z_prev)) > 0.8 * BOX_SIZE:
                pass_z -= np.sign(z - z_prev)
                z = row[4] + pass_z * BOX_SIZE

            disp_sq.append((x - x0) ** 2 + (y - y0) ** 2 + (z - z0) ** 2)

            x_prev = x
            y_prev = y
            z_prev = z

        ax.plot(disp_sq, c=colour)

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
### READ FROM DYN FILE ###

# Read the entire water and partcle file to get the simulation results for all time
if ANIMATE or PLOT_SQ_DISP_P:
    part_data = parse_whole_output("particles.out")
if ANIMATE_WATER or PLOT_SQ_DISP_W:
    wat_data  = parse_whole_output("water.out")

# Only read the last time step of water results
if PLOT_MAXBOLTZ or PLOT_CROSSFLOW:
   wat_data_last = parse_output_final_timestep("water.out", N_WATER)
    
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
### READ FROM VISC FILE ###

if PLOT_VISC:
    # Read visc file to get the visc results
    with open("shear.out","r") as s:
        shear_data = s.readlines()

    # Split by comma and grab all columns except the first, which is the time
    shears = np.array([shear_entry.split(',')[1:] for shear_entry in shear_data])
    shears = shears.astype(float)

    shear_xy = shears[:,0]
    shear_xz = shears[:,1]
    shear_yz = shears[:,2]

    times = np.arange(0, N_STEPS * TIME_STEP, TIME_STEP)

#------------------------------------------------------------
### READ FROM PRESSURE FILE ###

if PLOT_PRESSURE:
     # Read pressure file to get the pressure results
    with open("pressure.out","r") as p:
        pressure_data = p.readlines()

    # Split by space and grab the last part, which is the pressure, and then 
    # convert to float.
    pressure = [float(pressure_entry.split()[-1]) for pressure_entry in pressure_data]
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

    water = []
    centres_w = ()

    if ANIMATE_WATER:
        water = [Circle((0, 0), R_WATER, fc='blue') for _ in range(num_water_circles)]
        centres_w, = ax_animate.plot([], [], 'ko')
        for patch in water:
            ax_animate.add_patch(patch)

    particles = [Circle((0, 0), R_PARTICLE, fc='brown') for _ in range(num_particles_circles)]
    particles_shell = [Circle((0, 0), R_PARTICLE_SOFT, fc='brown', alpha=0.5) for _ in range(num_particles_circles)]
    particles_inner = [Circle((0, 0), R_PARTICLE_INNER, fc='yellow', alpha=0.5) for _ in range(num_particles_circles)]

    for patch in particles:
        ax_animate.add_patch(patch)
    for patch in particles_shell:
        ax_animate.add_patch(patch)
    for patch in particles_inner:
        ax_animate.add_patch(patch)

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

        if ANIMATE_WATER:
            centres_w.set_markersize(2) 
        centres_p.set_markersize(2)

        animated_items = tuple(particles) +  tuple(particles_shell) + tuple(particles_inner) + (centres_p, rect, time_text)

        if ANIMATE_WATER:
            # Don't use += so that water is in front and so plotted below
            animated_items = tuple(water) + animated_items
            animated_items += (centres_w,)
        if TILE:
            animated_items += (rect_top, rect_tl, rect_tl2, rect_tr, rect_bot, rect_bl, rect_bl2, rect_br)

        return animated_items


    def animate(i):
        """perform animation step"""

        xmin, xmax = ax_animate.get_xlim()
        range_x = xmax - xmin

        try:
            t = (i + part_data[0, 0]) * TIME_STEP
        except (IndexError, NameError):
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
            if ANIMATE_WATER:
                water_xy_lists = data_tiling(i, wat_data, N_WATER, range_x, xmin)
            particles_xy_lists = data_tiling(i, part_data, N_PARTICLES, range_x, xmin)

        else: # Only plots the main cell   
            # Calculate data of the water & particles using nontile
            if ANIMATE_WATER:
                water_xy_lists = data_nontile(i, wat_data, N_WATER)
            particles_xy_lists = data_nontile(i, part_data, N_PARTICLES)

        # Use zip (with unpacking of the tuple) to make it into ordered (x,y) 
        # pairs instead, and then call tuple() to convert from generator to 
        # tuple.
        if ANIMATE_WATER:
            water_xy = tuple(zip(*water_xy_lists))
        particles_xy = tuple(zip(*particles_xy_lists))

        # Assign these (x,y) pairs to the Circle objects
        if ANIMATE_WATER:
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
        if ANIMATE_WATER:
            centres_w.set_data(water_xy_lists)
        centres_p.set_data(particles_xy_lists)

        time_text.set_text('Time = {0} / {1}'.format(i, N_STEPS))
        
        animated_items = tuple(particles) +  tuple(particles_shell) + tuple(particles_inner) + (centres_p, rect, time_text)
        if ANIMATE_WATER:
            # Don't use += so that water is in front and so plotted below
            animated_items = tuple(water) + animated_items
            animated_items += (centres_w,)
        if TILE:
            animated_items += (rect_top, rect_tl, rect_tl2, rect_tr, rect_bot, rect_bl, rect_bl2, rect_br)

        return animated_items


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
    avg_quart_temp = np.mean(quart_temp)
    quart_temp_stdev = np.std(quart_temp)

    ax.text(0, 0, 'Average last quarter temp: {0:.5} $\\pm$ {1:.4}'.format(avg_quart_temp, quart_temp_stdev), fontsize=12)

if PLOT_VISC:
    plt.figure()
    ax = plt.subplot()

    # Chop each one to half its length to discard the first half
    shear_xy = shear_xy[len(shear_xy)//2:]
    shear_xz = shear_xz[len(shear_xz)//2:]
    shear_yz = shear_yz[len(shear_yz)//2:]
    times_chopped = times[len(times)//2:]

    # Number of time steps we are doing the integral over
    n_time_int = len(shear_xy)
    integrated_viscosity = 0

    # Delay time between left and right window (index)
    for t_d in range(n_time_int):
        # Start index for each window
        for t_start in range(n_time_int - t_d):
            integrated_viscosity += (shear_xy[t_start] * shear_xy[t_start + t_d]) / (n_time_int - t_d) 
            integrated_viscosity += (shear_xz[t_start] * shear_xz[t_start + t_d]) / (n_time_int - t_d) 
            integrated_viscosity += (shear_yz[t_start] * shear_yz[t_start + t_d]) / (n_time_int - t_d) 

    integrated_viscosity *= TIME_STEP * (BOX_SIZE ** 3)
    integrated_viscosity /= 3

    ax.plot(times_chopped, shear_xy)
    ax.plot(times_chopped, shear_xz)
    ax.plot(times_chopped, shear_yz)

    x_lim = ax.get_xlim()
    y_lim = ax.get_ylim()

    ax.text(x_lim[0], y_lim[0], 'Integrated viscosity over last {0} steps: {1:.4}'.format(n_time_int, integrated_viscosity), fontsize=12)

if PLOT_PRESSURE:
    plt.figure()
    ax = plt.subplot()

    quart_press = pressure[3*len(pressure)//4:]
    avg_quart_press = np.mean(quart_press)
    quart_press_stdev = np.std(quart_press)

    ax.plot(times, pressure)
    ax.plot(times, avg_quart_press * np.ones(len(times)), 'k', linewidth='3.0')
    ax.plot(times, (avg_quart_press + quart_press_stdev) * np.ones(len(times)), 'r--')
    ax.plot(times, (avg_quart_press - quart_press_stdev) * np.ones(len(times)), 'r--')

    x_lim = ax.get_xlim()
    y_lim = ax.get_ylim()

    ax.text(x_lim[0], y_lim[0], 'Average last quarter pressure: {0:.4} $\\pm$ {1:.4}'.format(avg_quart_press, quart_press_stdev), fontsize=12)

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
    num_boundaries = 11
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

    vx_theory = V_SHEAR * (y_bins / BOX_SIZE - 0.5)

    plt.figure()
    
    ax1 = plt.subplot(211)
    ax1.plot(y_points, xvel_avg, 'bo')
    ax1.plot(y_bins, vx_theory)
    plt.title("Average x velocity as a function of y position")

    ax2 = plt.subplot(212)
    ax2.plot(y_points, count_list, 'ro')
    ax2.plot(y_points, np.ones(len(y_points)) * sum(count_list) / len(count_list), 'b-')
    plt.title("Number of water as a function of y position")
    ax2.set_ylim([0, sum(count_list)/4])


if PLOT_SQ_DISP_W or PLOT_SQ_DISP_P:
    plt.figure()
    ax = plt.subplot()

    if PLOT_SQ_DISP_W:
        plot_sq_disp(wat_data, N_WATER, ax, 'blue')
    if PLOT_SQ_DISP_P:
        plot_sq_disp(part_data, N_PARTICLES, ax, 'brown')
    
plt.show()
plt.close()
