# a_first_md_program
This folder is a reference molecular dynamics folder about writing a first md program.

Outline of how this folder is organised.
The main reference for this organisation comes from "Understanding Molecular Simulation by Daan Frenkel and Bernard Smith ".
This folder is a reference molecular dynamics folder about writing a first md program.

**This is NOT a tutorial, but meant to be a companion to the book.**

note: The md_program as taught in the book is organised into a molecular dynamics look that performs the different steps as defined as subroutines. 

The book teaches them in pseudo algorithms such as to not make the commitment to specific language and writes it an algorithimic way that one can adapt into a code, no matter what the language choice is of the user. 
The code I am using is written in `python` just because it's what I am using the most at the time of writing this. 
Therefore the code is organised in a similar way such that one can use it as a reference code / companion code to the pseudo algorithms.


## md_loop

```python
## Running the Simulation MD Loop

# Creating the MD Object
sim1 = md()

# Initialisation
sim1.init()


while (sim1.t < sim1.tmax):
    # Determine the Forces
    sim1.force()

    # Intergrate Equations of Motion
    sim1.integrate()

    # Update time t
    print('t = ', sim1.t,' with tmax = ', sim1.tmax)
    sim1.t = sim1.t + sim1.dt
    
    # Sampling Data Output if needed

```

Notes: So this is the principle MD Loop that performs the simuation. The hierarchical organistion is as follows. I defined a class called "md" which is a class that defines a md_simulation object, where the sub_routines of the book are defined as internal functions of the class. Where here we created a "md object" called "sim1".  

## md_class initialization

IMPORTANT: here I am refering to the class initialization of creating the md object, and not the initialization subroutine of the md as refered to in the book which we call "init". 

```python
def __init__(self):
        
        if (self.ndim == 1): 
            print('Warning: ndim == 1!')

        # number of dimensions
        self.ndim = 2
        #number of particles
        self.npart = 50

        # Temperature
        self.temp = 1

        # Initialize the time at t=0
        self.t = 0

        self.dt = 0.001
        # sometimes the time step is also called 'delt' and sometimes 'dt'.
        self.delt = self.dt
        # tmax of the simuation
        self.tmax = 0.5

        # Size of Simulation Box
        self.box = 10

        # Mass of Particles
        self.mass = 1

        # LJ-cutoff
        self.rc = 4
        self.rc2 = self.rc**2

        # Boltzmann const
        self.k = 1

        # ecut = value of LJ potential at r=rc
        self.ecut = 4*((1/(self.rc2**6))-(1/(self.rc2**3)))

        # Positions x, Pervious Positions xm and Velocities v
        self.x = np.zeros((self.npart,self.ndim))
        self.xm = np.zeros((self.npart,self.ndim))
        self.v = np.zeros((self.npart,self.ndim))

        # Forces on Particles
        self.f = np.zeros((self.npart,self.ndim))

        # Energy
        self.en = 0
        self.etot = 0

        # note: sumv has the dimensions of a velocity
        self.sumv = np.zeros(self.ndim)
        self.sumv2 = 0
```        

We define the necessary variables that are needed to define our md_simulation. 
Small note: this part is not explicitely mentioned in the book of what or how to do this declaration of variables or what is needed, so this is an intepretation of what is explained. There are probably multiple ways to do or organise this that are equivalent or more optimized. 

##  Initialization of the MD Program

This is the subroutine that initializes the molecular dynamics program.

```python
##  Initialization of the MD Program
    def init(self, xswitch = 0, vswitch = 0):

        # xswitch = a switch that decides how we initialize positions
        # vswitch = a switch that decides how we initialize velocities.

        sumv = np.zeros(self.ndim)
        sumv2 = 0
        
        #####################################
        ## Particel positions and velocities:
        # Note, rand takes from a uniform dis of [0,1)

        ## Position Initialization Options: 
        
        # xswitch = 0 DEFAULT
        # Completely Random
        if(xswitch == 0):

            
            for k in range(self.ndim):
                k_axis = np.random.uniform(0,self.box, self.npart)
                for n in range(self.npart):
                    
                    self.x[n][k] = k_axis[n]
        
            
                                        
        # xswitch = 1
        # Build a ndim lattice, then randomly choose a position from the lattice for each particle
        if (xswitch == 1):

            for k in range(self.ndim):
                k_axis = np.linspace(0,self.box, self.npart)
                for n in range(self.npart):
                    self.x[n][k] = random.choice(k_axis)


        ## Velocity Initialization Options: 
        
        # vswitch = 0 DEFALUT
        # Boltzmann Distribution set at the Temperature

        if (vswitch == 0):
            factor = math.sqrt(self.k * self.temp / self.mass)
            self.v = np.random.normal(loc=0,scale=factor,size=(self.npart, self.ndim))

        # vxwitch = 1
        # Uniform Dis
        if (vswitch == 1):
            self.v = np.random.uniform(-0.5,0.5,(self.npart,self.ndim))

        #####################################

        # adding them up
        for i in range(self.npart):
            sumv = sumv + self.v[i]
            sumv2 = sumv2 + (np.linalg.norm(self.v[i]))**2

        # Dividing
        # Velocity center of mass and mean-squared velocity
        sumv = sumv/self.npart
        sumv2 = sumv2/self.npart

        # scale factor of velocity
        fs = math.sqrt(3*self.temp/sumv2)

        # Setting the desired kinetic energy
        # and set the velocity center of mass to zero

        # previous postions initial estimation 
        for i in range(self.npart):
            self.v[i] = (self.v[i] - sumv)*fs
            self.xm[i] = self.x[i] - self.v[i]*self.dt

        self.sumv = sumv
        self.sumv2 = sumv2
```

Comments on Initializing Positions and Velocities: 

xswitch and vswitch : so this is an internal decision to choose how to initialize the velocities and positions. Depending on the switch number you can decide which method to use. 

In the book : positions x is initialized on a lattice ; velocities v is initialized from a unifrom distribution.
In my code (by default switches = 0): positions x, first build a ndim lattice then randomly assign a position to each particle from the lattice ; velocities v initialized from the Boltzmann distribution. 

In theory there are many different ways one can intialize and many different reasons to do so.
One however needs to be careful, depending on the method about the treatment. 

For example : 
### Minimization

One possibile danger is that the particles find themselves too close to one another, and this leads to an explosion of the potential energy, so I iterate through the particles and minimize the potential energy by shifting the postiions such that potential energy is minimized for a specific number of steps. To do this, there is also a subroutine that calculated the potential energy *pe*, which is a LJ-potential with a cutoff specified by one of the variables of the simulation.

```python
def minimize(self, min_steps, max_dr = 0.15):
        for i in range(min_steps):
            p = np.random.randint(0,self.npart)
            dr = np.random.random_sample(size=self.ndim)*2*max_dr - max_dr
            pe_bef = self.pe()
            self.x[p] += dr
            pe_after = self.pe()
            if(pe_bef < pe_after):
                self.x[p] -= dr
        
        # Since the minimizing step changes the positions, we need to re-estimate the "previous positions"
        for i in range(self.npart):
            self.xm[i] = self.x[i] - self.v[i]*self.dt

def pe(self):

    pe = 0

    # loop over all pairs to calculate the force
    for i in range(self.npart - 1):
        for j in range(i+1,self.npart):
            xr = self.x[i] - self.x[j]
            # periodic BC
            xr = xr - self.box*np.rint(xr/self.box)
            r2 = (np.linalg.norm(xr))**2

            # test cutoff
            if (r2 < self.rc2):
                r2i = 1/r2
                r6i = r2i**3
                # LJ Potential
                ff = 48*r2i*r6i*(r6i-0.5)
                pe = pe + 4*r6i*(r6i - 1) - self.ecut
    return pe
```

Note: this process can be avoided completely if one decides to place them on a lattice directly, which one can easily implement if needed in the code, and setting a specific switch for it. 

## Drawing the Particles

It's always useful (and cool) to understand what is going on in a simulation, so this is a subroutine that draws the current positions of the particles.

```python
    def draw_particles(self):
    
        # 1D
        if (self.ndim == 1):

            print('Hard to draw 1D. :/')

        if (self.ndim == 2):
            # Determine Appropriate Size of Figure:
            plt.figure(figsize=(5,5))
            axis = plt.gca()
            
            axis.set_xlim(-10,self.box+10)
            axis.set_ylim(-10,self.box+10)


            for i in range(self.npart):
                axis.add_patch( plt.Circle(self.x[i], radius=0.5, linewidth=2, edgecolor='black') )
            plt.show()

        if (self.ndim == 3):

            fig = plt.figure()
            axis = fig.add_subplot(projection='3d')

            axis.set_xlim(-5,self.box+5)
            axis.set_ylim(-5,self.box+5)
            axis.set_zlim(-5,self.box+5)

            for i in range(self.npart):
                axis.scatter( self.x[i][0], self.x[i][1], self.x[i][2], marker="o", c='blue', linewidth=2, edgecolor='black')
            plt.show()
```

## Calculation of the Forces

A coded version of the pseudo algorithm that calculates and updates the forces on each particle.

```python3
## Force
    # Determine the force and energy
    # box = diameter of periodic box
    def force(self):

        en = 0
        
        # Forces set already to zero at initialization of the MD Object.
        # Reset Forces to Zero
        self.f = np.zeros((self.npart,self.ndim))
        
        # loop over all pairs to calculate the force
        for i in range(self.npart - 1):
            for j in range(i+1,self.npart):
                xr = self.x[i] - self.x[j]
                # periodic BC
                xr = xr - self.box*np.rint(xr/self.box)
                r2 = (np.linalg.norm(xr))**2
            
            
            
                # test cutoff
                if (r2 < self.rc2):
                    r2i = 1/r2
                    r6i = r2i**3
                    # LJ Potential
                    ff = 48*r2i*r6i*(r6i-0.5)
                    #Update force
                    self.f[i] = self.f[i] + ff*xr
                    self.f[j] = self.f[j] - ff*xr
                    en = en + 4*r6i*(r6i - 1) - self.ecut
        
        self.en = en
```

## Integrate Equations of Motion

We integrate the equations of motion numerically and one can choose which algorithm to us. Here we use the Verlet algorithm as is one in the book.
We also calculate the instantaneous temperature, as well as compute the energies.

```python
def integrate(self):
        sumv = np.zeros(self.ndim)
        sumv2 = 0

        # MD Loop
        for i in range(self.npart):
            # Verlet Algo
            xx = 2*self.x[i] - self.xm[i] + (self.delt**2)*self.f[i]
            # velocity
            self.v[i] = (xx-self.xm[i])/(2*self.delt)
            # velocity center of mass
            sumv = sumv + self.v[i]
            # total kin energy
            sumv2 = sumv2 + np.linalg.norm(self.v[i])**2

            # update positions previous time
            self.xm[i] = self.x[i]
            self.x[i] = xx
        
        # instantaneous temperature
        self.temp = sumv2/(3*self.npart)

        # total energy per particle
        self.etot = (self.en + 0.5*sumv2)/self.npart

        # kin en
        self.sumv = sumv
        self.sumv2 = sumv2
```

## Sampling

In theory, now you have everything you need to run the simulation, but almost always we want to write out useful information, such as the temeperature, energies, or whatever other information that we compute at each (or every few) timesteps. 

This is done in the md_loop and here is a simple exaple of writing the potential energy, kinetic energy, total energy per particle and temperature. 

```python
##################################
## Running the Simulation MD Loop

# Creating the MD Object
sim1 = md()

# Initialisation
sim1.init()
sim1.draw_particles()

# Minimize Potential Energy
sim1.minimize(2000)
sim1.draw_particles()

# Sampling Data Output 
time = [0]
temp_output = [sim1.temp]
en_output = [sim1.en]
etot_output = [sim1.etot]
sumv2_output = [sim1.sumv2]


while (sim1.t < sim1.tmax):
    # Determine the Forces
    sim1.force()

    # Intergrate Equations of Motion
    sim1.integrate()

    # Update time t
    print('t = ', sim1.t,' with tmax = ', sim1.tmax)
    sim1.t = sim1.t + sim1.dt

    ## Sampling Averages
    # For time file
    time = np.append(time, sim1.t)
    # For Temperature File
    temp_output = np.append(temp_output,sim1.temp)
    # For en
    en_output = np.append(en_output, sim1.en)
    # For etot i.e total energy per particle
    etot_output = np.append(etot_output, sim1.etot)
    # For Kin Energy
    sumv2_output = np.append(sumv2_output, sim1.sumv2)


# Writing Output Files
# Temperature Output
# [time] [temperature]
temp_data = np.column_stack([time, temp_output])
np.savetxt('temperature.dat', temp_data, delimiter='   ', newline='\r\n')

# etot output
# [time] [etot]
etot_data = np.column_stack([time, etot_output])
np.savetxt('etot.dat', etot_data, delimiter='   ', newline='\r\n')

# en output
# [time] [en]
en_data = np.column_stack([time, en_output])
np.savetxt('en.dat', en_data, delimiter='   ', newline='\r\n')

# sumv2 output
sumv2_data = np.column_stack([time, sumv2_output])
np.savetxt('sumv2.dat', sumv2_data, delimiter='   ', newline='\r\n')
#################################

sim1.draw_particles()

```

## Full Code

The full code of the topics dicussed so far is in the repository under the name `md_class.py`.
