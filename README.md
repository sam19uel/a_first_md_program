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
        # 2D: X- Positions Placed on Grid, Y taken from Uniform Dist
        if(xswitch == 0):
        
            x_axis = np.linspace(0,self.box, self.npart)
            y_axis = np.random.uniform(0,self.box, self.npart)
            
            for n in range(self.npart):
                self.x[n][0] = x_axis[n]
                self.x[n][1] = y_axis[n]

        # xswitch = 1
        # Completely Random
        if (xswitch == 1):
        
            self.x = np.random.uniform(0,self.box,(self.npart,self.ndim))


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
In my code (by default switches = 0): positions x, x-componenets on an equally spaced lattice, y-components randomnly ; velocities v initialized from the Boltzmann distribution. 

In theory there are many different ways one can intialize and many different reasons to do so.
One however needs to be careful, depending on the method about the treatment. 

For example : 
### Minimization

Because I do not choose to put the particles on a lattice, one possibility is that the particles find themselves too close to one another, and this leads to an explosion of the potential energy, so I iterate through the particles and minimize the potential energy by shifting the postiions such that potential energy is minimized for a specific number of steps. To do this, there is also a subroutine that calculated the potential energy *pe*, which is a LJ-potential with a cutoff specified by one of the variables of the simulation.

```python
def minimize(self, min_steps, max_dr = 0.15):
        for i in range(min_steps):
            p = np.random.randint(0,self.npart)
            dr = np.random.random_sample(size=2)*2*max_dr - max_dr
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


