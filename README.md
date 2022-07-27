# a_first_md_program
This folder is a reference molecular dynamics folder about writing a first md program.

Outline of how this folder is organised.
The main reference for this organisation comes from "Understanding Molecular Simulation by Daan Frenkel and Bernard Smith ".
This folder is a reference molecular dynamics folder about writing a first md program.
This is NOT a tutorial, but meant to be a companion to the book. 

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

def init(self):

        sumv = np.zeros(self.ndim)
        sumv2 = 0
        
        # Particel positions and velocities:
        # Note, rand takes from a uniform dis of [0,1)
        ## Positions
        # Random
        # self.x = np.random.uniform(0,self.box,(self.npart,self.ndim))

        # 2D: X- Positions Placed on Grid, Y taken from Uniform Dist
        x_axis = np.linspace(0,self.box, self.npart)
        y_axis = np.random.uniform(0,self.box, self.npart)
        
        for n in range(self.npart):
            self.x[n][0] = x_axis[n]
            self.x[n][1] = y_axis[n]


        ## Velocity Options: 
        ## Uniform Dis
        #self.v = np.random.uniform(-0.5,0.5,(self.npart,self.ndim))
        # Boltzmann
        factor = math.sqrt(self.k * self.temp / self.mass)
        self.v = np.random.normal(loc=0,scale=factor,size=(self.npart, self.ndim))

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


