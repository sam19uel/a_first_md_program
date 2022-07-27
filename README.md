# a_first_md_program
This folder is a reference molecular dynamics folder about writing a first md program.

Outline of how this folder is organised.
The main reference for this organisation comes from "Understanding Molecular Simulation by Daan Frenkel and Bernard Smith ".
This folder is a reference molecular dynamics folder about writing a first md program.

note: The md_program as taught in the book is organised into a molecular dynamics look that performs the different steps as defined as subroutines. 

The book teaches them in pseudo algorithms such as to not make the commitment to specific language and writes it an algorithimic way that one can adapt into a code, no matter what the language choice is of the user. 
The code I am using is written in $python$ just because it's what I am using the most at the time of writing this. 
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





