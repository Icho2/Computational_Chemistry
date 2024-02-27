def bond_lengths(particle_index, pos):
    r = [[0 for _ in range(len(pos))] for _ in range(len(pos))]
    for j in range(len(pos)):
        xij = pos[particle_index][0][-1] - pos[j][0][-1]
        yij = pos[particle_index][1][-1] - pos[j][1][-1]
        zij = pos[particle_index][2][-1] - pos[j][2][-1]
        r[particle_index][j] = np.sqrt(xij**2 + yij**2 + zij**2)
    return r    

def lj_3d(particle_index, pos, r, sigma=1.0, epsilon=1.0):
    #the force vector containing all 3D forces acting on a specified particle
    f = [0, 0, 0]
    for i in range(len(pos)):
        if i != particle_index and r[particle_index][i] != 0:
                f[0] += -(pos[particle_index][0][-1] - pos[i][0][-1])*((24*epsilon*(sigma**6))/(r[particle_index][i]**8))*(1 - 2*((sigma/r[particle_index][i])**6))
                f[1] += -(pos[particle_index][1][-1] - pos[i][1][-1])*((24*epsilon*(sigma**6))/(r[particle_index][i]**8))*(1 - 2*((sigma/r[particle_index][i])**6))
                f[2] += -(pos[particle_index][2][-1] - pos[i][2][-1])*((24*epsilon*(sigma**6))/(r[particle_index][i]**8))*(1 - 2*((sigma/r[particle_index][i])**6)) 
    return f

def multiple_atoms_Verlet(pos, vel, h, end_time, m):
    t = [0]
    while t[-1] < end_time:
        for i in range(len(pos)):
            #Compute all of our r_ij b/c we need them for computing F
            r = bond_lengths(i, pos)
            f = lj_3d(i, pos, r)
            
            # step 1: calculate
            pos[i][0].append(pos[i][0][-1] + (h*vel[i][0][-1]) + (((h**2)*f[0])/(2*m[i])))#x1(k+1)
            pos[i][1].append(pos[i][1][-1] + (h*vel[i][1][-1]) + (((h**2)*f[1])/(2*m[i])))#y1(k+1)
            pos[i][2].append(pos[i][2][-1] + (h*vel[i][2][-1]) + (((h**2)*f[2])/(2*m[i])))#z1(k+1)

        # step 2: evaluate        
        for i in range(len(pos)):
            # Because we have computed new positions for each particle, we must re-compute our new bond lengths.
            r = bond_lengths(i, pos)
            f_new = lj_3d(i, pos, r)
            # step 3: calculate
            vel[i][0].append(vel[i][0][-1] + ((h/(2*m[i])) * (f[0] + f_new[0])))
            vel[i][1].append(vel[i][1][-1] + ((h/(2*m[i])) * (f[1] + f_new[1])))
            vel[i][2].append(vel[i][2][-1] + ((h/(2*m[i])) * (f[2] + f_new[2])))
        t.append(t[-1] + h)
    return t, pos, vel

# def multiple_atoms_Verlet(pos, vel, h, end_time, m):
A = [[0.0], [0.0], [0.0]] 
B = [[2.0], [0.0], [0.0]] #format particle name = [[x], [y], [z]]
pos = [A, B]
A_v = [[-0.1], [0.0], [0.0]]
B_v = [[0.0], [0.0], [0.0]] #format {particle name}_v  = [[v_x], [v_y], [v_z]]
vel = [A_v, B_v]
h = 0.02 # homework step size diverges ?
m = [1, 1]
end_time = 20
t, pos, vel = multiple_atoms_Verlet(pos, vel, h, end_time, m)

A = np.array(pos[0])
B = np.array(pos[1])
vel = np.array(vel)
r_AB = np.sqrt((A - B)**2)
plt.plot(t, r_AB[0], label="Diatomic molecule");plt.xlabel("time (s)");plt.ylabel("$r_{AB}$");plt.legend();plt.title("Velocity Verlet")
plt.show()

# The second plot is here to confirm that our system is aligning with the potential.
energy_at_end = lj_potential(r_AB[0][-1])
energy_at_start = lj_potential(r_AB[0][0])
plt.scatter(r_AB[0][-1], energy_at_end)
plt.scatter(r_AB[0][0], energy_at_start)
plt.plot(x1, potential, label="Lennard-Jones Potential") 
plt.xlabel('r');plt.legend();plt.grid()
plt.show()

def energy_of_system(t, m, v, r_AB, epsilon = 1.0, sigma = 1.0):
    energy = list()
    for i in range(0, len(t)):
        energy_new = ((1/2)*m[0]*(v[0][0][i]**2)) + ((1/2)*m[1]*(v[1][0][i]**2)) + (4*epsilon*((sigma/r_AB[0][i])**12 - (sigma/r_AB[0][i])**6))
        #energy_new = ((1/2)*m[0]*(v[0][0][i]**2)) - ((1/2)*m[1]*(v[1][0][i]**2)) + (4*epsilon*((sigma/r_AB[0][i])**12 - (sigma/r_AB[0][i])**6))
        if i == 0:
            print("Initial energy of our system:",energy_new)
        energy.append(energy_new)    
    return energy

energy = energy_of_system(t, m, vel, r_AB)
print("minimum energy reached at time step",h,":", min(energy))
print(vel[0][0][-1])
print(vel[1][0][-1])
print(r_AB[0][-1])
plt.plot(t, energy);plt.xlabel("time (s)");plt.ylabel("Energy");plt.title("Conservation of energy");plt.show() 
