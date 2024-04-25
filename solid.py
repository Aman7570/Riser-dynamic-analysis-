from __future__ import division

import numpy as np
import csv
import os 
import precice

configuration_file_name = "../precice-config.xml"
participant_name = "Solid"
mesh_name_1 = "Solid-Mesh-1"
mesh_name_2 = "Solid-Mesh-2"
mesh_name_3 = "Solid-Mesh-3"
mesh_name_4 = "Solid-Mesh-4"
mesh_name_5 = "Solid-Mesh-5"
mesh_name_6 = "Solid-Mesh-6"
mesh_name_7 = "Solid-Mesh-7"
mesh_name_8 = "Solid-Mesh-8"
mesh_name_9 = "Solid-Mesh-9"
mesh_name_10 = "Solid-Mesh-10"


write_data_name_1 = 'Displacement-1'
read_data_name_1 = 'Force-1'
write_data_name_2 = 'Displacement-2'
read_data_name_2 = 'Force-2'
write_data_name_3 = 'Displacement-3'
read_data_name_3 = 'Force-3'
write_data_name_4 = 'Displacement-4'
read_data_name_4 = 'Force-4'
write_data_name_5 = 'Displacement-5'
read_data_name_5 = 'Force-5'
write_data_name_6 = 'Displacement-6'
read_data_name_6 = 'Force-6'
write_data_name_7 = 'Displacement-7'
read_data_name_7 = 'Force-7'
write_data_name_8 = 'Displacement-8'
read_data_name_8 = 'Force-8'
write_data_name_9 = 'Displacement-9'
read_data_name_9 = 'Force-9'
write_data_name_10 = 'Displacement-10'
read_data_name_10 = 'Force-10'


num_vertices =  1 # Number of vertices

solver_process_index = 0
solver_process_size = 1
participant= precice.Participant(participant_name, configuration_file_name,
                              solver_process_index, solver_process_size)

dimensions = participant.get_data_dimensions(mesh_name_1,read_data_name_1)
#dimensions_2 = participant.get_data_dimensions(mesh_name_2,read_data_name_1)
vertices_1= np.zeros((num_vertices, dimensions))
vertices_2= np.zeros((num_vertices, dimensions))
vertices_3= np.zeros((num_vertices, dimensions))
vertices_4= np.zeros((num_vertices, dimensions))
vertices_5= np.zeros((num_vertices, dimensions))
vertices_6= np.zeros((num_vertices, dimensions))
vertices_7= np.zeros((num_vertices, dimensions))
vertices_8= np.zeros((num_vertices, dimensions))
vertices_9= np.zeros((num_vertices, dimensions))
vertices_10= np.zeros((num_vertices, dimensions))



read_data_1 = np.zeros((num_vertices, dimensions))
write_data_1 = np.zeros((num_vertices, dimensions))
read_data_2 = np.zeros((num_vertices, dimensions))
write_data_2 = np.zeros((num_vertices, dimensions))
read_data_3 = np.zeros((num_vertices, dimensions))
write_data_3 = np.zeros((num_vertices, dimensions))
read_data_4 = np.zeros((num_vertices, dimensions))
write_data_4 = np.zeros((num_vertices, dimensions))
read_data_5 = np.zeros((num_vertices, dimensions))
write_data_5 = np.zeros((num_vertices, dimensions))
read_data_6 = np.zeros((num_vertices, dimensions))
write_data_6 = np.zeros((num_vertices, dimensions))
read_data_7 = np.zeros((num_vertices, dimensions))
write_data_7 = np.zeros((num_vertices, dimensions))
read_data_8 = np.zeros((num_vertices, dimensions))
write_data_8 = np.zeros((num_vertices, dimensions))
read_data_9 = np.zeros((num_vertices, dimensions))
write_data_9 = np.zeros((num_vertices, dimensions))
read_data_10 = np.zeros((num_vertices, dimensions))
write_data_10 = np.zeros((num_vertices, dimensions))




vertices_1= np.array([[0,0]])
vertices_2= np.array([[0,0]])
vertices_3= np.array([[0,0]])
vertices_4= np.array([[0,0]])
vertices_5= np.array([[0,0]])
vertices_6= np.array([[0,0]])
vertices_7= np.array([[0,0]])
vertices_8= np.array([[0,0]])
vertices_9= np.array([[0,0]])
vertices_10= np.array([[0,0]])


read_data_1 = vertices_1.copy()
write_data_1 = vertices_1.copy()
read_data_2 = vertices_2.copy()
write_data_2 = vertices_2.copy()
read_data_3 = vertices_3.copy()
write_data_3 = vertices_3.copy()
read_data_4 = vertices_4.copy()
write_data_4 = vertices_4.copy()
read_data_5 = vertices_5.copy()
write_data_5 = vertices_5.copy()
read_data_6 = vertices_1.copy()
write_data_6 = vertices_1.copy()
read_data_7 = vertices_2.copy()
write_data_7 = vertices_2.copy()
read_data_8 = vertices_3.copy()
write_data_8 = vertices_3.copy()
read_data_9 = vertices_4.copy()
write_data_9 = vertices_4.copy()
read_data_10 = vertices_5.copy()
write_data_10 = vertices_5.copy()








element_no = 60
n_s =10
w_o_s=0.0002
node_no = element_no + 1
dof = 2 * node_no
L = 13.12
l = L / element_no
Do = 0.028
Di = 0.027
D_r = 3000
E = 3.6 * (10 ** 9)
D_s = 1000
D_rf = 1000
Kbt = 0
Kbb = 0
#ftop = 1.3
#fb = 0.9
#f_b=1-fb
dt=0.0075
time=42
n = int(time / dt)
Ae = np.pi * (Do ** 2) / 4
Ai = np.pi * (Di ** 2) / 4
An = np.pi * ((Do ** 2) - (Di ** 2)) / 4
I = np.pi * ((Do ** 4) - (Di ** 4)) / 64
mf = D_rf * Ai
mr = D_r * An
m = (mr + mf)
g=9.81
#EI = E * I
EI=29.88
T=1600
#Elementary matrix
k = (EI / (l ** 3))*np.array([[12, 6 * l, - 12, 6 * l], [6 * l, 4 * l * l, -6 * l, 2 * l * l],
                  [-12, - 6 * l, 12, - 6 * l], [6 * l, 2 * l * l, - 6 * l, 4 * l * l]])

k1=EI*np.array([[12/l**3, 6/l**2, -12/l**3 ,6/l**2],
    [6/l**2 ,4/l ,-6/l**2 ,2/l],
    [-12/l**3 ,-6/l**2, 12/l**3, -6/l**2],
    [6/l**2 ,2/l, -6/l**2,4/l]])

kg = np.array([
    [6 / (5 * l), 1 / 10 ,- 6 / (5 * l), 1 / 10],
     [1 / 10,  2 * l / 15 ,- 1 / 10 ,- l / 30],
     [-6 / (5 * l), - 1 / 10,  6 / (5 * l), - 1 / 10],
     [1 / 10 ,- l / 30 ,- 1 / 10 , 2 * l / 15]])

M_l= m *np.array([[13*l/35, 11*(l**2)/210, 9*l/70, -13*(l**2)/420],
           [11*(l**2)/210 ,(l**3)/105, 13*(l**2)/420, -(l**3)/140 ],
           [9*l/70, 13*(l**2)/420, 13*l/35, -11*(l**2)/210 ],
          [-13*(l**2)/420 ,-(l**3)/140, -11*(l**2)/210 ,(l**3)/105]])

C_l=10*np.array([[13*l/35, 11*(l**2)/210, 9*l/70, -13*(l**2)/420],
           [11*(l**2)/210 ,(l**3)/105, 13*(l**2)/420, -(l**3)/140 ],
           [9*l/70, 13*(l**2)/420, 13*l/35, -11*(l**2)/210 ],
          [-13*(l**2)/420 ,-(l**3)/140, -11*(l**2)/210 ,(l**3)/105]])


f = np.array([6, l, 6, -l]).reshape(-1, 1)
f = ( l *1/ 12) * f
# Initialize arrays
q_f = np.zeros((element_no, 1))
M = np.zeros((dof, dof))
K = np.zeros((dof, dof))
C = np.zeros((dof, dof))

## matrix formulation
f_b=0.1
for i in range(1, element_no + 1):
    s = slice(2 * i - 2, 2 * i + 2)
    # t = (ftop * m * g * L) - (f_b * m * 9.81 * (L - (i * l)))
    t = T - (f_b * m * 9.81 * (L - (i * l)))
    K[s, s] = K[s, s] + k1+(t*kg)
    M[s, s] = M[s, s] + M_l
    C[s, s] = C[s, s] + C_l


##matrix reduction
b_c = [0,dof-2]
Kr = np.delete(K, b_c, axis=0)
Kr = np.delete(Kr, b_c, axis=1)
Mr = np.delete(M, b_c, axis=0)
Mr = np.delete(Mr, b_c, axis=1)
Cr = np.delete(C, b_c, axis=0)
Cr = np.delete(Cr, b_c, axis=1)

beta = 0.25
gamma = 0.5
a0 = 1 / (beta * dt * dt)
a1 = gamma / (beta * dt)
a2 = 1 / (beta * dt)
a3 = (1 / (2 * beta)) - 1
a4 = (gamma / beta) - 1
a5 = ((gamma / beta) - 2) * 0.5 * dt
a6 = (1 - gamma) * dt
a7 = gamma * dt
n = int(time / dt)
x = np.zeros((dof-2, n))
xd = np.zeros((dof-2, n))
xdd = np.zeros((dof-2, n))
x0 = np.zeros((dof - 2, 1))
v0 = np.zeros((dof - 2, 1))
k_cap = Kr + a0 * Mr + a1 * Cr
y = np.zeros((dof-2, n))
yd = np.zeros((dof-2, n))
ydd = np.zeros((dof-2, n))
y0 = np.zeros((dof - 2, 1))
vy0 = np.zeros((dof - 2, 1))
y[:, 0:1] = x0
yd[:, 0:1] = v0
Da = np.zeros((dof - 2, 1))
Va = np.zeros((dof - 2, 1))
Aa = np.zeros((dof - 2, 1))

def newmark_x(i, Fx):
    Da = x[:, i-1:i]
    Va = xd[:, i-1:i]
    Aa = xdd[:,i-1:i]
    p_cap = Fx + np.dot(Mr, (a0 * Da + a2 * Va + a3 * Aa)) + np.dot(Cr, (a1 * Da + a4 * Va + a5 * Aa))
    x[:, i:i+1] =np.dot(np.linalg.inv(k_cap), p_cap)
    xdd[:, i:i+1] = a0 * (x[:, i:i+1] - Da) - a2 * Va - a3 * Aa
    xd[:, i:i+1] = Va + a6 * Aa + a7 * xdd[:, i:i+1]

    # write_x = x[:, i]
def newmark_y(i, Fy):
    Da = y[:, i-1:i]
    Va = yd[:, i-1:i]
    Aa = ydd[:,i-1:i]
    p_cap = Fy + np.dot(Mr, (a0 * Da + a2 * Va + a3 * Aa)) + np.dot(Cr, (a1 * Da + a4 * Va + a5 * Aa))
    y[:, i:i+1] =np.dot(np.linalg.inv(k_cap), p_cap)
    ydd[:, i:i+1] = a0 * (y[:, i:i+1] - Da) - a2 * Va - a3 * Aa
    yd[:, i:i+1] = Va + a6 * Aa + a7 * ydd[:, i:i+1]
vertex_ids_1=participant.set_mesh_vertices(mesh_name_1, vertices_1)
vertex_ids_2=participant.set_mesh_vertices(mesh_name_2, vertices_2)
vertex_ids_3=participant.set_mesh_vertices(mesh_name_3, vertices_3)
vertex_ids_4=participant.set_mesh_vertices(mesh_name_4, vertices_4)
vertex_ids_5=participant.set_mesh_vertices(mesh_name_5, vertices_5)
vertex_ids_6=participant.set_mesh_vertices(mesh_name_6, vertices_6)
vertex_ids_7=participant.set_mesh_vertices(mesh_name_7, vertices_7)
vertex_ids_8=participant.set_mesh_vertices(mesh_name_8, vertices_8)
vertex_ids_9=participant.set_mesh_vertices(mesh_name_9, vertices_9)
vertex_ids_10=participant.set_mesh_vertices(mesh_name_10, vertices_10)
# vertex_ids_11=participant.set_mesh_vertices(mesh_name_11, vertices_11)
# vertex_ids_12=participant.set_mesh_vertices(mesh_name_12, vertices_12)
# vertex_ids_13=participant.set_mesh_vertices(mesh_name_13, vertices_13)
# vertex_ids_14=participant.set_mesh_vertices(mesh_name_14, vertices_14)
# vertex_ids_15=participant.set_mesh_vertices(mesh_name_15, vertices_15)
# vertex_ids_16=participant.set_mesh_vertices(mesh_name_16, vertices_16)
# vertex_ids_17=participant.set_mesh_vertices(mesh_name_17, vertices_17)
# vertex_ids_18=participant.set_mesh_vertices(mesh_name_18, vertices_18)
# vertex_ids_19=participant.set_mesh_vertices(mesh_name_19, vertices_19)
# vertex_ids_20=participant.set_mesh_vertices(mesh_name_20, vertices_20)
## disp initialize

TIME=0
itr=0
i=1
participant.initialize()

while participant.is_coupling_ongoing():
    if participant.requires_writing_checkpoint():
        print("DUMMY: Writing iteration checkpoint")

    
    dt = participant.get_max_time_step_size()
    
    TIME=TIME+dt
    ######
    #print("RECIEVED FORCE FROM OPENFOAM")
    read_data_1 = participant.read_data(mesh_name_1, read_data_name_1, [vertex_ids_1], dt)
    #print(read_data_1)
    read_data_2 = participant.read_data(mesh_name_2, read_data_name_2, [vertex_ids_2], dt)
    #print(read_data_2)
    read_data_3 = participant.read_data(mesh_name_3, read_data_name_3, [vertex_ids_3], dt)
    #print(read_data_3)
    read_data_4 = participant.read_data(mesh_name_4, read_data_name_4, [vertex_ids_4], dt)
    #print(read_data_4)
    read_data_5 = participant.read_data(mesh_name_5, read_data_name_5, [vertex_ids_5], dt)
    #print(read_data_5)
    read_data_6 = participant.read_data(mesh_name_6, read_data_name_6, [vertex_ids_6], dt)
    # print(read_data_6)
    read_data_7 = participant.read_data(mesh_name_7, read_data_name_7, [vertex_ids_7], dt)
    #print(read_data_7)
    read_data_8 = participant.read_data(mesh_name_8, read_data_name_8, [vertex_ids_8], dt)
    # print(read_data_8)
    read_data_9 = participant.read_data(mesh_name_9, read_data_name_9, [vertex_ids_9], dt)
    # print(read_data_9)
    read_data_10 = participant.read_data(mesh_name_10, read_data_name_10, [vertex_ids_10], dt)
    # print(read_data_10)
    # read_data_11 = participant.read_data(mesh_name_11, read_data_name_11, [vertex_ids_11], dt)
    # print(read_data_11)
    # read_data_12 = participant.read_data(mesh_name_12, read_data_name_12, [vertex_ids_12], dt)
    # print(read_data_12)
    # read_data_13 = participant.read_data(mesh_name_13, read_data_name_13, [vertex_ids_13], dt)
    # print(read_data_13)
    # read_data_14 = participant.read_data(mesh_name_14, read_data_name_14, [vertex_ids_14], dt)
    # print(read_data_14)
    # read_data_15 = participant.read_data(mesh_name_15, read_data_name_15, [vertex_ids_15], dt)
    # print(read_data_15)
    # read_data_16 = participant.read_data(mesh_name_16, read_data_name_16, [vertex_ids_16], dt)
    # print(read_data_16)
    # read_data_17 = participant.read_data(mesh_name_17, read_data_name_17, [vertex_ids_17], dt)
    # print(read_data_17)
    # read_data_18 = participant.read_data(mesh_name_18, read_data_name_18, [vertex_ids_18], dt)
    # print(read_data_18)
    # read_data_19 = participant.read_data(mesh_name_19, read_data_name_19, [vertex_ids_19], dt)
    # print(read_data_19)
    # read_data_20 = participant.read_data(mesh_name_20, read_data_name_20, [vertex_ids_20], dt)
    # print(read_data_20)




    if TIME >11:
        variables = []
        Fx = np.zeros((dof, 1))
        Fy = np.zeros((dof, 1))
        for j in range(1, n_s + 1):
            variables.append(locals()[f"read_data_{j}"])
        f1 = np.concatenate(variables)
        f1=f1/(w_o_s)
        print("Force")
        print(f1)
        qx = f1[:, 0:1]
        qy = f1[:, 1:2]

        fx = np.repeat(qx, (element_no / n_s), axis=0)
        fy = np.repeat(qy, (element_no / n_s), axis=0)
        h = 0
        for m in range(1, element_no + 1):
            s = slice(2 * m - 2, 2 * m + 2)
            Fx[s] = Fx[s] + (fx[h] * f)
            Fy[s] = Fy[s] + (fy[h] * f)
            h = h + 1
        Fx = np.delete(Fx, b_c, axis=0)
        Fy = np.delete(Fy, b_c, axis=0)
        if itr == 0:
            ax = (Fx - np.dot(Cr, v0) - np.dot(Kr, x0))
            ay = (Fy - np.dot(Cr, v0) - np.dot(Kr, x0))
            xdd[:, itr:itr + 1] = np.dot(np.linalg.inv(Mr), ax)
            ydd[:, itr:itr + 1] = np.dot(np.linalg.inv(Mr), ay)
        itr = itr + 1
        newmark_x(itr, Fx)
        newmark_y(itr, Fy)
        dis_x = x[:, i:i + 1].copy()
        dis_y = y[:, i:i + 1].copy() 
        strip_x_1 = np.zeros((dof, 1))
        strip_x = np.zeros((element_no, 1))
        strip_y_1 = np.zeros((dof, 1))
        strip_y = np.zeros((element_no, 1))
        #strip_x_1[2:dof - 2] = dis_x
        #strip_y_1[2:dof - 2] = dis_y
        strip_x_1[1:dof - 2,:] = dis_x[0:dof - 3,:]
        strip_x_1[dof - 1,:] = dis_x[dof - 3,:]
        strip_y_1[1:dof - 2,:] = dis_y[0:dof - 3,:]
        strip_y_1[dof - 1,:] = dis_y[dof - 3,:]
        strip_x_1 = strip_x_1[::2]
        strip_y_1 = strip_y_1[::2]
        for e in range(1, element_no):
            strip_x[e, :] = 0.5 * (strip_x_1[e, :] + strip_x_1[e + 1, :])
            strip_y[e, :] = 0.5 * (strip_y_1[e, :] + strip_y_1[e + 1, :])
        write_data = np.zeros((n_s, 2))
        print("elemental_disp")
        print(strip_x)
        factor = int(element_no / n_s)
        st_x = np.zeros((n_s, 1))
        st_y = np.zeros((n_s, 1))
        if int(element_no / n_s) % 2 == 0:
            for ii in range(0, n_s):
                a = int((factor / 2) - 1 + (factor * ii))
                b = int((factor / 2) + (factor * ii))
                st_x[ii, :] = 0.5 * (strip_x[a, :] + strip_x[b, :])
                st_y[ii, :] = 0.5 * (strip_y[a, :] + strip_y[b, :])
            write_data[:, 0:1] = st_x
            write_data[:, 1:2] = st_y


        else:
            strip_x = strip_x[int(0.5 * element_no / n_s)::int(element_no / n_s)]
            strip_y = strip_y[int(0.5 * element_no / n_s)::int(element_no / n_s)]

            write_data[:, 0:1] = strip_x
            write_data[:, 1:2] = strip_y
        print("write_data")    
        print(write_data)
        for j, row in enumerate(write_data, start=1):
            locals()[f"write_data_{j}"] = row
        i += 1
        write_data_1=write_data_1.reshape(1, 2)
        # print(write_data_1)
        write_data_2=write_data_2.reshape(1, 2)
        # print(write_data_2)
        write_data_3= write_data_3.reshape(1, 2)
        # print(write_data_3)
        write_data_4= write_data_4.reshape(1, 2)
        # print(write_data_4)
        write_data_5=write_data_5.reshape(1, 2)
        # print(write_data_5)
        write_data_6=write_data_6.reshape(1, 2)
        write_data_7=write_data_7.reshape(1, 2)
        write_data_8=write_data_8.reshape(1, 2)
        write_data_9=write_data_9.reshape(1, 2)
        write_data_10=write_data_10.reshape(1, 2)
        # write_data_11=write_data_11.reshape(1, 2)
        # write_data_12=write_data_12.reshape(1, 2)
        # write_data_13=write_data_13.reshape(1, 2)
        # write_data_14=write_data_14.reshape(1, 2)
        # write_data_15=write_data_15.reshape(1, 2)
        # write_data_16=write_data_16.reshape(1, 2)
        # write_data_17=write_data_17.reshape(1, 2)
        # write_data_18=write_data_18.reshape(1, 2)
        # write_data_19=write_data_19.reshape(1, 2)
        # write_data_20=write_data_20.reshape(1, 2)
    else:
        write_data_1 = np.zeros((1,2))
        write_data_2 = np.zeros((1,2))
        write_data_3 = np.zeros((1,2))
        write_data_4 = np.zeros((1,2))
        write_data_5 = np.zeros((1,2))
        write_data_6 = np.zeros((1,2))
        write_data_7 = np.zeros((1,2))
        write_data_8 = np.zeros((1,2))
        write_data_9 = np.zeros((1,2))
        write_data_10 = np.zeros((1,2))
        # write_data_11 = np.zeros((1,2))
        # write_data_12 = np.zeros((1,2))
        # write_data_13 = np.zeros((1,2))
        # write_data_14 = np.zeros((1,2))
        # write_data_15 = np.zeros((1,2))
        # write_data_16 = np.zeros((1,2))
        # write_data_17 = np.zeros((1,2))
        # write_data_18 = np.zeros((1,2))
           

    	
   

    participant.write_data(mesh_name_1, write_data_name_1, [vertex_ids_1], write_data_1)
    participant.write_data(mesh_name_2, write_data_name_2, [vertex_ids_2], write_data_2)
    participant.write_data(mesh_name_3, write_data_name_3, [vertex_ids_3], write_data_3)
    participant.write_data(mesh_name_4, write_data_name_4, [vertex_ids_4], write_data_4)
    participant.write_data(mesh_name_5, write_data_name_5, [vertex_ids_5], write_data_5)
    participant.write_data(mesh_name_6, write_data_name_6, [vertex_ids_6], write_data_6)
    participant.write_data(mesh_name_7, write_data_name_7, [vertex_ids_7], write_data_7)
    participant.write_data(mesh_name_8, write_data_name_8, [vertex_ids_8], write_data_8)
    participant.write_data(mesh_name_9, write_data_name_9, [vertex_ids_9], write_data_9)
    #participant.write_data(mesh_name_10, write_data_name_10, [vertex_ids_10], write_data_10)
    # participant.write_data(mesh_name_11, write_data_name_11, [vertex_ids_11], write_data_11)
    # participant.write_data(mesh_name_12, write_data_name_12, [vertex_ids_12], write_data_12)
    # participant.write_data(mesh_name_13, write_data_name_13, [vertex_ids_13], write_data_13)
    # participant.write_data(mesh_name_14, write_data_name_14, [vertex_ids_14], write_data_14)
    # participant.write_data(mesh_name_15, write_data_name_15, [vertex_ids_15], write_data_15)
    # participant.write_data(mesh_name_16, write_data_name_16, [vertex_ids_16], write_data_16)
    # participant.write_data(mesh_name_17, write_data_name_17, [vertex_ids_17], write_data_17)
    # participant.write_data(mesh_name_18, write_data_name_18, [vertex_ids_18], write_data_18)
    # participant.write_data(mesh_name_19, write_data_name_19, [vertex_ids_19], write_data_19)
    # participant.write_data(mesh_name_20, write_data_name_20, [vertex_ids_20], write_data_20)
    
     #########
    print("DUMMY: Advancing in time")
    participant.advance(dt)


    if participant.requires_reading_checkpoint():
        print("DUMMY: Reading iteration checkpoint")

participant.finalize()
print("DUMMY: Closing python solver dummy...")
#output
script_dir=os.path.dirname(os.path.abspath(__file__))
file_names=[os.path.join(script_dir ,'x.csv'),
            os.path.join(script_dir ,'xd.csv'),
            os.path.join(script_dir ,'xdd.csv'),
            os.path.join(script_dir ,'y.csv'),
            os.path.join(script_dir ,'yd.csv'),
            os.path.join(script_dir ,'ydd.csv')]
arrays=[x,xd,xdd,y,yd,ydd]
for file_name, array in zip(file_names, arrays):
    with open(file_name ,'w',newline='') as file:
        csv_write=csv.writer(file,delimiter=',')
        csv_write.writerows(array)
