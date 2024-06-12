import numpy as np
import time
import matplotlib.pyplot as plt
from scipy.sparse import csr_matrix


iter = 10  #iteration
st = time.time()  # Start time
ste = time.time()  # End time, for example

def read_data_from_file(file_path):
    bus_data = []
    branch_data = []

    with open(file_path, 'r') as file:
        lines = file.readlines()
        is_bus_data = False
        is_branch_data = False

        for line in lines:
            # Check if it's the start of bus data
            if "BUS DATA FOLLOWS" in line:
                is_bus_data = True
                is_branch_data = False
                continue
            # Check if it's the start of branch data
            elif "BRANCH DATA FOLLOWS" in line:
                is_bus_data = False
                is_branch_data = True
                continue
            # Check if it's the end of bus or branch data
            elif "-999" in line:
                is_bus_data = False
                is_branch_data = False
                continue

            # Process bus data
            if is_bus_data:
                # Check if it's a valid line (excluding characters between 6 to 16)
                if len(line) >= 17:
                    bus_data.append(line[:5] + line[17:].strip())

            # Process branch data
            #elif is_branch_data:
            #    # Check if it's a valid line (excluding characters between 6 to 16)
            #    if len(line) >= 16:
            #        branch_data.append(line[:5] + line[7:].strip())
            elif is_branch_data:
                branch_data.append(line.strip())
    return bus_data, branch_data

def create_matrix(data):
    matrix = []
    for item in data:
        row = item.split()
        matrix.append(row)
    return matrix

# Path to the text file
file_path = "ieee300cdf.txt"
bus_data, branch_data = read_data_from_file(file_path)

# Convert extracted data into matrices without spaces
bus_matrix = create_matrix(bus_data)
branch_matrix = create_matrix(branch_data)
'''''
for row in branch_matrix:
    for element in row:
        # Convert the element to a string
        element_str = str(element)
        # Remove leading zeros and the decimal point
        element_str = element_str.lstrip('0').replace('.', '')
        # Check if the element has more than 6 digits
        if len(element_str) > 7:
            print("Element with more than 6 digits:", element)
            
quit()            
'''''


#print(bus_matrix[280])
#quit()
#print("Bus Data Matrix:")
#for row in bus_matrix:
#    print(row)

#print("\nBranch Data Matrix:")
#for row in branch_matrix:
#   print(row)

# Convert matrices to numpy arrays
branch_matrix_np = np.array(branch_matrix)
bus_matrix_np = np.array(bus_matrix)



frombus = branch_matrix_np[:, 0].astype(int)
tobus = branch_matrix_np[:, 1].astype(int)


r = branch_matrix_np[:, 6].astype(float)
x = branch_matrix_np[:, 7].astype(float)
b = branch_matrix_np[:, 8].astype(float)


tap = np.where(branch_matrix_np[:, 14].astype(float) == 0, 1.0, branch_matrix_np[:, 14].astype(float)) #doğru
phase_shift = branch_matrix_np[:, 15].astype(float)  # Assuming column 15 contains phase shift angles in degrees
#print(tap)
#quit()
shuntdataG = bus_matrix_np[:, 14].astype(float)  ####doğru
shuntdataB = bus_matrix_np[:, 15].astype(float)  ####doğru

#print(shuntdataG)
#print(shuntdataB)


nb = len(bus_matrix_np)     # number of buses = 300
nl = len(frombus)           # number of branches (line) = 411

bus_indices = {}
for i, bus_number in enumerate(bus_matrix_np[:, 0].astype(int), start=1):
    bus_indices[bus_number] = i - 1  # Adjust for 0-based indexing




i = 1j
z = r + x * i
y = 1. / z   #series admitance of the lines
chargingSusceptance = b  ## total line charging = B, B/2 yi mi kullanacağız??

# Initialize the yBus and yShunt matrices
yBus = np.zeros((nb, nb), dtype=complex)
yShunt = np.zeros((nb, nb), dtype=complex)

# Fill the yBus and yShunt matrices
for idx in range(nl):
    from_bus_number = frombus[idx]
    to_bus_number = tobus[idx]
    m = bus_indices[from_bus_number]
    n = bus_indices[to_bus_number]
    a = tap[idx]
    
    theta = np.deg2rad(phase_shift[idx])
    y_series = y[idx]
    y_shunt = chargingSusceptance[idx]*1j

    yBus[m, m] += y_series / a**2 + y_shunt / 2
    yBus[m, n] -= y_series / a
    yBus[n, m] -= y_series / a
    yBus[n, n] += y_series + y_shunt / 2

    yShunt[m, n] += y_shunt / 2 + y_series * (1 - a) / a**2
    yShunt[n, m] += y_shunt / 2 + y_series * (a - 1) / a

# Include shunt data in the yBus matrix
for i in range(nb):
    yBus[i, i] += shuntdataG[i] + shuntdataB[i] * 1j


# Extract the real and imaginary parts
G = np.real(yBus)
B = np.imag(yBus)

yBus_sparse = csr_matrix(yBus)

# Plotting the sparsity pattern
plt.figure(figsize=(10, 8))
plt.spy(yBus_sparse, markersize=1)
plt.title('Sparsity Pattern of yBus')
plt.xlabel('Bus Number')
plt.ylabel('Bus Number')
plt.show()
#quit()

# Example output
#print("Ybus Matrix:")
#print(yBus)    #57*57 matrice
#quit()
#print("G Matrix (Real part of Ybus):")
#print(G)       #57*57 matrice
#print("B Matrix (Imaginary part of Ybus):")
#print(B)       #57*57 matrice(checked)


#şuanda değişkenler yBus içinde G ve B, nb ve nl uzunluklar, shuntdataG, shuntdataB, bus_matrix_np, branch_matrix_np
# P, Q, V, theta ve bus_type empty variables inside def power flow.

def power_flow_newton_raphson(bus_matrix_np, G, B, yBus, tol=1e-6, max_iter=20):
  
    P = np.zeros(nb)
    Q = np.zeros(nb)
    V = np.ones(nb)  # Flat start with voltage magnitudes of 1 per unit
    theta = np.zeros(nb)  # Flat start with phase angles of 0 degrees
    
    # Identifying bus types
    bus_type = np.zeros(nb)  # 0 = PQ_bus_idx, 1 = PQ_bus_idx, 2 = PV, 3 = Slack
    for idx, bus in enumerate(bus_matrix_np):
        bus_type[idx] = int(bus[3])  # Assuming bus type information is in the fourth column

        # Assign initial values of P and Q from CDF file
        P_generation = float(bus[8])  # Real power generation from CDF file
        P_load = float(bus[6])  # Real power load from CDF file
        Q_generation = float(bus[9])  # Reactive power generation from CDF file
        Q_load = float(bus[7])  # Reactive power load from CDF file
        P[idx] = P_generation - P_load  # Injected real power
        Q[idx] = Q_generation - Q_load  # Injected reactive power   # P VE Q DEĞERLERİ BURAYA KADAR DOĞRU
    
    

    # Identify the slack, PQ_bus_idx, and PV buses based on the bus type information
    slack_bus_idx = np.where(bus_type == 3)[0][0]  # Assuming there's one slack bus
    PQ_bus_idx = np.where(bus_type == 0)[0]  # Assuming bus type 0 corresponds to PQ_bus_idx bus
    PV_bus_idx = np.where(bus_type == 2)[0]  # Assuming bus type 1 corresponds to PV bus

    # Calculate the number of slack, PQ_bus_idx, and PV buses
    nslack = 1  # Assuming there's one slack bus
    nPQ = len(PQ_bus_idx)
    nPV = len(PV_bus_idx)           #PV BUS AND PQ_bus_idx BUS NUMBERS ARE CORRECT.

    
    
    for iteration in range(max_iter):
        # Calculate power injections
        P_calc = np.zeros(nb)
        Q_calc = np.zeros(nb)
        for k in range(nb):
            for m in range(nb):
                P_calc[k] += V[k] * V[m] * (G[k, m] * np.cos(theta[k] - theta[m]) + B[k, m] * np.sin(theta[k] - theta[m]))
                Q_calc[k] += V[k] * V[m] * (G[k, m] * np.sin(theta[k] - theta[m]) - B[k, m] * np.cos(theta[k] - theta[m]))
        
        
        dP = P - P_calc
        dQ = Q - Q_calc
        
        # Build the mismatch Q vector
        mismatch_Q = np.zeros(nPQ)
        for k in range(nPQ):
            n = PQ_bus_idx[k]
            mismatch_Q[k] = dQ[n]

        # Build the mismatch P vector non-including the slack bus
        mismatch_P = dP[1:nb]
        
        # Mismatch vector [dP; dQ]
        mismatch = np.concatenate((mismatch_P, mismatch_Q))  #mismatch_PQ yu mismatch yaptım.
        
        #print(len(mismatch))        ##### 5 doğru
        #quit()              

        # Apply tolerance check
        if np.max(np.abs(mismatch)) < tol:
            break
        
        # Initialize Jacobian matrix
        J1 = np.zeros((nb-1, nb-1))
        J2 = np.zeros((nb-1, nPQ))
        J3 = np.zeros((nPQ, nb-1))
        J4 = np.zeros((nPQ, nPQ))
        
        # J1 [dPk/dTheta_j]
        for k in range(1, nb):  # k, position of the row, slack bus is not included so the for is initialized from 2
            for j in range(1, nb):  # j, column position
                if j == k:  # Diagonal element
                    for m in range(nb):  # Sum of each element
                        J1[k-1, j-1] += V[k] * V[m] * (-G[k, m] * np.sin(theta[k] - theta[m]) + B[k, m] * np.cos(theta[k] - theta[m]))
                    J1[k-1, j-1] -= V[k]**2 * B[k, k]  # Subtracts the element when m = k, since its derivative is zero
                else:  # Non-diagonal element
                    J1[k-1, j-1] = V[k] * V[j] * (G[k, j] * np.sin(theta[k] - theta[j]) - B[k, j] * np.cos(theta[k] - theta[j]))

        # J2 [dPk/dV_j]
        for k in range(1, nb):  # k, position of the row, slack bus is not included so the for is initialized from 2
            for j in range(nPQ):  # j, column position
                m = PQ_bus_idx[j]
                if m == k:  # Diagonal element
                    for m in range(nb):  # Sum of each element
                        J2[k-1, j] += V[m] * (G[k, m] * np.cos(theta[k] - theta[m]) + B[k, m] * np.sin(theta[k] - theta[m]))
                    J2[k-1, j] += V[k] * G[k, k]  # Add an element when m = k, when m=k its derivative contains a term of 2V(k)
                else:  # Non-diagonal element
                    J2[k-1, j] = V[k] * (G[k, m] * np.cos(theta[k] - theta[m]) + B[k, m] * np.sin(theta[k] - theta[m]))

        # J3 [dQk/dTheta_j]
        for k in range(nPQ):  # k, position of the row, slack bus is not included so the for is initialized from 2
            n = PQ_bus_idx[k]
            for j in range(1, nb):  # j, column position
                if j == n:  # Diagonal element
                    for m in range(nb):  # Sum of each element
                        J3[k, j-1] += V[n] * V[m] * (G[n, m] * np.cos(theta[n] - theta[m]) + B[n, m] * np.sin(theta[n] - theta[m]))
                    J3[k, j-1] -= V[n]**2 * G[n, n]  # Add an element when m = k, when m=k its derivative contains a term of 2V(k)
                else:  # Non-diagonal element
                    J3[k, j-1] = V[n] * V[j] * (-G[n, j] * np.cos(theta[n] - theta[j]) - B[n, j] * np.sin(theta[n] - theta[j]))

        # J4 [dQk/dV_j]
        for k in range(nPQ):  # k, position of the row, slack bus is not included so the for is initialized from 2
            n = PQ_bus_idx[k]
            for j in range(nPQ):  # j, column position
                m = PQ_bus_idx[j]
                if m == n:  # Diagonal element
                    for m in range(nb):  # Sum of each element
                        J4[k, j] += V[m] * (G[n, m] * np.sin(theta[n] - theta[m]) - B[n, m] * np.cos(theta[n] - theta[m]))
                    J4[k, j] -= V[n] * B[n, n]  # Add an element when m = k, when m=k its derivative contains a term of 2V(k)
                else:  # Non-diagonal element
                    J4[k, j] = V[n] * (G[n, m] * np.sin(theta[n] - theta[m]) - B[n, m] * np.cos(theta[n] - theta[m]))

        J = np.vstack((np.hstack((J1, J2)), np.hstack((J3, J4))))
        
        #print(J)       ######## test amaçlı
        #quit()         ######## Jacobian doğru !!!!
    
        # Solve for delta
        #delta = np.linalg.solve(J, mismatch)
        
        # Update voltage magnitudes and phase angles
        #theta += delta[nb-1]
        #V += delta[nb-1]

        mismatch_X = np.linalg.inv(J) @ mismatch

        # Extract dTeth and dV
        dTeth = mismatch_X[:nb-1]
        dV = mismatch_X[nb-1:]              #buraya -1 koyunca çözüldü dimension problemi.

        # Update Theta
        theta[1:nb] += dTeth
        
        # Update V for PQ buses
        for i in range(nPQ):
            n = PQ_bus_idx[i-1]
            V[n] += V[i-1]
    
    return V, theta



V, theta = power_flow_newton_raphson(bus_matrix_np, G, B, yBus)

#print("Voltage Magnitudes:", np.abs(V))
#print("Voltage Angles:", (180 / np.pi) * theta)
#quit()

# Create a mapping from bus number to index
bus_numbers = bus_matrix_np[:, 0].astype(int)
bus_index_map = {bus_num: idx for idx, bus_num in enumerate(bus_numbers)}

# Map the frombus and tobus arrays to indices
frombus_idx = np.array([bus_index_map[bus] for bus in frombus])
tobus_idx = np.array([bus_index_map[bus] for bus in tobus])

# Ensure y and chargingSusceptance have the correct shapes
y = np.array([complex(r[i], x[i]) for i in range(len(r))])
chargingSusceptance = np.array([complex(0, b[i]) for i in range(len(b))])

# Calculate complex power flows
Ss = V[frombus_idx] * np.conj((V[frombus_idx] - V[tobus_idx]) * y + V[frombus_idx] * chargingSusceptance)
Sr = V[tobus_idx] * np.conj((V[tobus_idx] - V[frombus_idx]) * y + V[tobus_idx] * chargingSusceptance)

# Define active and reactive power flows
Pij = np.real(Ss)
Qij = np.imag(Ss)
Pji = np.real(Sr)
Qji = np.imag(Sr)

# Active power losses
P_loss = np.sum(Pij + Pji)

# Reactive power losses
Q_loss = np.sum(Qij + Qji)

#print(P_loss)



print('                      Newton Raphson Load-Flow Study')
print('                    Report of Power Flow Calculations ')
print(' ')
print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
print(f'Number of iterations       : {iter}')
print(f'Solution time              : {time.time() - st} sec.')
print(f'Total real power losses    : {P_loss}.')
print(f'Total reactive power losses: {Q_loss}.\n')
#print('                                      Generation             Load')
#print('      Bus      Volts     Angle      Real  Reactive      Real  Reactive ')
print('      Bus      Volts     Angle     ')
abcdef=bus_matrix_np[:, 0].astype(int)
abcdef=np.array(abcdef)
ywz = np.column_stack((abcdef.astype(int), np.abs(V), (180 / np.pi) * theta,))
print(ywz)

print('                      Line Flows                     ')
print('    #Line    From Bus   To Bus     Real    Reactive   ')
l = np.arange(1, len(frombus) + 1)
xy = np.column_stack((l, frombus, tobus, Pij, Qij))
yx = np.column_stack((l, tobus, frombus, Pji, Qji))
print(xy)
print(yx)


#buses_at_limits = []
# Check which buses are at their limits
#for i in range(len(P)):
#    if P[i] >= Pmax[i] or P[i] <= Pmin[i] or Q[i] >= Qmax[i] or Q[i] <= Qmin[i]:
#        buses_at_limits.append(i + 1)  # Append bus index (adding 1 to convert to 1-based index)

# Print buses at limits
#print(f"Buses at PQ limits: {buses_at_limits}")