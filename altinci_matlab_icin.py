import numpy as np
import sys

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
            elif is_branch_data:
                # Check if it's a valid line (excluding characters between 6 to 16)
                if len(line) >= 16:
                    branch_data.append(line[:5] + line[7:].strip())

    return bus_data, branch_data

def create_matrix(data):
    matrix = []
    for item in data:
        row = item.split()
        matrix.append(row)
    return matrix

# Path to the text file
file_path = "ieee14cdf.txt"
bus_data, branch_data = read_data_from_file(file_path)

# Convert extracted data into matrices without spaces
bus_matrix = create_matrix(bus_data)
branch_matrix = create_matrix(branch_data)


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


nbuses = len(bus_matrix_np)     # number of buses = 57
nl = len(frombus)           # number of branches (line) = 80

i = 1j
z = r + x * i
y = 1. / z   #series admitance of the lines
chargingSusceptance = b  ## total line charging = B, B/2 yi mi kullanacağız??

# Initialize the Y and yShunt matrices
Y = np.zeros((nbuses, nbuses), dtype=complex)
yShunt = np.zeros((nbuses, nbuses), dtype=complex)

# Fill the Y and yShunt matrices
for idx in range(nl):
    m = frombus[idx] - 1  # Adjust for 0-based indexing
    n = tobus[idx] - 1  # Adjust for 0-based indexing
    a = tap[idx]
    Theta = np.deg2rad(phase_shift[idx])
    y_series = y[idx]
    y_shunt = chargingSusceptance[idx]

    Y[m, m] += y_series / a**2 + y_shunt / 2
    Y[m, n] -= y_series / a * np.exp(-i * Theta)
    Y[n, m] -= y_series / a * np.exp(i * Theta)
    Y[n, n] += y_series + y_shunt / 2

    yShunt[m, n] += y_shunt / 2 + y_series * (1 - a) / a**2
    yShunt[n, m] += y_shunt / 2 + y_series * (a - 1) / a

# Include shunt data in the Y matrix
for i in range(nbuses):
    Y[i, i] += shuntdataG[i] + shuntdataB[i] * 1j


# Extract the real and imaginary parts
G = np.real(Y)
B = np.imag(Y)

# Example output
#print("Ybus Matrix:")
#print(Y)    #57*57 matrice
#print("G Matrix (Real part of Ybus):")
#print(G)       #57*57 matrice
#print("B Matrix (Imaginary part of Ybus):")
#print(B)       #57*57 matrice(checked)


#şuanda değişkenler Y içinde G ve B, nbuses ve nl uzunluklar, shuntdataG, shuntdataB, bus_matrix_np, branch_matrix_np
# P_neta, Q_neta, V, Theta ve bus_type empty variables inside def power flow.

def power_flow_newton_raphson(bus_matrix_np, branch_matrix_np, Y, G, B, tol=1e-6, max_iter=100):
  
    P_neta = np.zeros(nbuses)
    Q_neta = np.zeros(nbuses)
    V = np.ones(nbuses)  # Flat start with voltage magnitudes of 1 per unit
    Theta = np.zeros(nbuses)  # Flat start with phase angles of 0 degrees
    
    # Identifying bus types
    bus_type = np.zeros(nbuses)  # 0 = PQ, 1 = PQ, 2 = PV, 3 = Slack
    for idx, bus in enumerate(bus_matrix_np):
        bus_type[idx] = int(bus[3])  # Assuming bus type information is in the fourth column

        # Assign initial values of P_neta and Q_neta from CDF file
        P_generation = float(bus[8])  # Real power generation from CDF file
        P_load = float(bus[6])  # Real power load from CDF file
        Q_generation = float(bus[9])  # Reactive power generation from CDF file
        Q_load = float(bus[7])  # Reactive power load from CDF file
        P_neta[idx] = P_generation - P_load  # Injected real power
        Q_neta[idx] = Q_generation - Q_load  # Injected reactive power   # P_neta VE Q_neta DEĞERLERİ BURAYA KADAR DOĞRU
    
    

    # Identify the slack, PQ, and PV buses based on the bus type information
    slack_bus_idx = np.where(bus_type == 3)[0][0]  # Assuming there's one slack bus
    PQ = np.where(bus_type == 0)[0]  # Assuming bus type 0 corresponds to PQ bus
    PV = np.where(bus_type == 2)[0]  # Assuming bus type 1 corresponds to PV bus

    # Calculate the number of slack, PQ, and PV buses
    nslack = 1  # Assuming there's one slack bus
    nPQ = len(PQ)
    nPV = len(PV)           #PV BUS AND PQ BUS NUMBERS BURAYA KADAR DOĞRU

    
    
    for iteration in range(max_iter):
        # Calculate power injections
        P_cal = np.zeros(nbuses)
        Q_cal = np.zeros(nbuses)
        for k in range(nbuses):
            for m in range(nbuses):
                P_cal[k] += V[k] * V[m] * (G[k, m] * np.cos(Theta[k] - Theta[m]) + B[k, m] * np.sin(Theta[k] - Theta[m]))
                Q_cal[k] += V[k] * V[m] * (G[k, m] * np.sin(Theta[k] - Theta[m]) - B[k, m] * np.cos(Theta[k] - Theta[m]))
        
        
        dP = P_neta - P_cal
        dQ = Q_neta - Q_cal
        
        # Build the mismatch_PQ Q_neta vector
        mismatch_Q = np.zeros(nPQ)  
        for k in range(nPQ):
            n = PQ[k]
            mismatch_Q[k] = dQ[n]

        # Build the mismatch_PQ P_neta vector non-including the slack bus
        mismatch_P = dP[1:nbuses]
        
        # Mismatch vector [dP; dQ]
        mismatch_PQ = np.concatenate((mismatch_P, mismatch_Q))  #mismatch_PQ yu mismatch_PQ yaptım.

        # Apply tolerance check
        if np.max(np.abs(mismatch_PQ)) < tol:
            break
        
        # Initialize Jacobian matrix
        # J = np.zeros((2*nbuses, 2*nbuses))
        
        J1 = np.zeros((nbuses-1, nbuses-1))
        J2 = np.zeros((nbuses-1, nPQ))
        J3 = np.zeros((nPQ, nbuses-1))
        J4 = np.zeros((nPQ, nPQ))
        
        # J1 [dPk/dTheta_j]
        for k in range(1, nbuses):  # k, position of the row, slack bus is not included so the for is initialized from 2
            for j in range(1, nbuses):  # j, column position
                if j == k:  # Diagonal element
                    for m in range(nbuses):  # Sum of each element
                        J1[k-1, j-1] += V[k] * V[m] * (-G[k, m] * np.sin(Theta[k] - Theta[m]) + B[k, m] * np.cos(Theta[k] - Theta[m]))
                    J1[k-1, j-1] -= V[k]**2 * B[k, k]  # Subtracts the element when m = k, since its derivative is zero
                else:  # Non-diagonal element
                    J1[k-1, j-1] = V[k] * V[j] * (G[k, j] * np.sin(Theta[k] - Theta[j]) - B[k, j] * np.cos(Theta[k] - Theta[j]))

        # J2 [dPk/dV_j]
        for k in range(1, nbuses):  # k, position of the row, slack bus is not included so the for is initialized from 2
            for j in range(nPQ):  # j, column position
                m = PQ[j]
                if m == k:  # Diagonal element
                    for m in range(nbuses):  # Sum of each element
                        J2[k-1, j] += V[m] * (G[k, m] * np.cos(Theta[k] - Theta[m]) + B[k, m] * np.sin(Theta[k] - Theta[m]))
                    J2[k-1, j] += V[k] * G[k, k]  # Add an element when m = k, when m=k its derivative contains a term of 2V(k)
                else:  # Non-diagonal element
                    J2[k-1, j] = V[k] * (G[k, m] * np.cos(Theta[k] - Theta[m]) + B[k, m] * np.sin(Theta[k] - Theta[m]))

        # J3 [dQk/dTheta_j]
        for k in range(nPQ):  # k, position of the row, slack bus is not included so the for is initialized from 2
            n = PQ[k]
            for j in range(1, nbuses):  # j, column position
                if j == n:  # Diagonal element
                    for m in range(nbuses):  # Sum of each element
                        J3[k, j-1] += V[n] * V[m] * (G[n, m] * np.cos(Theta[n] - Theta[m]) + B[n, m] * np.sin(Theta[n] - Theta[m]))
                    J3[k, j-1] -= V[n]**2 * G[n, n]  # Add an element when m = k, when m=k its derivative contains a term of 2V(k)
                else:  # Non-diagonal element
                    J3[k, j-1] = V[n] * V[j] * (-G[n, j] * np.cos(Theta[n] - Theta[j]) - B[n, j] * np.sin(Theta[n] - Theta[j]))

        # J4 [dQk/dV_j]
        for k in range(nPQ):  # k, position of the row, slack bus is not included so the for is initialized from 2
            n = PQ[k]
            for j in range(nPQ):  # j, column position
                m = PQ[j]
                if m == n:  # Diagonal element
                    for m in range(nbuses):  # Sum of each element
                        J4[k, j] += V[m] * (G[n, m] * np.sin(Theta[n] - Theta[m]) - B[n, m] * np.cos(Theta[n] - Theta[m]))
                    J4[k, j] -= V[n] * B[n, n]  # Add an element when m = k, when m=k its derivative contains a term of 2V(k)
                else:  # Non-diagonal element
                    J4[k, j] = V[n] * (G[n, m] * np.sin(Theta[n] - Theta[m]) - B[n, m] * np.cos(Theta[n] - Theta[m]))

        J = np.vstack((np.hstack((J1, J2)), np.hstack((J3, J4))))
        
        #print(J[0]) test amaçlı
        #quit()
        
        # Solve for delta
        ''''
        delta = np.linalg.solve(J, mismatch_PQ)
     
        # Update voltage magnitudes and phase angles
        Theta += delta[:nbuses]
        V += delta[nbuses:]
        '''
        
        mismatch_X = np.linalg.inv(J) @ mismatch_PQ

        # Extract dTeth and dV
        dTeth = mismatch_X[:nbuses-1]
        dV = mismatch_X[nbuses-1:]              #buraya -1 koyunca çözüldü dimension problemi.

        # Update Theta
        Theta[1:nbuses] += dTeth

        # Update V for PQ buses
        for i in range(nPQ):
            n = PQ[i-1]
            V[n] += dV[i-1]
        
        
        # Complex voltage in per unit (PU)
        Vm = V * np.cos(Theta) + 1j * V * np.sin(Theta)

        # Currents at each node
        I = np.dot(Y, Vm)

        # Power at each node
        S = Vm * np.conj(I)

    return V, Theta, S, I



V, Theta, S, I = power_flow_newton_raphson(bus_matrix_np, branch_matrix_np, Y, G, B)


print("Voltage Magnitudes:", V)
print("Voltage Angles:", Theta)
