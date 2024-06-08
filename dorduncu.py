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


nb = len(bus_matrix_np)     # number of buses = 57
nl = len(frombus)           # number of branches (line) = 80

i = 1j
z = r + x * i
y = 1. / z   #series admitance of the lines
chargingSusceptance = b  ## total line charging = B, B/2 yi mi kullanacağız??

# Initialize the yBus and yShunt matrices
yBus = np.zeros((nb, nb), dtype=complex)
yShunt = np.zeros((nb, nb), dtype=complex)

# Fill the yBus and yShunt matrices
for idx in range(nl):
    m = frombus[idx] - 1  # Adjust for 0-based indexing
    n = tobus[idx] - 1  # Adjust for 0-based indexing
    a = tap[idx]
    theta = np.deg2rad(phase_shift[idx])
    y_series = y[idx]
    y_shunt = chargingSusceptance[idx]

    yBus[m, m] += y_series / a**2 + y_shunt / 2
    yBus[m, n] -= y_series / a * np.exp(-i * theta)
    yBus[n, m] -= y_series / a * np.exp(i * theta)
    yBus[n, n] += y_series + y_shunt / 2

    yShunt[m, n] += y_shunt / 2 + y_series * (1 - a) / a**2
    yShunt[n, m] += y_shunt / 2 + y_series * (a - 1) / a

# Include shunt data in the yBus matrix
for i in range(nb):
    yBus[i, i] += shuntdataG[i] + shuntdataB[i] * 1j


# Extract the real and imaginary parts
G = np.real(yBus)
B = np.imag(yBus)

# Example output
#print("Ybus Matrix:")
#print(yBus)    #57*57 matrice
#print("G Matrix (Real part of Ybus):")
#print(G)       #57*57 matrice
#print("B Matrix (Imaginary part of Ybus):")
#print(B)       #57*57 matrice(checked)




def power_flow_newton_raphson(bus_matrix_np, branch_matrix_np, yBus, tol=1e-6, max_iter=10):
    nb = len(bus_matrix_np)
    P = np.zeros(nb)
    Q = np.zeros(nb)
    V = np.ones(nb)  # Flat start with voltage magnitudes of 1 per unit
    theta = np.zeros(nb)  # Flat start with phase angles of 0 degrees
    
    # Identifying bus types
    bus_type = np.zeros(nb)  # 0 = PQ, 1 = PQ, 2 = PV, 3 = Slack
    for idx, bus in enumerate(bus_matrix_np):
        bus_type[idx] = int(bus[3])  # Assuming bus type information is in the fourth column

        # Assign initial values of P and Q from CDF file
        P_generation = float(bus[8])  # Real power generation from CDF file
        P_load = float(bus[6])  # Real power load from CDF file
        Q_generation = float(bus[9])  # Reactive power generation from CDF file
        Q_load = float(bus[7])  # Reactive power load from CDF file
        P[idx] = P_generation - P_load  # Injected real power
        Q[idx] = Q_generation - Q_load  # Injected reactive power
    
    slack_bus_idx = np.where(bus_type == 3)[0][0]  # Assuming there's one slack bus
    
    for iteration in range(max_iter):
        # Calculate power injections
        P_calc = np.zeros(nb)
        Q_calc = np.zeros(nb)
        for i in range(nb):
            for k in range(nb):
                P_calc[i] += V[i] * V[k] * (yBus[i, k].real * np.cos(theta[i] - theta[k]) + yBus[i, k].imag * np.sin(theta[i] - theta[k]))
                Q_calc[i] += V[i] * V[k] * (yBus[i, k].real * np.sin(theta[i] - theta[k]) - yBus[i, k].imag * np.cos(theta[i] - theta[k]))
        
        # Compute mismatches
        P_mismatch = P - P_calc
        Q_mismatch = Q - Q_calc
        mismatch = np.concatenate((P_mismatch, Q_mismatch))

        #print(mismatch)
        #quit()

        # Apply tolerance check
        if np.max(np.abs(mismatch)) < tol:
            break
        
        # Initialize Jacobian matrix
        J = np.zeros((2*nb, 2*nb))
        for i in range(nb):
            for k in range(nb):
                # Fill Jacobian matrix
                if i == k:
                    J[i, k] = -Q_calc[i] - (V[i]**2 * yBus[i, i].imag)
                    J[i, k + nb] = P_calc[i] / V[i] + V[i] * yBus[i, i].real
                    J[i + nb, k] = P_calc[i] - (V[i]**2 * yBus[i, i].real)
                    J[i + nb, k + nb] = Q_calc[i] / V[i] - V[i] * yBus[i, i].imag
                else:
                    J[i, k] = V[i] * V[k] * (yBus[i, k].imag * np.cos(theta[i] - theta[k]) - yBus[i, k].real * np.sin(theta[i] - theta[k]))
                    J[i, k + nb] = V[i] * (yBus[i, k].real * np.cos(theta[i] - theta[k]) + yBus[i, k].imag * np.sin(theta[i] - theta[k]))
                    J[i + nb, k] = -V[i] * V[k] * (yBus[i, k].real * np.cos(theta[i] - theta[k]) + yBus[i, k].imag * np.sin(theta[i] - theta[k]))
                    J[i + nb, k + nb] = V[i] * (yBus[i, k].real * np.sin(theta[i] - theta[k]) - yBus[i, k].imag * np.cos(theta[i] - theta[k]))
        
        #print(J[0]) test amaçlı
        #quit()
        
        # Solve for delta
        delta = np.linalg.solve(J, mismatch)
        
        # Update voltage magnitudes and phase angles
        theta += delta[:nb]
        V += delta[nb:]
    
    return V, theta



V, theta = power_flow_newton_raphson(bus_matrix_np, branch_matrix_np, yBus)

print("Voltage Magnitudes:", V)
print("Voltage Angles:", theta)
