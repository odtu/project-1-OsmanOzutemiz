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


#şuanda değişkenler yBus içinde G ve B, nb ve nl uzunluklar, shuntdataG, shuntdataB, bus_matrix_np, branch_matrix_np
# P, Q, V, theta ve bus_type empty variables inside def power flow.

P = np.zeros(nb)
Q = np.zeros(nb)
bus_type = np.zeros(nb)  # 0 = PQ, 1 = PQ, 2 = PV, 3 = Slack


for idx, bus in enumerate(bus_matrix_np):
    bus_type[idx] = int(bus[3])  # Assuming bus bus_type information is in the fourth column

    # Assign initial values of P and Q from CDF file
    P_generation = float(bus[8])  # Real power generation from CDF file
    P_load = float(bus[6])  # Real power load from CDF file
    Q_generation = float(bus[9])  # Reactive power generation from CDF file
    Q_load = float(bus[7])  # Reactive power load from CDF file
    P[idx] = P_generation - P_load  # Injected real power
    Q[idx] = Q_generation - Q_load  # Injected reactive power   # P VE Q DEĞERLERİ BURAYA KADAR DOĞRU



# Extract data from busdata
bus = bus_matrix_np[:, 0].astype(int) #1 den 14 e
V = np.ones(nb)
a = tap


#print(bus)
#quit()

# Initialize variables
Vprev = V.copy()
P = (P_generation - P_load) / 100
Q = (Q_generation - Q_load) / 100

P = np.array(P)
Q = np.array(Q)

print(P)
quit()
toler = 1
iteration = 0
A = 1.6

# Newton-Raphson iteration
while toler > 0.001:
    for i in range(0, nb):
        sumyv = 0
        for k in range(nb):
            if i != k:
                sumyv += yBus[i, k] * V[k]

        if bus_type[i] == 2:  # PV bus
            Q[i] = -np.imag(np.conj(V[i]) * (sumyv + yBus[i, i] * V[i]))
            bus_matrix_np[i, 6] = Q[i]  # Update reactive power in busdata

        V[i] = (1 / yBus[i, i]) * ((P[i] - 1j * Q[i]) / np.conj(V[i]) - sumyv)

        if bus_type[i] == 2:  # PV bus
            Vprev[i] = abs(Vprev[i]) * (V[i] / abs(V[i]))
            bus_matrix_np[i, 1] = Vprev[i]
            V[i] = Vprev[i]

        if bus_type[i] == 3:  # Slack bus
            Vacc = Vprev[i] + A * (V[i] - Vprev[i])
            bus_matrix_np[i, 1] = Vacc
            V[i] = Vacc

    iteration += 1
    toler = np.max(np.abs(np.abs(V) - np.abs(Vprev)))
    Vprev = V.copy()

# Results
Q_result = Q
iteration_result = iteration
sumyv_result = sumyv
V_result = V

# Output z
z = np.zeros((nb, 3))
z[:, 0] = busdata[:, 0]
z[:, 2] = np.angle(V) * (180 / np.pi)
z[:, 1] = np.abs(V)

print(z)
