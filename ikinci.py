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
                if len(line) >= 16:
                    bus_data.append(line[:5] + line[16:].strip())

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
file_path = "ieee57cdf.txt"
bus_data, branch_data = read_data_from_file(file_path)

# Convert extracted data into matrices without spaces
bus_matrix = create_matrix(bus_data)
branch_matrix = create_matrix(branch_data)

# Example output
print("Bus Data Matrix:")
for row in bus_matrix:
    print(row)

print("\nBranch Data Matrix:")
for row in branch_matrix:
    print(row)


# Extract columns 5, 6, and 7 from bus data
ybus = []
for item in bus_data:
    row = item.split()
    ybus.append(row[4:7])

# Print Ybus Variable
print("Ybus Variable:")
for row in ybus:
    print(row)
