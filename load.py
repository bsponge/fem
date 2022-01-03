import fem

def load_grid(filename):
    options = {}
    grid = fem.Grid.empty()
    with open(filename) as f:
        lines = f.readlines()
        options["SimulationTime"] = float(find_in_list("SimulationTime", lines).split(" ")[1])
        options["SimulationStepTime"] = float(find_in_list("SimulationStepTime", lines).split(" ")[1])
        options["Conductivity"] = float(find_in_list("Conductivity", lines).split(" ")[1])
        options["Alfa"] = float(find_in_list("Alfa", lines).split(" ")[1])
        options["Tot"] = float(find_in_list("Tot", lines).split(" ")[1])
        options["InitialTemp"] = float(find_in_list("InitialTemp", lines).split(" ")[1])
        options["Density"] = float(find_in_list("Density", lines).split(" ")[1])
        options["SpecificHeat"] = float(find_in_list("SpecificHeat", lines).split(" ")[1])
        options["Nodes number"] = float(find_in_list("Nodes number", lines).split(" ")[2])
        options["Elements number"] = float(find_in_list("Elements number", lines).split(" ")[2])

        idx = lines.index("*Node\n")

        nodes = []

        while '.' in lines[idx+1]:
            arr = lines[idx+1].split(" ")
            arr = list(filter(lambda x: x != '', arr))[1:]
            arr = [float(arr[0].replace(',', '')), float(arr[1].replace(',', ''))]
            nodes.append(arr)
            idx += 1

        options["*Node"] = nodes

        idx = find_index("*Element", lines)
        indices = []

        while ',' in lines[idx+1]:
            arr = lines[idx+1].split(',')[1:]
            for i in range(len(arr)):
                arr[i] = int(arr[i])-1
            indices.append(arr)
            idx += 1

        options["*Element"] = indices

        idx = find_index("*BC", lines)
        bcs = lines[idx+1].split(',')
        for i in range(len(bcs)):
            bcs[i] = int(bcs[i])-1

        options["*BC"] = bcs

    grid.nB = 0.0
    grid.nH = 0.0
    grid.H = 0.0
    grid.B = 0.0

    grid.nodes = []
    grid.elements = []

    for node in options["*Node"]:
        node = fem.Node(node[0], node[1])
        node.t = options["InitialTemp"]
        grid.nodes.append(node)

    for bc in bcs:
        grid.nodes[bc].BC = True

    for element in options["*Element"]:
        grid.elements.append(fem.Element(element))

    grid.alpha = options["Alfa"]
    grid.ro = options["Density"]
    grid.k = options["Conductivity"]
    grid.c = options["SpecificHeat"]

    grid.init_aggregated_matrices()
    grid.calculate_jacobians()
    grid.calculate_dNs()
    grid.calculate_values()
    grid.calculate_H_and_C()
    grid.load_sides()
    grid.calculate_H_BC()
    grid.calculate_P()
    grid.aggregate()
        
    return (grid, options)


def find_index(phrase, l):
    for i in range(len(l)):
        if phrase in l[i]:
            return i
    return None
        

def find_in_list(phrase, l):
    for s in l:
        if phrase in s:
            return s
    return None