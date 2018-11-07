import sys


def get_data(filename):
    with open(filename, "r") as file:
        kmers = []
        for line in file:
            kmers.append(line.strip())
    kmers.sort()
    return kmers


def prefix(text):
    return text[:-1]


def suffix(text):
    return text[1::]


def db_graph(kmers):
    adj = {}
    for kmer in kmers:
        pre = prefix(kmer)
        if pre not in adj.keys():
            adj[pre] = [suffix(kmer)]
        else:
            adj[pre].append(suffix(kmer))
    return adj


def calc_degrees(graph):
    degrees = {}
    for node, neighbors in graph.items():
        degrees[node] = (0, len(neighbors))

    for _, neighbors in graph.items():
        for node in neighbors:
            if node in degrees:
                degrees[node] = (degrees[node][0] + 1, degrees[node][1])

    return degrees


def find_start_node(degrees):
    start_node = None
    for node, degree in degrees.items():
        if degree[0] < degree[1]:
            start_node = node
    return start_node


def find_eulerian_path(node, graph, degrees, path):
    path += [node]  # == cycle.append(node)

    if node not in degrees or degrees[node][1] == 0:
        return path

    while len(graph[node]) > 0:
        temp_node = graph[node][0]
        graph[node].remove(temp_node)
        sub_path = find_eulerian_path(temp_node, graph, degrees, [])
        path = path[:1] + sub_path + path[1:]
    return path


def find_non_branching(graph, degrees):
    paths = []
    for v, _ in graph.items():
        if degrees[v] != (1,1):
            if degrees[v][1] > 0:
                for w in graph[v]:
                    edge = w
                    non_branching = v + edge[-1]
                    while edge in degrees.keys() and degrees[edge] == (1,1):
                        non_branching += graph[edge][0][-1]
                        edge = graph[edge][0]
                    paths.append(non_branching)
    return paths


def print_output(outset):
    with open("output.txt", "w") as file:
        for x in outset:
            file.write(x + ' ')
    print(outset)


if __name__ == '__main__':
    kmers = get_data(sys.argv[1])
    graph = db_graph(kmers)
    degrees = calc_degrees(graph)
    paths = find_non_branching(graph, degrees)
    print_output(paths)
