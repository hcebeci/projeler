import numpy as np
from scipy.stats import binom
import random
import time

def duzenleme(sequence):
    vertecies_degree=dict()
    vertecies_adjacent=dict()
    for i,j in enumerate(sequence):
        vertecies_degree[i]=j
        vertecies_adjacent[i]=[]
    return vertecies_degree,vertecies_adjacent

def check_if_graph(data):
    data.sort(reverse=True)
    data.insert(0,0)
    graph = True
    for x in range(1,len(data)):
        if (sum(data) % 2) == 1:
            graph = False
            break
        first_value = sum(data[1:x+1])
        second_value = x * (x - 1)
        for y in range(x + 1, len(data)):
            second_value += min(x, data[y])
        if first_value <= second_value:
            graph = True
        else:
            graph = False
            break
    data.pop(0)
    return graph

def HH_high_to_high(degrees,adjacent):
    vertex=max(degrees,key=degrees.get)
    if degrees[vertex]==0:  #if all verticies have reached its limit, STOP
        return adjacent
    v=degrees.copy()
    v.pop(vertex)
    d= sorted(v.items(), key=lambda x: x[1], reverse=True)

    #adding edges to the highest degree vertex
    for i in range(degrees[vertex]):
        x=d[i][0]
        degrees[x]-=1
        degrees[vertex] -= 1
        adjacent[x].append(vertex)
        adjacent[vertex].append(x)

    return  HH_high_to_high(degrees,adjacent)

def HH_random_to_high(degrees,adjacent):
    k = len(degrees)
    verticies = random.sample(range(k),k)

    for vertex in verticies:
        v = degrees.copy()
        v.pop(vertex)   #Selected vertex will not match with itself
        d = sorted(v.items(), key=lambda x: x[1], reverse=True)
        for i in range(degrees[vertex]):        #determining the adjacent verticies and decreasing the degrees for the next sequence.
            x=d[i][0]
            degrees[x]-=1
            degrees[vertex] -= 1
            adjacent[x].append(vertex)
            adjacent[vertex].append(x)
    return adjacent

def HH_small_to_high(degrees,adjacent):
    temp = list()
    for values in degrees:
        if degrees[values] == 0:
            temp.append(values)
    for p in temp:
        degrees.pop(p)    #deletes the values which are 0 to find the other min values other than 0.
    if len(degrees) == 0:   #If all the verticies have reached its degree limit STOP and return adjacent list
        return adjacent

    vertex = min(degrees, key=degrees.get)

    v = degrees.copy()
    v.pop(vertex)
    d = sorted(v.items(), key=lambda x: x[1], reverse=True)

    for i in range(degrees[vertex]):
        x = d[i][0]
        degrees[x] -= 1
        degrees[vertex] -= 1
        adjacent[x].append(vertex)
        adjacent[vertex].append(x)

    return HH_small_to_high(degrees, adjacent)

def HH_small_to_small(degrees,adjacent):
    temp = list()
    for values in degrees:
        if degrees[values] == 0:
            temp.append(values)  #finding the nodes that achieved its limit
    for p in temp:
        degrees.pop(p)   #removing the nodes that achieved its limit
    if len(degrees) == 0:   #TRUE when all nodes achived its limit
        return adjacent

    vertex = min(degrees, key=degrees.get)
    v = degrees.copy()
    v.pop(vertex)
    d = sorted(v.items(), key=lambda x:x[1])

    #determining the adjacent verticies and decreasing the degrees for the next sequence.
    for i in range(degrees[vertex]):
        x = d[i][0]
        degrees[x] -= 1
        degrees[vertex] -= 1
        adjacent[x].append(vertex)
        adjacent[vertex].append(x)

    return HH_small_to_small(degrees, adjacent)

def DFS_connected(node,verticies_adjacent,visited):
    if node not in visited:     #adds the starting node to the list
        visited.append(node)
    for i in verticies_adjacent[node]:      #checks all the possible ways
        if i not in visited:
            visited.append(i)       #if it is not visited adds it to the visited list
            DFS_connected(i, verticies_adjacent, visited)       #explores the new path
    if len(visited) == len(verticies_adjacent.keys()):
        visited.sort()
        return True , visited
    visited.sort()
    return False , visited

def subgraph(verticies_degree,verticies_adjacent):
    verticieslist = list(verticies_degree.keys())
    total = list()
    start = 0
    n_sub = 0
    graph_sub = list()
    while total != verticieslist:
        n_sub += 1
        tf, visited = DFS_connected(start, verticies_adjacent, [])  #returns the visited list which is a subraph
        total += visited    #check the total visited verticies
        graph_sub.append(visited)   #adds the subgraphs
        total.sort()
        for x in verticieslist:
            if x not in total:  #if x not visited, starts from x and finds subgraph starting from x
                start = x
                break
    return graph_sub

def make_connected(verticies_degree,verticies_adjacent):
    graph_sub = subgraph(verticies_degree,verticies_adjacent)
    changable = list()
    for i in range(len(graph_sub)):
        for y in range(len(graph_sub[i])):
            for z in range(y+1,len(graph_sub[i])):
                y_node = graph_sub[i][y]
                z_node = graph_sub[i][z]
                y_list = verticies_adjacent[y_node]
                z_list = verticies_adjacent[z_node]
                if y_node in z_list and z_node in y_list:   #checks if they are connected

                    y_loc_z = z_list.index(y_node)          #finds the location in adjacent list
                    z_loc_y = y_list.index(z_node)
                    verticies_adjacent[y_node].pop(z_loc_y)     #removes the edges
                    verticies_adjacent[z_node].pop(y_loc_z)

                    sayi = None
                    for x in graph_sub[i]:
                        if len(verticies_adjacent[x]) > 0:      #finds a vertex that has degree > 0 after the removal of edges.
                            sayi = x
                            break
                    if sayi == None:        #if it cannot find a vertex that has degree > 0 in the subgraph, adds the edge back
                        verticies_adjacent[y_node].append(z_node)
                        verticies_adjacent[z_node].append(y_node)
                        continue
                    elif DFS_connected(sayi,verticies_adjacent,[])[1]==graph_sub[i]:       #if graph still gives the same subgraph after removal of edge, this means this edge can be removed
                                                                                            #therefore we can connect this verticies to the other subgraphs.
                        changable.append([y_node,z_node])
                    else:
                        verticies_adjacent[y_node].append(z_node)       #if graph does not give the same subgraph, it means it is a bridge. We cannot take it.3
                        verticies_adjacent[z_node].append(y_node)
                        continue
    if changable == []:
        return 0 , None
    else:
        if len(changable) >= len(graph_sub)-1 : #checks if there are enough edges to be removed to connect the subgraphs.
            return 1 , len(changable)
        else:
            return 0 , 0

def sequential_algorithm(verticies_degree,adjacent,E,index=0,):
    m_degeri = sum(verticies_degree.values())
    while index < m_degeri:
        prob = list()

        c_degree = verticies_degree.copy()
        temp = list()
        for values in c_degree:
            if c_degree[values] == 0:
                temp.append(values)  # finding the nodes that achieved its limit
        for p in temp:
            c_degree.pop(p)  # removing the nodes that achieved its limit

        #check if all the verticies have reached its limit, if so stop and give the adjacent list.
        max_vertex = max(verticies_degree, key=verticies_degree.get)
        if verticies_degree[max_vertex] == 0:
            for a in range(len(E)):
                first_vertex = E[a][0]
                second_vertex = E[a][1]
                adjacent[first_vertex].append(second_vertex)
                adjacent[second_vertex].append(first_vertex)
            return adjacent
        vertex = min(c_degree, key=verticies_degree.get)

        #defining probabilities
        v = c_degree.copy()
        v.pop(vertex)    #choosen vertex will have 0 prob of to be choosen again

        #determining probabilities
        toplam=0
        for x in v:
            toplam += v[x]
        for i in v:
            prob.append(v[i] / toplam)


        deger = False
        visited_E = list()
        while deger==False:
            y = random.choices(list(v), weights=prob, k=1)
            if [vertex, y[0]] not in E:
                if y[0] not in visited_E:
                    E.append([vertex, y[0]])
                    for s in E[index]:
                        verticies_degree[s] -= 1
                        c_degree[s] -= 1
                    # checking if the remaning sequence is graphical
                    sequence = list(c_degree.values())
                    deger = check_if_graph(sequence)
                    if deger == False:
                        for t in E[index]:
                            verticies_degree[t] += 1
                            c_degree[t] += 1
                        E.pop(index)
                        visited_E.append(y[0])
                    else:
                        break
        #controlling if the remaning sequence is graphical. If not restart the process. If graphial continue with new edge.
        if deger:
            index+=1


def sequence_powerlaw(seq,n):
    alfa = -5/2
    prob = list()
    for i in range(1,n):
        prob.append(i**alfa)
    sequence = list()
    sequence = random.choices(population=range(1,n),weights=prob[:n],k=n)    #creates random sequence
    if sum(sequence) < 2*n-2:   #this is min requirement for a graph to be connected.
                                #for example [1,2,2,2,1] sum = 8, total node = 5
        seq = sequence_powerlaw(seq,n)
    elif check_if_graph(sequence):  #if the sequence is graphical returns the sequence
        return sequence
    else:       #if sequence is not graphical re-enters the function.
        seq = sequence_powerlaw(seq,n)
    return seq

def sequence_binomial(p,seq,n):
    prob = binom.pmf(range(1,n), n-1, p)
    sequence = list()
    sequence = random.choices(population=range(1,n),weights=prob[:n-1],k=n)  # creates random sequence
    if sum(sequence) < 2 * n - 2:  # this is min requirement for a graph to be connected.
        # for example [1,2,2,2,1] sum = 8, total node = 5
        seq = sequence_binomial(p,seq,n)
    elif check_if_graph(sequence):  # if the sequence is graphical returns the sequence
        return sequence
    else:  # if sequence is not graphical re-enters the function.
        seq = sequence_binomial(p,seq,n)
    return seq


seq_sequences = list()
seq_algorithm_time = list()
seq_connectivity = list()
seq_n_pairwise = list()
seq_adjacent_list = list()

HTH_sequences = list()
HTH_algorithm_time = list()
HTH_connectivity = list()
HTH_n_pairwise = list()
HTH_adjacent_list = list()

STH_sequences = list()
STH_algorithm_time = list()
STH_connectivity = list()
STH_n_pairwise = list()
STH_adjacent_list = list()

RTH_sequences = list()
RTH_algorithm_time = list()
RTH_connectivity = list()
RTH_n_pairwise = list()
RTH_adjacent_list = list()


for dene in range(36):
    sequence = list()
    print(dene)
    if dene < 3:
        sequence = sequence_powerlaw([],40)
    elif dene < 6:
        sequence = sequence_powerlaw([],100)
    elif dene < 9:
        sequence = sequence_powerlaw([], 250)
    elif dene < 12:
        sequence = sequence_binomial(0.2,[],40)
    elif dene < 15:
        sequence = sequence_binomial(0.2,[],100)
    elif dene < 18:
        sequence = sequence_binomial(0.2,[],250)
    elif dene < 21:
        sequence = sequence_binomial(0.5,[],40)
    elif dene < 24:
        sequence = sequence_binomial(0.5,[],100)
    elif dene < 27:
        sequence = sequence_binomial(0.5,[],250)
    elif dene < 30:
        sequence = sequence_binomial(0.7,[],40)
    elif dene < 33:
        sequence = sequence_binomial(0.7,[],100)
    elif dene < 36:
        sequence = sequence_binomial(0.7,[],250)

    seq_sequences.append(sequence)

    ###Sequential algorithm
    deg, adj = duzenleme(sequence)
    start = time.time()
    sonuc = sequential_algorithm(deg,adj,[])
    seq_adjacent_list.append(sonuc)
    if DFS_connected(0,sonuc,[])[0]:
        seq_connectivity.append(1)
        seq_n_pairwise.append(0)
    else:
        x , number = make_connected(deg,sonuc)
        seq_connectivity.append(0)
        seq_n_pairwise.append(number)
    end=time.time()
    seq_algorithm_time.append(end-start)


    lines = list()
    lines.append(str(seq_connectivity[dene]))
    lines.append(str(seq_n_pairwise[dene]))
    lines.append(str(seq_algorithm_time[dene]))
    for y in seq_adjacent_list[dene]:
        lines.append(' '.join([str(elem) for elem in seq_adjacent_list[dene][y]]))
    with open('Group16-' + str(len(seq_sequences[dene])) + '-' + str((sum(seq_sequences[dene]) // 2)) + '-Input-' + str(dene+1) + '-Output-' + str(dene+1) + '-Sequential.txt', 'w') as f:
        for line in lines:
            f.write(line)
            f.write('\n')

    input = list()
    input.append(' '.join([str(elem) for elem in seq_sequences[dene]]))
    with open('Group16-' + str(len(seq_sequences[dene])) + '-' + str(str((sum(seq_sequences[dene]) // 2))) + '-Input-' + str(dene+1) + '.txt',
              'w') as f:
        for line in input:
            f.write(line)
            f.write('\n')
            


    ###SMALL TO  HIGH
    deg, adj = duzenleme(sequence)
    start = time.time()
    STH_sonuc = HH_small_to_high(deg, adj)
    STH_adjacent_list.append(STH_sonuc)
    if DFS_connected(0, STH_sonuc, [])[0]:
        STH_connectivity.append(1)
        STH_n_pairwise.append(0)
    else:
        x, number = make_connected(deg, STH_sonuc)
        STH_connectivity.append(0)
        STH_n_pairwise.append(number)
    end = time.time()
    STH_algorithm_time.append(end-start)


    lines = list()
    lines.append(str(STH_connectivity[dene]))
    lines.append(str(STH_n_pairwise[dene]))
    lines.append(str(STH_algorithm_time[dene]))
    for y in STH_adjacent_list[dene]:
        lines.append(' '.join([str(elem) for elem in STH_adjacent_list[dene][y]]))
    with open('Group16-' + str(len(seq_sequences[dene])) + '-' + str((sum(seq_sequences[dene]) // 2)) + '-Input-' + str(dene + 1) + '-Output-' + str(dene + 1) + '-HH_small_to_high.txt', 'w') as f:
        for line in lines:
            f.write(line)
            f.write('\n')


    ###RANDOM TO HIGH
    deg, adj = duzenleme(sequence)
    start = time.time()
    RTH_sonuc = HH_random_to_high(deg, adj)
    RTH_adjacent_list.append(RTH_sonuc)
    if DFS_connected(0, RTH_sonuc, [])[0]:
        RTH_connectivity.append(1)
        RTH_n_pairwise.append(0)
    else:
        x, number = make_connected(deg, RTH_sonuc)
        RTH_connectivity.append(0)
        RTH_n_pairwise.append(number)
    end = time.time()
    RTH_algorithm_time.append(end-start)


    lines = list()
    lines.append(str(RTH_connectivity[dene]))
    lines.append(str(RTH_n_pairwise[dene]))
    lines.append(str(RTH_algorithm_time[dene]))
    for y in RTH_adjacent_list[dene]:
        lines.append(' '.join([str(elem) for elem in RTH_adjacent_list[dene][y]]))
    with open('Group16-' + str(len(seq_sequences[dene])) + '-' + str((sum(seq_sequences[dene]) // 2)) + '-Input-' + str(dene + 1) + '-Output-' + str(dene + 1) + '-HH_random_to_high.txt', 'w') as f:
        for line in lines:
            f.write(line)
            f.write('\n')

    
    ###HIGH TO HIGH
    deg, adj = duzenleme(sequence)
    start = time.time()
    HTH_sonuc = HH_high_to_high(deg, adj)
    HTH_adjacent_list.append(HTH_sonuc)
    if DFS_connected(0, HTH_sonuc, [])[0]:
        HTH_connectivity.append(1)
        HTH_n_pairwise.append(0)
    else:
        x, number = make_connected(deg, HTH_sonuc)
        HTH_connectivity.append(0)
        HTH_n_pairwise.append(number)
    end = time.time()
    HTH_algorithm_time.append(end-start)

/
    lines = list()
    lines.append(str(HTH_connectivity[dene]))
    lines.append(str(HTH_n_pairwise[dene]))
    lines.append(str(HTH_algorithm_time[dene]))
    for y in HTH_adjacent_list[dene]:
        lines.append(' '.join([str(elem) for elem in HTH_adjacent_list[dene][y]]))
    with open('Group16-' + str(len(seq_sequences[dene])) + '-' + str((sum(seq_sequences[dene]) // 2)) + '-Input-' + str(dene + 1) + '-Output-' + str(dene + 1) + 'HH_high_to_high.txt', 'w') as f:
        for line in lines:
            f.write(line)
            f.write('\n')


print('bitti')

print(seq_sequences)
print("--")
print(seq_algorithm_time)
print(seq_connectivity)
print(seq_n_pairwise)
print("--")
print(HTH_algorithm_time)
print(HTH_connectivity)
print(HTH_n_pairwise)
print("--")
print(RTH_algorithm_time)
print(RTH_connectivity)
print(RTH_n_pairwise)
print("--")
print(STH_algorithm_time)
print(STH_connectivity)
print(STH_n_pairwise)
print("--")




