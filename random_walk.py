import numpy
from scipy import linalg

# Adapted from the GUILD software
# Guney E, Oliva B (2012) Exploiting Protein-Protein Interaction Networks for Genome-Wide Disease-Gene Prioritization. PLoS ONE 7(9): e43557. doi:10.1371/journal.pone.0043557

restart = 0.75
max_nb_iter = 50
convergence_cutoff = 1e-6

def get_adjacency_matrix(edge_file, nodeList):

    # Creation of adjacency matrix
    nbNodes = len(nodeList)
    adjacency = numpy.tile(1e-9, [nbNodes, nbNodes])

    # File reading and adjacency matrix modification
    with open(edge_file, 'r', encoding = 'utf-8') as input:
        for line in input:
            lineList = line.split(' ')
            id1Index = nodeList.index(lineList[0].strip())
            id2Index = nodeList.index(lineList[2].strip())
            adjacency[id1Index, id2Index] = lineList[1]
            adjacency[id2Index, id1Index] = lineList[1]

    return(adjacency)

def get_node_score(node_file):
    
    # Creation of an empty set of nodes
    nodeDic = {}
    
    #File reading and nodeList creation
    with open(node_file, 'r', encoding = 'utf-8') as input:
        for line in input:
            lineList = line.split(' ')
            nodeDic[lineList[0].strip()] = lineList[1]
 
    # Elimination of redondancy
    nodeList = list(nodeDic.keys())
    nodeList.sort()
    nbNodes = len(nodeList)
    
    # Creation of the scores vector
    scores = numpy.zeros((nbNodes,1))
    
    # File reading and score vector modification
    for node in nodeList:
        idIndex = nodeList.index(node)
        scores[idIndex, 0] = nodeDic[node]
            
    return(scores,nodeList)

def random_walk_with_restart(adjacency, p_0, restart, max_nb_iter, convergence_cutoff):

    # matrice d'adjacence normalisée par colonne
    for row in range(0, adjacency.shape[1]):
        divisor = sum(adjacency[:,row])
        for line in range(0, adjacency.shape[0]):
            adjacency[line,row] = adjacency[line,row]/divisor
    
    # Assignation des probabilités
    p_0 = p_0/sum(p_0)

    # Initialisation du vecteur p_t
    p_t = p_0
    
    # Itération pour calculer le vecteur de probabilité
    for iteration in range(1,max_nb_iter+1):
        p_tx = (1-restart)*numpy.dot(adjacency, p_t) + restart*p_0
        if linalg.norm(p_tx-p_t) < convergence_cutoff:
            break
        p_t = p_tx
   
    return(p_tx)

def run_random_walk_with_restart(node_file, edge_file, out_file):
    (p_0, nodeList) = get_node_score(node_file)
    adjacency = get_adjacency_matrix(edge_file, nodeList)
    p_tx = random_walk_with_restart(adjacency, p_0, restart, max_nb_iter, convergence_cutoff)
    with open(out_file, 'w', encoding = 'utf-8') as output:
        for line in range(0, p_tx.shape[0]):
            output.write(nodeList[line]+' '+str(p_tx[line,0])+'\n')