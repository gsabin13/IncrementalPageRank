import networkx as nx

G = nx.DiGraph()

f = open("edgelist.in")

first_line = f.readline();
num_vertices, num_edges = [int(x) for x in first_line.split()]

print "num_vertices, num_edges", num_vertices, num_edges

for line in f:
	a, b = [int(x) for x in line.split()]
	G.add_edge(a,b)

f.close()

li = G.adjacency_list()

print li

#G=nx.gnp_random_graph(50,0.01,directed=True)

MAX_ITER =300

y = nx.pagerank(G,tol=1e-10,max_iter=MAX_ITER)

print ""

print y

y = nx.pagerank_scipy(G,tol=1e-10,max_iter=MAX_ITER)

print ""
print y

y = nx.pagerank(G,max_iter=MAX_ITER)
print ""
print y