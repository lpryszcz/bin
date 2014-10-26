"""
Some graph functions. Usefull, to work with orthologs.
"""
import sys
import numpy as np

class HomologsGraph(object):
    """Undirected Graph class for homologs representation.
    Note: score (sos) is given as int. maybe change to float,
    but handling floats consumes 25% more memory!
    And when compares to binary score handling, it's almost 200%.
    """
    def __init__(self,verbose=0):
        """Construct an empty graph."""
        #store taxa and protein info
        self.taxa      = set()
        self.prot2taxa = {}
        #store relationships and scores
        self.homologs  = {}
        self.scores    = {}
        self.verbose   = verbose
        #count homologies
        self.homologiesCount  = 0
        
    def __len__(self):
        return len(self.taxa)
        
    def __str__(self):
        """Produce string representation of the graph"""
        out  = '%s proteins from %s taxa %s\n' % (len(self.prot2taxa),len(self.taxa),str(self.taxa))  
        out += '%s homologs derived from %s homologies\n' % (len(self.scores),self.homologiesCount)
        out += 'homologs pair\tscore\tcount\n'
        i = 0
        for pair,intscore in self.scores.items():
            v1,v2 = pair.split('-')
            if i > 10:
                break
            i += 1
            pair = "%s-%s" % (v1,v2)
            #unload scores: convert int (10) to binary 0b1010 and to scores (0,1,0)
            scores = [ int(x) for x in bin(intscore)[3:] ]
            #scores = self.scores[pair]
            out += '%s - %s\t%.3f\t%s\n' % (v1,v2,np.mean(scores),len(scores))
        return out

    def _update_homologs(self,v):
        """Update homologs"""
        if v not in self.homologs:
            self.homologs[v] = set()        

    def _update_connections(self,scores):
        """Update scores."""
        for pair,intscore in self.scores.items():
            v1,v2 = pair.split('-')
            #add v1 and v2 if not in neighbours
            self._update_homologs(v1)
            self._update_homologs(v2)
            #store score
            if pair not in self.scores:
                #add connection if not yet stored
                self.homologs[v1].add(v2)
                self.homologs[v2].add(v1)
                self.scores[pair] = intscore
            else:
                #or merge current with previous scores
                self.scores[pair] = int( bin(self.scores[pair])+bin(iscore)[3:],2 ) 
            
    def __radd__(self,anotherHomologsGraph):
        self.__add__(anotherHomologsGraph)
                
    def __add__(self,anotherHomologsGraph):
        """Merge two homologs graphs.
        No error checking."""
        #add taxa
        self.taxa.update(anotherHomologsGraph.taxa) 
        #add prot2taxa
        self.prot2taxa.update(anotherHomologsGraph.taxa)
        #add connections
        self._update_connections(anotherHomologsGraph.scores)
        #update homologies count
        self.homologiesCount += anotherHomologsGraph.homologiesCount

    def _add_protein(self,v,taxid):
        """Add taxid and protein to taxa"""
        self.taxa.add(taxid)
        self.prot2taxa[v] = taxid

    def _add_homologs(self,v1,v2,sos):
        """Add pair to neighbours and their connection score (INT!)"""
        #define pair
        if v1 < v2:
            pair = "%s-%s" % (v1,v2)
        else:
            pair = "%s-%s" % (v2,v1)
        #add v1 and v2 if not in neighbours
        self._update_homologs(v1)
        self._update_homologs(v2)
        #check if connection exists
        if pair not in self.scores:
            #add connection if not yet stored
            self.homologs[v1].add(v2)
            self.homologs[v2].add(v1)
            self.scores[pair] = int('0b1',2) 
        #store score
        self.scores[pair] = int( bin(self.scores[pair]) + str(sos),2 ) #.append(sos)
        #count
        self.homologiesCount += 1

    def add_homologs(self,taxid1,v1,taxid2,v2,sos):
        """Add connection between pair.
        Update taxa, homologs and scores
        """
        sos = int(sos)
        #update taxa info
        self._add_protein(v1,taxid1)
        self._add_protein(v2,taxid2)
        #add connection
        self._add_homologs(v1,v2,sos)
        
    def get_orthologs(self,CSth=0.5,out=sys.stdout):
        """Report orthologs to output stream.
        By default to stdout."""
        sys.stderr.write( "Reporting orthologs to %s\n" % out.name )
        io = 0
        ih = len(self.scores)
        for pair in self.scores:
            v1,v2 = pair.split('-')
            #unload scores
            scores = [ int(x) for x in bin(self.scores[pair])[3:] ]
            #homologies store 0 for orthologs and 1 for paralogs
            if np.mean(scores) > CSth:
                continue
            #store homologs
            io += 1
            out.write( "%s\t%s\n" % (v1,v2) )

        sys.stderr.write( " %s homologs parsed. %s orthologs having CS >= %s reported.\n" % (ih,io,CSth) )

    def get_all_homologs(self,v,processed=set()):
        """Return homologs (and homologs of homologs) for v.
        UNFINISHED!"""
        selected = set()
        processed.add(v)
        for v2 in self.homologs[v]:
            selected += self.homologs[v2]

class DirectedGraph(object):
    """Directed graph. Keeps track of direction and number of
    connections between elements.
    Apply for contig joining.
    """
    def __init__(self,vertices):
        """Construct a graph with the given vertices, but no lines
        """
        self._neighbours  = {}
        self._maxConnects = 1
        for vertex in vertices:
            self._neighbours[vertex] = {}

    def __str__(self):
        """Produce string representation of the graph
        """
        out = 'Vertices: %s \nConnections:\n' % self._neighbours.keys() 
        for vertex in self._neighbours:
            for neighbour,strength in self._neighbours[vertex].iteritems():
                out += '\t%s -> %s (%s)' % (vertex,neighbour,strength)
        return out
    
    def add_connection(self,v1,v2):
        """Add directed connection from v1 to v2 (with no error checking!)
        If such connection already present, increase it's strengh.
        """
        if v2 not in self._neighbours[v1]:
            self._neighbours[v1][v2] = 1
        else:
            self._neighbours[v1][v2]+= 1

    def _check(self,vertices):
        if not frozenset(vertices).issubset(self._neighbours):
            raise ValueError('%s contains non-vertices' % vertices)    
        
class MyGraph(object):
    """Undirected Graph class.
    """
    ### from Grap.py (J.Cussens, York, 2008)
    def __init__(self, vertices=[]):
        """Construct a graph with the given vertices, but no lines
        """
        self._neighbours = {}
        self._features = {}
        for vertex in vertices:
            self._neighbours[vertex] = set()
            self._features[vertex] = {}

    def add_vertice(self, vertex, key="", value=""):
        """Add vertice to vertices.
        If already present, just update vertice features."""
        if vertex not in self._neighbours:
            self._neighbours[vertex] = set()
            self._features[vertex] = {}
        if key:
            self._features[vertex][key] = value

    def add_vertices(self, vertices, key="", value=""):
        """Add vertice to vertices.
        If already present, just update vertice features."""
        for vertex in vertices:
            if vertex not in self._neighbours:
                self._neighbours[vertex] = set()
                self._features[vertex] = {}
            if key:
                self._features[vertex][key] = value

    def add_line(self, v1, v2):
        """Add a line from v1 to v2 (with no error checking!)
        """
        self._neighbours[v1].add(v2)
        self._neighbours[v2].add(v1)

    def __str__(self):
        """Produce string representation of the graph
        """
        out = 'Vertices: %s \nLines:\n' % self._neighbours.keys() 
        for vertex, neighbours in self._neighbours.items():
            for neighbour in neighbours:
                if vertex < neighbour:
                    out += '\t%s - %s\n' % (vertex,neighbour)
        return out

    #get interconnected
    def get_clusters(self, key="", value=""):
        """Return connected nodes having given feature (key: value)."""
        processed = set()
        clusters = []
        for v in self._neighbours:
            if v not in processed:
                clusters.append([v])
                processed.add(v)
            i = 0
            while i<len(clusters[-1]): #for v2 in self._neighbours[v]:
                v = clusters[-1][i]
                for v2 in self._neighbours[v]:
                    if v2 in processed:
                        continue
                    clusters[-1].append(v2)
                    processed.add(v2)
                i += 1
        return clusters
        
    ### end of Grap.py (J.Cussens, York, 2008)  
    ### by LPryszcz, York, 2008
        
    def _check(self, vertices):
        """Check if vertice in vertices.
        """
        if not frozenset(vertices).issubset(self._neighbours):
            raise ValueError('%s contains non-vertices' % vertices)

    def delete_line(self,v1,v2):
        """Delete line between v1 and v2 if any. 
        Return error if no line between them.
        """
        self._check((v1,v2))
        try:
            self._neighbours[v1].remove(v2)
            self._neighbours[v2].remove(v1)
        except KeyError:
            raise ValueError('No line between %s and %s' % (v1,v2))

    def independent(self,vertices):
        """Return True only if no two vertices in vertices 
        are connected. Otherwise return False.
        """
        self._check(vertices)
        vertices = tuple(vertices)
        for i,v in enumerate(vertices):
            v_neighbours = self._neighbours[v]
            for w in vertices[i+1:]:
                if w in v_neighbours:
                    return False
        return True

    def complete(self,vertices):
        """Return True only if all vertices in vertices
        are connected. Otherwise return False.
        If vertices contain element that is not a vertex 
        in graph self, return err message. 
        """
        self._check(vertices)
        vertices = tuple(vertices)
        for i, v in enumerate(vertices):
            vnbrs = self._neighbours[v]
            for w in vertices[i+1:]:
                if w not in vnbrs: 
                    return False
        return True
  
    def complete_with_inparalogs(self,vertices,id2group):
        """Return True only if all vertices in vertices
        not belonging to same group accordinly to id2group
        are connected with line. Otherwise return False.
        If vertices contain element that is not a vertex 
        in graph self, return err message. 
        """
        self._check(vertices)
        vertices = tuple(vertices)
        for i, v in enumerate(vertices):
            vnbrs = self._neighbours[v]
            for w in vertices[i+1:]:
                if w not in vnbrs: 
                    if id2group(v)!=id2group(w): 
                        return False
        return True

    def clique(self,vertices):
        """Return true only if all vertices in verices 
        are connected and there is no vertex 
        in a graph self which is not in vertices, but to 
        which all vertices from vertices are connected. 
        Return error if vertices contain element that is not 
        a vertex in a graph self.
        """
        self._check(vertices)
        vertices = tuple(vertices)
        if not self.complete(vertices):
            return False
        for other in self._neighbours[vertices[0]] - frozenset(vertices):
            othernbrs = self._neighbours[other]
            for v in vertices:
                if v not in othernbrs:
                    break
                else:
                    return False
        return True
  
    def unique(self,vertices):
        """Return true only if all vertices from vertices 
        are connected by a line, and none of vertices is 
        connected with vertices from outside vertices 
        to vertice which name starts with the same 3-letter 
        code that is arleady present in vertices.
        Return error if vertices contain element that is not 
        a vertex in a graph self.
        """
        self._check(vertices)
        vertices = tuple(vertices)
        if not self.complete(vertices):
            return False
        vspecies = map(lambda x: x[:3], vertices)
        for v in vertices:
            for w in self._neighbours[v]:
                if w not in vertices and w[:3] in vspecies:
                    return False
        return True

    def _get_spcode( l ):
        return l.split('_')[-1]

    def nfamily_one2one(self,species,n=2,_get_spcode=_get_spcode,verbose=False):
        """Return groups of one2one orthologs only if:
        - number of vertices for every specie from species is equal to n
        - all vertices from each specie from species are connected to one
        and only one vertice from other species from species.
        Species defining function can be provided. By default, it's last
        part of vertice name ie Phy0008594_HUMAN -> HUMAN.
        ie. 1_A - 1_B; 2_A - 2_B
        """
        #check if one2one connections between all elements from all species
        groups = []
        for i in range( len(species)-1 ):
            sp1 = species[i]
            vertices1 = filter( lambda x: _get_spcode(x)==sp1,self._neighbours.keys() )      
            #check if correct number of vertices for every specie
            if len(vertices1)!=n:        
                return False
            #iterate through other species
            for j in range( i+1,len(species) ):
                #copy vertices1
                #_vertices = set( vertices1 )
                sp2 = species[j]
                vertices2 = filter( lambda x: _get_spcode(x)==sp2,self._neighbours.keys() )
                if verbose:
                    print vertices1, vertices2
                #check if correct number of vertices for every specie      
                if len(vertices2)!=n:
                    return False
                #check connections
                k=0
                for v1 in vertices1:
                    if len(groups)<=k:
                        groups.append( [v1] )
                    #get only connected neigbours
                    neighbours = filter( lambda x: _get_spcode(x)==sp2,self._neighbours[v1] )
                    if verbose:
                        print v1,"-",",".join( neighbours )#self._neighbours[v1],neighbours
                    #check if connected to only one and vertices2 wasn't previously connected to this
                    if len(neighbours)!=1 or neighbours[0] not in vertices2:
                        return False
                    else:
                        vertices2.remove( neighbours[0] )
                        if i==0:
                            groups[k].append( neighbours[0] )
                    k+=1
        
        return groups
    