
def contained(t1,t2):
    if len(t1)>len(t2):
        return False
    elif len(t1)==len(t2):
        return t1==t2
    else:
        return t1==t2[:len(t1)]

class EventTree:

    nodes = None # list of nodes
    probas = dict()
    values = dict()

    def children(self, s):
        # returns a list of all children of s
        if len(s)==len(self):
            return [s]
        else:
            return [e for e in self.nodes if (s,e) in self.probas]

    def parent(self, s):
        # returns parent of s
        return s[:-1]

    def history(self, s):
        return [s[:i] for i in range(1,len(s)+1)]

    def __len__(self):
        return max([len(e) for e in self.nodes])

    def graph(self, plot=False):
        import networkx as nx
        G = nx.DiGraph()
        for n1,n2 in self.probas:
            G.add_edge(n1,n2,{'p':self.probas[n1,n2]})
        if plot:
            nx.draw_spring(G)
        else:
            return G



class DeterministicTree(EventTree):
    def __init__(self, N):

        etree = self
        etree.nodes = [(0,)*t for t in range(1,N+1)]
        for s in etree.nodes:
            etree.values[s] = 0
            if len(s)>1:
                etree.probas[(etree.parent(s), s)] = 1
            if len(s)==len(etree):
                etree.probas[(s,s)] = 1


class BranchTree(EventTree):

    def __init__(self, N, T, p):
        '''
        p: [p_0, p_1]
        '''

        T = T+1
        etree = self
        nodes = [(0,)*t for t in range(1,T)]
        new_paths = []
        for i in range(2):
            nn = nodes[-1] + (i,)
            new_paths.append([nn+(0,)*t for t in range(0,N-T+1)])
        for i in range(len(new_paths[0])):
            for b in range(2):
                nodes.append(new_paths[b][i])
        etree.nodes = nodes

        for s in etree.nodes:
            etree.values[s] = 0
            if len(s)<len(etree):
                if len(etree.children(s))>1:
                    etree.probas[(s, etree.children(s)[0])] = p[0]
                    etree.probas[(s, etree.children(s)[1])] = p[1]
                else:
                    etree.probas[(s, etree.children(s)[0])] = 1
            if len(s)==len(etree):
                etree.probas[(s,s)] = 1


class DeathTree(EventTree):

    def __init__(self, N, p):
        '''
        p: [p_0, p_1]
        '''

        self.nodes = [(0,)]
        self.probas = dict()
        for t in range(1, N-1):
            nn = self.nodes[-1]
            nn0 = nn + (1,)
            nn1 = nn + (1, 1)
            self.probas[(nn, nn0)] = p
            self.probas[(nn0, nn1)] = 1.0
            self.probas[(nn1, nn1)] = 1.0
            self.nodes.extend([nn0, nn1])
            if t < N-2:
                nn2 = nn + (0,)
                self.nodes.append(nn2)
                self.probas[(nn, nn2)] = 1-p
            else:
                self.probas[(nn, nn)] = 1-p


class PersistentDeathTree(EventTree):

    def __init__(self, N, p, k=5):
        '''
        p: [p_0, p_1]
        '''

        self.nodes = [(0,)]
        self.probas = dict()
        for t in range(1, N-1):
            nn = self.nodes[-1]
            nn0 = nn + (1,)
            self.probas[(nn, nn0)] = p
            for j in range(k):
                nn1 = nn0 + (1,)
            # nn0 = nn +
            nn1 = nn + (1, 1)
            self.probas[(nn0, nn1)] = 1.0
            self.probas[(nn1, nn1)] = 1.0
            self.nodes.extend([nn0, nn1])
            if t < N-2:
                nn2 = nn + (0,)
                self.nodes.append(nn2)
                self.probas[(nn, nn2)] = 1-p
            else:
                self.probas[(nn, nn)] = 1-p


        # nodes = [(0,)*t for t in range(1,T)]
        # new_paths = []
        # for i in range(2):
        #     nn = nodes[-1] + (i,)
        #     new_paths.append([nn+(0,)*t for t in range(0,N-T+1)])
        # for i in range(len(new_paths[0])):
        #     for b in range(2):
        #         nodes.append(new_paths[b][i])
        # etree.nodes = nodes
        #
        # for s in etree.nodes:
        #     etree.values[s] = 0
        #     if len(s)<len(etree):
        #         if len(etree.children(s))>1:
        #             etree.probas[(s, etree.children(s)[0])] = p[0]
        #             etree.probas[(s, etree.children(s)[1])] = p[1]
        #         else:
        #             etree.probas[(s, etree.children(s)[0])] = 1
        #     if len(s)==len(etree):
        #         etree.probas[(s,s)] = 1


def get_ts(etree, sol, vname, terminal_state=None, ts_ind=0):
    import numpy
    if terminal_state is None:
        terminal_states = [e for e in etree.nodes if len(e)==len(etree)]
        terminal_state = terminal_states[ts_ind]
    history = etree.history(terminal_state)
    his_inds = [etree.nodes.index(e) for e in history]
    vals = [sol["{}_{}".format(vname, h)] for h in his_inds]
    return numpy.array(vals).astype(dtype=float)
