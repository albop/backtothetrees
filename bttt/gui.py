from bttt.trees import *
from bqplot import pyplot as plt
from ipywidgets import *

from collections import OrderedDict
def get_coordinates(tree):
    if isinstance( tree, DeathTree):
        coords = OrderedDict()
        for n in tree.nodes:
            x = len(n)
            y = 0.0
            if sum(n)==1:
                x += 0.3
                y = 0.5
            elif sum(n)==2:
                x += 0.6
                y = 1.0
            coords[n] = (x,y)
    elif isinstance( tree, DeterministicTree):
        coords = OrderedDict()
        for i,n in enumerate(tree.nodes):
            coords[n] = (float(i),0.0)
    else:
        raise Exception("Unknown tree type")
    return coords

class VTree:
    def __init__(self, tree):
        self.tree = tree
        self.coords = get_coordinates(tree)
        self.create_figure()
        self.selected = []
    def create_figure(self):
        cc = self.coords
        tree = self.tree
        fig = plt.figure(figsize=(8,0.5))
        lines = OrderedDict()
        for n in tree.nodes:
            c0 = cc[n]
            for nn in tree.children(n):
                c1 = cc[nn]
                l = plt.plot([c0[0],c1[0]],[c0[1],c1[1]],'b')
                lines[(n,nn)] = l
        xv = [c[0] for c in (cc).values()]
        yv = [c[1] for c in (cc).values()]
        points = plt.plot(xv,yv,'o')
        self.fig = fig
        self.lines = lines
        self.points = points
    def deselect(self):
        N = len(self.points.colors)
        self.points.colors = ['steelblue']*N
        self.selected = []
    def select(self, n):
        if isinstance(n,int):
            i = n
            n = self.tree.nodes[i]
        else:
            i = self.tree.nodes.index(n)
        N = len(self.tree.nodes)
        cc =  ['steelblue']*N
        cc[i] = 'red'
        for s in self.tree.children(n):
            i = self.tree.nodes.index(s)
            cc[i] = 'green'
        par = self.tree.parent(n)
        if len(par)>0:
            i = self.tree.nodes.index(par)
            cc[i] = 'pink'
        self.points.colors = cc
        self.selected=[n]

class EquationsWidget:
    def __init__(self):
        self.widget = HBox()
        self.equations = ()
        self.conditions = []

    def set_equations(self, eqs, conditions=None):
        self.equations = tuple(eqs)
        if conditions is None:
            conditions = [("", True) for eq in self.equations]
        self.conditions = conditions
        cond_widgets = []
        eq_widgets = []
        for i,eq in enumerate(self.equations):
            cc,bo = self.conditions[i]
            cc = cc.replace("<", "&lt")
            cc = cc.replace(">", "&gt")
            teqw = Text(eq)
            if not bo:
                teqw.disabled=True
            else:
                teqw.disabled=False

            eq_widgets.append(teqw)
            if bo:
                w = HTML('<div style="background-color:#FFFF99">'+cc+"</div>")
            else:
                w = HTML(cc)
            cond_widgets.append(w)
        self.widget.children = (VBox(cond_widgets), VBox(eq_widgets))


class ModelGUI:
    def __init__(self, model):
        self.model = model
        self.vtree = VTree(model.etree)
        self.eqw = EquationsWidget()
        self.widget = HBox([self.vtree.fig, self.eqw.widget])
        self.select(0)
        self.make_connections()

    # def get_equations(self):
    #     equations = []
    #     conditions = []
    #     import re
    #     regxeq = re.compile("((.*):|)(.*)⟂(.*)")
    #     for l in self.model.source['equations']:
    #         m = regxeq.match(l)
    #         cond, eq, comp = m.groups()[1:]
    #         equations.append(eq)
    #         conditions.append(cond)
    #     return (equations, conditions)
    #

    def select(self,i):
        self.vtree.select(i)
        n = self.model.etree.nodes[i]

        equations = self.model.equations
        conditions_ = self.model.conditions
        conditions = []
        etree = self.model.etree
        for c in conditions_:
            context = {"s":n, "t":len(n)-1,"T": len(etree.nodes)-1 }
            if c is None:
                c = "_"
                cond = True
            else:
                cond = eval(c, context)
            conditions.append((c, cond))
        equations = [str(e) for e in self.model.equations]
        self.eqw.set_equations(equations, conditions)

    def make_connections(self):
        self.vtree.points.on_element_click(lambda a,b: self.select(b["data"]["index"]))
